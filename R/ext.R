#' Calculate location and scale parameters for a specified distribution so that
#' it matches two specified quantiles
#'
#' @param ps vector of two probability levels at which the distribution's
#'   quantiles are distinct
#' @param qs vector of two distinct quantile values corresponding to the
#'   probability levels in ps
#' @param dist the probability distribution to use for extrapolation. This
#'   distribution should be in a location-scale family, such as `"norm"`al or
#'   `"cauchy"`
#' @param dup_tol numeric tolerance for identifying duplicated values indicating
#'   a discrete component of the distribution. If there is a run of values where
#'   each consecutive pair is closer together than the tolerance, all are
#'   labeled as duplicates even if not all values in the run are within the
#'   tolerance.
#' @param zero_tol numeric tolerance for identifying values in `qs` that are
#'   (approximately) zero.
#' @param zero_discrete boolean indicating whether qs of zero should always
#'   indicate the presence of a point mass at 0.
#' @param is_lower boolean indicator of whether we are extending into the lower
#'   or upper tail
#'
#' @param named list with entries `"a"`, the location parameter, and `"b"`, the
#'   scale parameter
calc_loc_scale_params <- function(ps, qs, dist, dup_tol, zero_tol,
                                  zero_discrete, is_lower) {
    if (length(ps) != 2L || length(qs) != 2L) {
        stop("`ps` and `qs` must be of length 2.")
    }
    if (dist == "lnorm" && any(qs < 0.0)) {
        stop("For dist = 'lnorm', all qs must be non-negative")
    }

    is_zero <- (abs(qs) < zero_tol)
    if (zero_discrete && any(is_zero)) {
        if (all(is_zero) || (is_zero[1] && is_lower) ||
                (is_zero[2] && !is_lower)) {
            # point mass at 0. We detect this when any of:
            #  * both qs are equal to 0 (either tail)
            #  * first q is 0 and we're extending into the lower tail
            #  * second q is 0 and we're extending into the upper tail
            return(list(a = ifelse(dist == "lnorm", -Inf, 0.0), b = 0.0))
        } else {
            # point mass at non-zero value. We detect this when:
            #  * first q != 0, second q == 0, lower tail (mass at first q)
            #  * first q == 0, second q != 0, upper tail (mass at second q)
            # note: point mass at non-zero value detected from duplicated
            # non-zero qs is handled below by calculating b = 0
            result <- list(a = qs[2 - is_lower], b = 0.0)
            if (dist == "lnorm") {
                result$a <- log(result$a)
            }
            return(result)
        }
    }

    if (dist == "lnorm") {
        qs <- log(qs)
        qdst <- qnorm
    } else {
        qdst <- get(paste0("q", dist))
    }

    if (any(duplicated_tol(qs, dup_tol))) {
        b <- 0.0
    } else {
        b <- (qs[2] - qs[1]) / (qdst(ps[2]) - qdst(ps[1]))
    }
    a <- qs[1] - b * qdst(ps[1])

    return(list(a = a, b = b))
}


#' Extrapolate density function in a location-scale family matching specified
#' quantiles.
#'
#' @param ps vector of two probability levels at which the distribution's
#'   quantiles are distinct
#' @param qs vector of two distinct quantile values corresponding to the
#'   probability levels in ps
#' @param dist the probability distribution to use for extrapolation. This
#'   distribution should be in a location-scale family, such as "normal" or
#'   "Cauchy"
#'
#' @return a function with parameters `x` and `log` that can be used to
#'   evaluate the density function (or its log) of the distribution in the
#'   specified location-scale family that has quantiles matching those in `ps`
#'   and `qs`
d_ext_factory <- function(ps, qs, dist, dup_tol, zero_tol, zero_discrete,
                          is_lower) {
    c(a, b) %<-% calc_loc_scale_params(ps, qs, dist, dup_tol, zero_tol,
                                       zero_discrete, is_lower)
    if (b == 0) {
        stop("Detected a point mass; cannot create density function.")
    }

    if (dist == "lnorm") {
        d_ext <- function(x, log = FALSE) {
            return(dlnorm(x, meanlog = a, sdlog = b, log = log))
        }
    } else {
        ddst <- get(paste0("d", dist))
        d_ext <- function(x, log = FALSE) {
            result <- ddst((x - a) / b, log = TRUE) - log(b)
            if (log) {
                return(result)
            } else {
                return(exp(result))
            }
        }
    }

    return(d_ext)
}


#' Extrapolate cumulative distribution function in a location-scale family
#' matching specified quantiles.
#'
#' @param ps vector of two probability levels at which the distribution's
#'   quantiles are distinct
#' @param qs vector of two distinct quantile values corresponding to the
#'   probability levels in ps
#' @param dist the probability distribution to use for extrapolation. This
#'   distribution should be in a location-scale family, such as "normal" or
#'   "Cauchy"
#'
#' @return a function with parameter `x` and `log.p` that can be used to
#'   evaluate the cumulative distribution function (or its log) of the
#'   distribution in the specified location-scale family that has quantiles
#'   matching those in `ps` and `qs`
p_ext_factory <- function(ps, qs, dist, dup_tol, zero_tol, zero_discrete,
                          is_lower) {
    c(a, b) %<-% calc_loc_scale_params(ps, qs, dist, dup_tol, zero_tol,
                                       zero_discrete, is_lower)

    if (dist == "lnorm") {
        if (b == 0) {
            p_ext <- function(q, log.p = FALSE) {
                result <- as.numeric(q >= exp(a))
                if (log.p) {
                    return(log(result))
                } else {
                    return(result)
                }
            }
        } else {
            p_ext <- function(q, log.p = FALSE) {
                return(plnorm(q, meanlog = a, sdlog = b, log.p = log.p))
            }
        }
    } else {
        if (b == 0) {
            p_ext <- function(q, log.p = FALSE) {
                result <- as.numeric(q >= a)
                if (log.p) {
                    return(log(result))
                } else {
                    return(result)
                }
            }
        } else {
            pdst <- get(paste0("p", dist))

            p_ext <- function(q, log.p = FALSE) {
                return(pdst((q - a) / b, log.p = log.p))
            }
        }
    }

    return(p_ext)
}


#' Extrapolate quantile function in a location-scale family matching specified
#' quantiles.
#'
#' @param ps vector of two probability levels at which the distribution's
#'   quantiles are distinct
#' @param qs vector of two distinct quantile values corresponding to the
#'   probability levels in ps
#' @param dist the probability distribution to use for extrapolation. This
#'   distribution should be in a location-scale family, such as "normal" or
#'   "Cauchy"
#'
#' @return a function with parameter `p` that can be used to evaluate the
#'   quantile function of the distribution in the specified location-scale
#'   family that has quantiles matching those in `ps` and `qs`
q_ext_factory <- function(ps, qs, dist, dup_tol, zero_tol, zero_discrete,
                          is_lower) {
    c(a, b) %<-% calc_loc_scale_params(ps, qs, dist, dup_tol, zero_tol,
                                       zero_discrete, is_lower)

    if (dist == "lnorm") {
        if (b == 0) {
            q_ext <- function(p) {
                return(rep(exp(a), length(p)))
            }
        } else {
            q_ext <- function(p) {
                return(qlnorm(p, meanlog = a, sdlog = b))
            }
        }
    } else {
        qdst <- get(paste0("q", dist))

        if (b == 0) {
            q_ext <- function(p) {
                return(rep(a, length(p)))
            }
        } else {
            q_ext <- function(p) {
                return(a + b * qdst(p))
            }
        }
    }

    return(q_ext)
}


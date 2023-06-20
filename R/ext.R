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
#' @param is_lower boolean indicator of whether we are extending into the lower
#'   or upper tail
#' @param lnorm_zero_buffer boolean or numeric specifying how to handle zero
#'   quantiles in the upper tail when `dist = "lnorm"`. If `FALSE`, an error is
#'   raised if the first element of `qs` is zero, the second element of `qs` is
#'   non-zero, and `is_lower` is `FALSE`. Otherwise, must be a positive numeric
#'   value, and in this situation the first element of `qs` is replaced by
#'   `min(lnorm_zero_buffer, qs[2] / 2)`.
#'
#' @param named list with entries `"a"`, the location parameter, and `"b"`, the
#'   scale parameter
calc_loc_scale_params <- function(ps, qs, dist, is_lower, lnorm_zero_buffer) {
    if (dist == "lnorm") {
        if (any(qs < 0.0)) {
            stop("For dist = 'lnorm', all qs must be non-negative")
        } else if (any(qs == 0)) {
            if (all(qs == 0.0) || (qs[1] == 0.0 && is_lower)) {
                # point mass at 0. We detect this when either:
                #  * both qs are equal to 0 (either tail)
                #  * first q is 0 and we're extending into the lower tail
                return(list(a = 0.0, b = 0.0))
            } else {
                # (qs[1] == 0.0 && !is_lower)
                # first q is 0 and we're extending into the upper tail
                # don't have enough information to identify the location-scale
                # parameters
                qs <- apply_lnorm_zero_buffer(qs, lnorm_zero_buffer)
            }
        }
        qs <- log(qs)
        qdst <- qnorm
    } else {
        qdst <- get(paste0("q", dist))
    }
    b <- (qs[2] - qs[1]) / (qdst(ps[2]) - qdst(ps[1]))
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
#' @param is_lower boolean indicator of whether we are extending into the lower
#'   or upper tail
#' @param lnorm_zero_buffer boolean or numeric specifying how to handle zero
#'   quantiles in the upper tail when `dist = "lnorm"`. If `FALSE`, an error is
#'   raised if the first element of `qs` is zero, the second element of `qs` is
#'   non-zero, and `is_lower` is `FALSE`. Otherwise, must be a positive numeric
#'   value, and in this situation the first element of `qs` is replaced by
#'   `min(lnorm_zero_buffer, qs[2] / 2)`.
#'
#' @return a function with parameters `x` and `log` that can be used to
#'   evaluate the density function (or its log) of the distribution in the
#'   specified location-scale family that has quantiles matching those in `ps`
#'   and `qs`
d_ext_factory <- function(ps, qs, dist, is_lower, lnorm_zero_buffer) {
    c(a, b) %<-% calc_loc_scale_params(ps, qs, dist, is_lower,
                                       lnorm_zero_buffer)

    if (dist == "lnorm") {
        if (b == 0) {
            stop("Detected a point mass")
        }
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
#' @param is_lower boolean indicator of whether we are extending into the lower
#'   or upper tail
#' @param lnorm_zero_buffer boolean or numeric specifying how to handle zero
#'   quantiles in the upper tail when `dist = "lnorm"`. If `FALSE`, an error is
#'   raised if the first element of `qs` is zero, the second element of `qs` is
#'   non-zero, and `is_lower` is `FALSE`. Otherwise, must be a positive numeric
#'   value, and in this situation the first element of `qs` is replaced by
#'   `min(lnorm_zero_buffer, qs[2] / 2)`.
#'
#' @return a function with parameter `x` and `log.p` that can be used to
#'   evaluate the cumulative distribution function (or its log) of the
#'   distribution in the specified location-scale family that has quantiles
#'   matching those in `ps` and `qs`
p_ext_factory <- function(ps, qs, dist, is_lower, lnorm_zero_buffer) {
    c(a, b) %<-% calc_loc_scale_params(ps, qs, dist, is_lower,
                                       lnorm_zero_buffer)

    if (dist == "lnorm") {
        p_ext <- function(q, log.p = FALSE) {
            return(plnorm(q, meanlog = a, sdlog = b, log.p = log.p))
        }
    } else {
        pdst <- get(paste0("p", dist))

        p_ext <- function(q, log.p = FALSE) {
            return(pdst((q - a) / b, log.p = log.p))
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
#' @param is_lower boolean indicator of whether we are extending into the lower
#'   or upper tail
#' @param lnorm_zero_buffer boolean or numeric specifying how to handle zero
#'   quantiles in the upper tail when `dist = "lnorm"`. If `FALSE`, an error is
#'   raised if the first element of `qs` is zero, the second element of `qs` is
#'   non-zero, and `is_lower` is `FALSE`. Otherwise, must be a positive numeric
#'   value, and in this situation the first element of `qs` is replaced by
#'   `min(lnorm_zero_buffer, qs[2] / 2)`.
#'
#' @return a function with parameter `p` that can be used to evaluate the
#'   quantile function of the distribution in the specified location-scale
#'   family that has quantiles matching those in `ps` and `qs`
q_ext_factory <- function(ps, qs, dist, is_lower, lnorm_zero_buffer) {
    c(a, b) %<-% calc_loc_scale_params(ps, qs, dist, is_lower,
                                       lnorm_zero_buffer)

    if (dist == "lnorm") {
        if (b == 0) {
            q_ext <- function(p) {
                rep(a, length(p))
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
                rep(a, length(p))
            }
        } else {
            q_ext <- function(p) {
                return(a + b * qdst(p))
            }
        }
    }

    return(q_ext)
}

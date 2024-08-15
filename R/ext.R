a <- b <- NULL

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
#'
#' @importFrom stats qnorm
#'
#' @return named list with entries `"a"`, the location parameter, and `"b"`, the
#'   scale parameter
calc_loc_scale_params <- function(ps, qs, dist) {
  if (dist == "lnorm") {
    if (any(qs <= 0.0)) {
      stop("For dist = 'lnorm', all qs must be positive")
    }
    qs <- log(qs)
    qdst <- qnorm
  } else {
    qdst <- get(paste0("q", dist))
  }
  b <- suppressWarnings((qs[2] - qs[1]) / (qdst(ps[2]) - qdst(ps[1])))
  a <- suppressWarnings(qs[1] - b * qdst(ps[1]))
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
#' @importFrom stats dlnorm
#'
#' @return a function with parameters `x` and `log` that can be used to
#'   evaluate the density function (or its log) of the distribution in the
#'   specified location-scale family that has quantiles matching those in `ps`
#'   and `qs`
d_ext_factory <- function(ps, qs, dist) {
  c(a, b) %<-% calc_loc_scale_params(ps, qs, dist)

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
#' @importFrom stats plnorm
#'
#' @return a function with parameter `x` and `log.p` that can be used to
#'   evaluate the cumulative distribution function (or its log) of the
#'   distribution in the specified location-scale family that has quantiles
#'   matching those in `ps` and `qs`
p_ext_factory <- function(ps, qs, dist) {
  c(a, b) %<-% calc_loc_scale_params(ps, qs, dist)

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
#'
#' @importFrom stats qlnorm
#'
#' @return a function with parameter `p` that can be used to evaluate the
#'   quantile function of the distribution in the specified location-scale
#'   family that has quantiles matching those in `ps` and `qs`
q_ext_factory <- function(ps, qs, dist) {
  c(a, b) %<-% calc_loc_scale_params(ps, qs, dist)

  if (dist == "lnorm") {
    q_ext <- function(p) {
      return(qlnorm(p, meanlog = a, sdlog = b))
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

#' @importFrom stats runif
#' @importFrom utils head tail
#' @importFrom zeallot %<-%
NULL

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
d_ext_factory <- function(ps, qs, dist) {
    ddst <- get(paste0("d", dist))
    qdst <- get(paste0("q", dist))
    b <- (qs[2] - qs[1]) / (qdst(ps[2]) - qdst(ps[1]))
    a <- qs[1] - b * qdst(ps[1])

    d_ext <- function(x, log) {
        result <- ddst((x - a) / b, log = TRUE) - log(b)
        if (log) {
            return(result)
        } else {
            return(exp(result))
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
p_ext_factory <- function(ps, qs, dist) {
    pdst <- get(paste0("p", dist))
    qdst <- get(paste0("q", dist))
    b <- (qs[2] - qs[1]) / (qdst(ps[2]) - qdst(ps[1]))
    a <- qs[1] - b * qdst(ps[1])

    p_ext <- function(q, log.p = FALSE) {
        return(pdst((q - a) / b, log.p = log.p))
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
q_ext_factory <- function(ps, qs, dist) {
    qdst <- get(paste0("q", dist))
    b <- (qs[2] - qs[1]) / (qdst(ps[2]) - qdst(ps[1]))
    a <- qs[1] - b * qdst(ps[1])

    q_ext <- function(p) {
        return(a + b * qdst(p))
    }

    return(q_ext)
}


#' Clean up ps and qs provided by user: handle missing and unsorted values
#' 
#' @param ps vector of probability levels
#' @param qs vector of quantile values correponding to ps
#' 
#' @return named list with entries `ps` and `qs`
clean_ps_and_qs <- function(ps, qs) {
    # drop missing values for qs
    na_idx <- is.na(qs) | is.na(ps)
    if (any(na_idx)) {
        ps <- ps[!na_idx]
        qs <- qs[!na_idx]
    }

    # sort ps and qs
    ps <- sort(ps)
    qs <- sort(qs)

    return(list(ps=ps, qs=qs))
}


#' Creates a function that evaluates the probability density function of an
#' approximation to a distribution obtained by interpolating and extrapolating
#' from a set of quantiles of the distribution.
#'
#' @param ps vector of probability levels
#' @param qs vector of quantile values correponding to ps
#' @param interior_method method for monotonic cubic spline fit to the quantiles
#'   used for approximating the distribution on the interior of the quantiles.
#'   See [stats::splinefun()]. One of `"hyman"` or `"monoH.FC"`.
#' @param lower_tail_dist name of parametric distribution for the lower tail
#' @param upper_tail_dist name of parametric distribution for the upper tail
#'
#' @return a function with arguments `x` and `log` that can be used to evaluate
#'   the approximate density function (or its `log`) at the points `x`.
#' @export
make_d_fn <- function(ps, qs, interior_method = c("hyman", "monoH.FC"),
                      lower_tail_dist = "norm", upper_tail_dist = "norm") {
    interior_method <- match.arg(interior_method)

    c(ps, qs) %<-% clean_ps_and_qs(ps, qs)

    # throw an error if there are duplicated qs: the distribution is not
    # continuous
    if (any(duplicated(qs))) {
        stop("make_d_fn requires all values in qs to be unique")
    }

    # fit a monotonic spline to the qs and ps to approximate the distribution
    # on the interior
    interior_pdf <- stats::splinefun(qs, ps, method = interior_method)

    # approximate the pdf in the lower tail by extrapolating from the two
    # lowest quantiles within a location-scale family
    lower_pdf <- d_ext_factory(head(ps, 2), head(qs, 2), lower_tail_dist)

    # approximate the pdf in the upper tail by extrapolating from the two
    # largest quantiles within a location-scale family
    upper_pdf <- d_ext_factory(tail(ps, 2), tail(qs, 2), upper_tail_dist)

    d_fn <- function(x, log=FALSE) {
        # short-circuit if less than two unique value in qs
        if (length(unique(qs)) < 2) {
            return(rep(qs, length(x)))
        }

        # instantiate result
        result <- rep(NA_real_, length(x))

        # interior points
        interior_idx <- (x >= qs[1]) & (x <= tail(qs, 1))
        if (any(interior_idx)) {
            result[interior_idx] <- interior_pdf(x[interior_idx], deriv = 1)
            if (log) {
                result[interior_idx] <- log(result[interior_idx])
            }
        }

        # lower points
        lower_idx <- (x < qs[1])
        if (any(lower_idx)) {
            result[lower_idx] <- lower_pdf(x[lower_idx], log = log)
        }

        # upper points
        upper_idx <- (x > tail(qs, 1))
        if (any(upper_idx)) {
            result[upper_idx] <- upper_pdf(x[upper_idx], log = log)
        }

        return(result)
    }

    return(d_fn)
}



#' Creates a function that evaluates the cumulative distribution function of an
#' approximation to a distribution obtained by interpolating and extrapolating
#' from a set of quantiles of the distribution.
#'
#' @param ps vector of probability levels
#' @param qs vector of quantile values correponding to ps
#' @param interior_method method for monotonic cubic spline fit to the quantiles
#'   used for approximating the distribution on the interior of the quantiles.
#'   See [stats::splinefun()]. One of `"hyman"` or `"monoH.FC"`.
#' @param lower_tail_dist name of parametric distribution for the lower tail
#' @param upper_tail_dist name of parametric distribution for the upper tail
#'
#' @return a function with arguments `q` and `log.p` that can be used to
#'   evaluate the approximate cumulative distribution function (or its `log`)
#'   at the points `q`.
#' @export
make_p_fn <- function(ps, qs, interior_method = c("hyman", "monoH.FC"),
                      lower_tail_dist = "norm", upper_tail_dist = "norm") {
    interior_method <- match.arg(interior_method)

    c(ps, qs) %<-% clean_ps_and_qs(ps, qs)

    # fit a monotonic spline to the qs and ps to approximate the distribution
    # on the interior
    interior_cdf <- stats::splinefun(qs, ps, method = interior_method)

    # approximate the cdf in the lower tail by extrapolating from the two
    # lowest quantiles within a location-scale family
    lower_cdf <- p_ext_factory(head(ps, 2), head(qs, 2), lower_tail_dist)

    # approximate the pdf in the upper tail by extrapolating from the two
    # largest quantiles within a location-scale family
    upper_cdf <- p_ext_factory(tail(ps, 2), tail(qs, 2), upper_tail_dist)

    p_fn <- function(q, log.p=FALSE) {
        # short-circuit if less than two unique value in qs
        if (length(unique(qs)) < 2) {
            return(rep(qs[1], length(q)))
        }

        # instantiate result
        result <- rep(NA_real_, length(q))

        # interior points
        interior_idx <- (q >= qs[1]) & (q <= tail(qs, 1))
        if (any(interior_idx)) {
            result[interior_idx] <- interior_cdf(q[interior_idx], deriv = 0)
        }

        # lower points
        lower_idx <- (q < qs[1])
        if (any(lower_idx)) {
            result[lower_idx] <- lower_cdf(q[lower_idx])
        }

        # upper points
        upper_idx <- (q > tail(qs, 1))
        if (any(upper_idx)) {
            result[upper_idx] <- upper_cdf(q[upper_idx])
        }

        return(result)
    }

    return(p_fn)
}



#' Creates a function that evaluates the quantile function of an approximation
#' to a distribution obtained by interpolating and extrapolating from a set of
#' quantiles of the distribution.
#'
#' @param ps vector of probability levels
#' @param qs vector of quantile values correponding to ps
#' @param interior_method method for monotonic cubic spline fit to the quantiles
#'   used for approximating the distribution on the interior of the quantiles.
#'   See [stats::splinefun()]. One of `"hyman"` or `"monoH.FC"`.
#' @param lower_tail_dist name of parametric distribution for the lower tail
#' @param upper_tail_dist name of parametric distribution for the upper tail
#'
#' @return a function with argument `p` that can be used to calculate quantiles
#'   of the approximated distribution at the probability levels `p`.
#'
#' @export
make_q_fn <- function(ps, qs, interior_method = c("hyman", "monoH.FC"),
                      lower_tail_dist = "norm", upper_tail_dist = "norm") {
    interior_method <- match.arg(interior_method)

    c(ps, qs) %<-% clean_ps_and_qs(ps, qs)

    # fit a monotonic spline to the ps and qs to approximate the distribution
    # on the interior
    interior_qf <- stats::splinefun(ps, qs, method = interior_method)

    # approximate the quantile function in the lower tail by extrapolating from
    # the two lowest quantiles within a location-scale family
    lower_qf <- q_ext_factory(head(ps, 2), head(qs, 2), lower_tail_dist)

    # approximate the quantile function in the upper tail by extrapolating from
    # the two largest quantiles within a location-scale family
    upper_qf <- q_ext_factory(tail(ps, 2), tail(qs, 2), upper_tail_dist)

    q_fn <- function(p) {
        # short-circuit if less than two unique value in qs
        if (length(unique(qs)) < 2) {
            return(rep(qs, length(p)))
        }

        # instantiate result
        result <- rep(NA_real_, length(p))

        # interior points
        interior_idx <- (p >= ps[1]) & (p <= tail(ps, 1))
        if (any(interior_idx)) {
            result[interior_idx] <- interior_qf(p[interior_idx])
        }

        # lower points
        lower_idx <- (p < ps[1])
        if (any(lower_idx)) {
            result[lower_idx] <- lower_qf(p[lower_idx])
        }

        # upper points
        upper_idx <- (p > tail(ps, 1))
        if (any(upper_idx)) {
            result[upper_idx] <- upper_qf(p[upper_idx])
        }

        return(result)
    }

    return(q_fn)
}


#' Creates a function that evaluates the quantile function of an approximation
#' to a distribution obtained by interpolating and extrapolating from a set of
#' quantiles of the distribution.
#'
#' @param ps vector of probability levels
#' @param qs vector of quantile values correponding to ps
#' @param interior_method method for monotonic cubic spline fit to the quantiles
#'   used for approximating the distribution on the interior of the quantiles.
#'   See [stats::splinefun()]. One of `"hyman"` or `"monoH.FC"`.
#' @param lower_tail_dist name of parametric distribution for the lower tail
#' @param upper_tail_dist name of parametric distribution for the upper tail
#'
#' @return a function with argument `n` that can be used to generate a sample of
#'   size `n` from the approximated distribution.
#' @export
make_r_fn <- function(ps, qs, interior_method = c("hyman", "monoH.FC"),
                      lower_tail_dist = "norm", upper_tail_dist = "norm") {
    interior_method <- match.arg(interior_method)
    q_fn <- make_q_fn(ps, qs, interior_method, lower_tail_dist, upper_tail_dist)

    r_fn <- function(n) {
        u <- runif(n)
        return(q_fn(u))
    }

    return(r_fn)
}

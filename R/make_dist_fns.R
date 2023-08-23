#' @importFrom stats runif
#' @importFrom utils head tail
#' @importFrom zeallot %<-%
#' @importFrom splines backSpline
NULL


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

    return(list(ps = ps, qs = qs))
}


#' Creates a function that evaluates the probability density function of an
#' approximation to a distribution obtained by interpolating and extrapolating
#' from a set of quantiles of the distribution.
#'
#' @param ps vector of probability levels
#' @param qs vector of quantile values correponding to ps
#' @param interior_method method for interpolating the distribution on the
#'    interior of the provided `qs`. This package provides one method for this,
#'    `"spline_cdf"`. The user may also provide a custom function; see the
#'    details for more.
#' @param interior_args an optional named list of arguments that are passed
#'    on to the `interior_method`
#' @param lower_tail_dist name of parametric distribution for the lower tail
#' @param upper_tail_dist name of parametric distribution for the upper tail
#' @param dup_tol numeric tolerance for identifying duplicated values indicating
#'   a discrete component of the distribution. If there is a run of values where
#'   each consecutive pair is closer together than the tolerance, all are
#'   labeled as duplicates even if not all values in the run are within the
#'   tolerance.
#' @param zero_tol numeric tolerance for identifying values in `qs` that are
#'   (approximately) zero.
#'
#' @details The default `interior_method`, `"spline_cdf"`, represents the
#'    distribution as a sum of a discrete component at any points where there
#'    are duplicated `qs` for multiple different `ps` and a continuous component
#'    that is estimated by using a monotonic cubic spline that interpolates the
#'    provided `(q, p)` pairs as an estimate of the cdf. The density function is
#'    then obtained by differentiating this estimate of the cdf.
#' 
#'    Optionally, the user may provide another function that accepts arguments
#'    `ps`, `qs`, `lower_tail_dist`, `upper_tail_dist`, and `fn_type` (which
#'    will be either `"d"`, `"p"`, or `"q"`), and optionally additional named
#'    arguments to be specified via `interior_args`. This function should return
#'    a function with arguments `x`, `log` that evaluates the pdf or its
#'    logarithm.
#'
#' @return a function with arguments `x` and `log` that can be used to evaluate
#'   the approximate density function (or its `log`) at the points `x`.
#' @export
make_d_fn <- function(ps, qs,
                      interior_method = "spline_cdf",
                      interior_args = list(),
                      lower_tail_dist = "norm", upper_tail_dist = "norm",
                      dup_tol = 1e-6,
                      zero_tol = 1e-12) {
    interior_method <- match.arg(interior_method)

    c(ps, qs) %<-% clean_ps_and_qs(ps, qs)

    # throw an error if there are duplicated qs, or fewer than 1 distinct qs:
    # the distribution is not continuous
    if (any(duplicated_tol(qs, tol = dup_tol))) {
        stop("make_d_fn requires all values in qs to be unique")
    }
    if (length(qs) < 2) {
        stop("make_d_fn requires at least two unique qs")
    }

    if ((lower_tail_dist == "lnorm" || upper_tail_dist == "lnorm") &&
            any(qs < 0.0)) {
        stop("For lognormal tail distributions, all qs must be non-negative")
    }

    # approximate the pdf on the interior by interpolating quantiles
    interior_args <- c(
        list(ps = ps, qs = qs, lower_tail_dist = lower_tail_dist,
             upper_tail_dist = upper_tail_dist, fn_type = "d"),
        interior_args)
    if (!(dup_tol %in% names(interior_args))) {
        interior_args <- c(
            interior_args,
            list(dup_tol = dup_tol)
        )
    }
    if (!(zero_tol %in% names(interior_args))) {
        interior_args <- c(
            interior_args,
            list(zero_tol = zero_tol)
        )
    }
    interior_pdf <- do.call(interior_method, args = interior_args)

    # approximate the pdf in the lower/upper tail by extrapolating from the two
    # lowest/largest quantiles within a location-scale family
    zero_discrete <- lower_tail_dist == "lnorm" || upper_tail_dist == "lnorm"
    lower_pdf <- d_ext_factory(head(ps, 2), head(qs, 2), lower_tail_dist,
                               dup_tol, zero_tol, zero_discrete, TRUE)
    upper_pdf <- d_ext_factory(tail(ps, 2), tail(qs, 2), upper_tail_dist,
                               dup_tol, zero_tol, zero_discrete, FALSE)

    d_fn <- function(x, log=FALSE) {
        # instantiate result
        result <- rep(NA_real_, length(x))

        # interior points
        interior_idx <- (x >= qs[1]) & (x <= tail(qs, 1))
        if (any(interior_idx)) {
            result[interior_idx] <- interior_pdf(x[interior_idx], log = log)
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
#' @param interior_method method for interpolating the distribution on the
#'    interior of the provided `qs`. This package provides one method for this,
#'    `"spline_cdf"`. The user may also provide a custom function; see the
#'    details for more.
#' @param interior_args an optional named list of arguments that are passed
#'    on to the `interior_method`
#' @param lower_tail_dist name of parametric distribution for the lower tail
#' @param upper_tail_dist name of parametric distribution for the upper tail
#' @param dup_tol numeric tolerance for identifying duplicated values indicating
#'   a discrete component of the distribution. If there is a run of values where
#'   each consecutive pair is closer together than the tolerance, all are
#'   labeled as duplicates even if not all values in the run are within the
#'   tolerance.
#' @param zero_tol numeric tolerance for identifying values in `qs` that are
#'   (approximately) zero.
#'
#' @details The default `interior_method`, `"spline_cdf"`, represents the
#'    distribution as a sum of a discrete component at any points where there
#'    are duplicated `qs` for multiple different `ps` and a continuous component
#'    that is estimated by using a monotonic cubic spline that interpolates the
#'    provided `(q, p)` pairs as an estimate of the cdf.
#'
#'    Optionally, the user may provide another function that accepts arguments
#'    `ps`, `qs`, `lower_tail_dist`, `upper_tail_dist`, and `fn_type` (which
#'    will be either `"d"`, `"p"`, or `"q"`), and optionally additional named
#'    arguments to be specified via `interior_args`. This function should return
#'    a function with arguments `x`, `log` that evaluates the pdf or its
#'    logarithm.
#'
#' @return a function with arguments `q` and `log.p` that can be used to
#'   evaluate the approximate cumulative distribution function (or its `log`)
#'   at the points `q`.
#' @export
make_p_fn <- function(ps, qs,
                      interior_method = "spline_cdf",
                      interior_args = list(),
                      lower_tail_dist = "norm", upper_tail_dist = "norm",
                      dup_tol = 1e-6,
                      zero_tol = 1e-12) {
    interior_method <- match.arg(interior_method)

    c(ps, qs) %<-% clean_ps_and_qs(ps, qs)

    if ((lower_tail_dist == "lnorm" || upper_tail_dist == "lnorm") &&
            any(qs < 0.0)) {
        stop("For lognormal tail distributions, all qs must be non-negative")
    }

    # short-circuit if less than two unique value in qs
    uq <- unique_tol(qs, tol = dup_tol)
    if (length(uq) == 1) {
        p_fn <- function(q, log.p = FALSE) {
            result <- as.numeric(q >= uq)
            if (log.p) {
                return(log(result))
            } else {
                return(result)
            }
        }
        return(p_fn)
    }

    # approximate the cdf on the interior by interpolating quantiles
    interior_args <- c(
        list(ps = ps, qs = qs, lower_tail_dist = lower_tail_dist,
             upper_tail_dist = upper_tail_dist, fn_type = "p"),
        interior_args)
    if (!(dup_tol %in% names(interior_args))) {
        interior_args <- c(
            interior_args,
            list(dup_tol = dup_tol)
        )
    }
    if (!(zero_tol %in% names(interior_args))) {
        interior_args <- c(
            interior_args,
            list(zero_tol = zero_tol)
        )
    }
    interior_cdf <- do.call(interior_method, args = interior_args)

    # approximate the cdf in the lower/upper tail by extrapolating from the two
    # lowest/largest quantiles within a location-scale family
    zero_discrete <- lower_tail_dist == "lnorm" || upper_tail_dist == "lnorm"
    lower_cdf <- p_ext_factory(head(ps, 2), head(qs, 2), lower_tail_dist,
                               dup_tol, zero_tol, zero_discrete, TRUE)
    upper_cdf <- p_ext_factory(tail(ps, 2), tail(qs, 2), upper_tail_dist,
                               dup_tol, zero_tol, zero_discrete, FALSE)

    # TODO: if zero_discrete, determination of interior/lower/upper points
    # should be by comparison to the non-zero qs. If we fit spline to non-zero
    # qs, we don't want to extrapolate beyond them.
    # this means we need to handle point masses "in the tails" as well!
    # if there is a pm at 0, we need to add it in here.
    # we should:
    #  - split disc/cont ps and qs out here
    #  - if there is any part of the distribution that is continuous, get an
    #    estimate of that based on cont_ps and cont_qs, including interior and tails
    #  - add discrete masses wherever they are, tails or interior
    # basically move logic related to adjusting for point masses to here from spline_cdf
    #
    # OR
    #
    # special one-off handling of point mass at 0
    p_fn <- function(q, log.p = FALSE) {
        # instantiate result
        result <- rep(NA_real_, length(q))

        # interior points
        interior_idx <- (q >= qs[1]) & (q <= tail(qs, 1))
        if (any(interior_idx)) {
            result[interior_idx] <- interior_cdf(q[interior_idx], log.p = log.p)
        }

        # lower points
        lower_idx <- (q < qs[1])
        if (any(lower_idx)) {
            result[lower_idx] <- lower_cdf(q[lower_idx], log.p = log.p)
        }

        # upper points
        upper_idx <- (q > tail(qs, 1))
        if (any(upper_idx)) {
            result[upper_idx] <- upper_cdf(q[upper_idx], log.p = log.p)
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
#' @param interior_method method for interpolating the distribution on the
#'    interior of the provided `qs`. This package provides one method for this,
#'    `"spline_cdf"`. The user may also provide a custom function; see the
#'    details for more.
#' @param interior_args an optional named list of arguments that are passed
#'    on to the `interior_method`
#' @param lower_tail_dist name of parametric distribution for the lower tail
#' @param upper_tail_dist name of parametric distribution for the upper tail
#' @param dup_tol numeric tolerance for identifying duplicated values indicating
#'   a discrete component of the distribution. If there is a run of values where
#'   each consecutive pair is closer together than the tolerance, all are
#'   labeled as duplicates even if not all values in the run are within the
#'   tolerance.
#' @param zero_tol numeric tolerance for identifying values in `qs` that are
#'   (approximately) zero.
#'
#' @details The default `interior_method`, `"spline_cdf"`, represents the
#'    distribution as a sum of a discrete component at any points where there
#'    are duplicated `qs` for multiple different `ps` and a continuous component
#'    that is estimated by using a monotonic cubic spline that interpolates the
#'    provided `(q, p)` pairs as an estimate of the cdf. The quantile function
#'    is then obtained by inverting this estimate of the cdf.
#'
#'    Optionally, the user may provide another function that accepts arguments
#'    `ps`, `qs`, `lower_tail_dist`, `upper_tail_dist`, and `fn_type` (which
#'    will be either `"d"`, `"p"`, or `"q"`), and optionally additional named
#'    arguments to be specified via `interior_args`. This function should return
#'    a function with argument `p` that evaluates the quantile function.
#'
#' @return a function with argument `p` that can be used to calculate quantiles
#'   of the approximated distribution at the probability levels `p`.
#'
#' @export
make_q_fn <- function(ps, qs,
                      interior_method = "spline_cdf",
                      interior_args = list(),
                      lower_tail_dist = "norm", upper_tail_dist = "norm",
                      dup_tol = 1e-6,
                      zero_tol = 1e-12) {
    interior_method <- match.arg(interior_method)

    c(ps, qs) %<-% clean_ps_and_qs(ps, qs)

    if ((lower_tail_dist == "lnorm" || upper_tail_dist == "lnorm") &&
            any(qs < 0.0)) {
        stop("For lognormal tail distributions, all qs must be non-negative")
    }

    # short-circuit if less than two unique value in qs
    uq <- unique_tol(qs, tol = dup_tol)
    if (length(uq) == 1) {
        q_fn <- function(p) {
            return(rep(uq, length(p)))
        }
        return(q_fn)
    }

    # approximate the pdf on the interior by interpolating quantiles
    interior_args <- c(
        list(ps = ps, qs = qs, lower_tail_dist = lower_tail_dist,
             upper_tail_dist = upper_tail_dist, fn_type = "q"),
        interior_args)
    if (!(dup_tol %in% names(interior_args))) {
        interior_args <- c(
            interior_args,
            list(dup_tol = dup_tol)
        )
    }
    if (!(zero_tol %in% names(interior_args))) {
        interior_args <- c(
            interior_args,
            list(zero_tol = zero_tol)
        )
    }
    interior_qf <- do.call(interior_method, args = interior_args)

    # approximate the quantile function in the lower/upper tail by extrapolating
    # from the two lowest/largest quantiles within a location-scale family
    zero_discrete <- lower_tail_dist == "lnorm" || upper_tail_dist == "lnorm"
    lower_qf <- q_ext_factory(head(ps, 2), head(qs, 2), lower_tail_dist,
                              dup_tol, zero_tol, zero_discrete, TRUE)
    upper_qf <- q_ext_factory(tail(ps, 2), tail(qs, 2), upper_tail_dist,
                              dup_tol, zero_tol, zero_discrete, TRUE)

    q_fn <- function(p) {
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


#' Creates a function that generates random deviates from an approximation
#' to a distribution obtained by interpolating and extrapolating from a set of
#' quantiles of the distribution.
#'
#' @param ps vector of probability levels
#' @param qs vector of quantile values correponding to ps
#' @param interior_method method for interpolating the distribution on the
#'    interior of the provided `qs`. This package provides one method for this,
#'    `"spline_cdf"`. The user may also provide a custom function; see the
#'    details for more.
#' @param interior_args an optional named list of arguments that are passed
#'    on to the `interior_method`
#' @param lower_tail_dist name of parametric distribution for the lower tail
#' @param upper_tail_dist name of parametric distribution for the upper tail
#' @param dup_tol numeric tolerance for identifying duplicated values indicating
#'   a discrete component of the distribution. If there is a run of values where
#'   each consecutive pair is closer together than the tolerance, all are
#'   labeled as duplicates even if not all values in the run are within the
#'   tolerance.
#' @param zero_tol numeric tolerance for identifying values in `qs` that are
#'   (approximately) zero.
#'
#' @details The default `interior_method`, `"spline_cdf"`, represents the
#'    distribution as a sum of a discrete component at any points where there
#'    are duplicated `qs` for multiple different `ps` and a continuous component
#'    that is estimated by using a monotonic cubic spline that interpolates the
#'    provided `(q, p)` pairs as an estimate of the cdf. The quantile function
#'    is then obtained by inverting this estimate of the cdf.
#'
#'    Optionally, the user may provide another function that accepts arguments
#'    `ps`, `qs`, `lower_tail_dist`, `upper_tail_dist`, and `fn_type` (which
#'    will be either `"d"`, `"p"`, or `"q"`), and optionally additional named
#'    arguments to be specified via `interior_args`. This function should return
#'    a function with argument `p` that evaluates the quantile function.
#'
#' @return a function with argument `n` that can be used to generate a sample of
#'   size `n` from the approximated distribution.
#' @export
make_r_fn <- function(ps, qs,
                      interior_method = "spline_cdf",
                      interior_args = list(),
                      lower_tail_dist = "norm", upper_tail_dist = "norm",
                      dup_tol = 1e-6,
                      zero_tol = 1e-12) {
    interior_method <- match.arg(interior_method)
    q_fn <- make_q_fn(ps, qs, interior_method, interior_args, lower_tail_dist,
                      upper_tail_dist, dup_tol, zero_tol)

    r_fn <- function(n) {
        u <- runif(n)
        return(q_fn(u))
    }

    return(r_fn)
}

#' @importFrom stats runif
#' @importFrom utils head tail
#' @importFrom zeallot %<-%
#' @importFrom splines backSpline
NULL

disc_weight <- disc_ps <- disc_qs <- cont_ps <- cont_qs <- disc_ps_range <- NULL

#' Clean up ps and qs provided by user: handle missing and unsorted values
#'
#' @param ps vector of probability levels
#' @param qs vector of quantile values correponding to ps
#'
#' @return named list with entries `ps` and `qs`
clean_ps_and_qs <- function(ps, qs) {
  checkmate::assert_numeric(ps, lower = 0, upper = 1)
  checkmate::assert_numeric(qs)

  if (length(ps) != length(qs)) {
    stop("'ps' and 'qs' must have the same length.")
  }

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
#'   interior of the provided `qs`. This package provides one method for this,
#'   `"spline_cdf"`. The user may also provide a custom function; see the
#'   details for more.
#' @param interior_args an optional named list of arguments that are passed
#'   on to the `interior_method`
#' @param tail_dist name of parametric distribution for the tails
#' @param dup_tol numeric tolerance for identifying duplicated values indicating
#'   a discrete component of the distribution. If there is a run of values where
#'   each consecutive pair is closer together than the tolerance, all are
#'   labeled as duplicates even if not all values in the run are within the
#'   tolerance.
#' @param zero_tol numeric tolerance for identifying values in `qs` that are
#'   (approximately) zero.
#'
#' @details The default `interior_method`, `"spline_cdf"`, represents the
#'   distribution as a sum of a discrete component at any points where there
#'   are duplicated `qs` for multiple different `ps` and a continuous component
#'   that is estimated by using a monotonic cubic spline that interpolates the
#'   provided `(q, p)` pairs as an estimate of the cdf. The density function is
#'   then obtained by differentiating this estimate of the cdf.
#'
#'   Optionally, the user may provide another function that accepts arguments
#'   `ps`, `qs`, `tail_dist`, and `fn_type` (which will be either `"d"`, `"p"`,
#'   or `"q"`), and optionally additional named arguments to be specified via
#'   `interior_args`. This function should return a function with arguments
#'   `x`, `log` that evaluates the pdf or its logarithm.
#'
#' @return a function with arguments `x` and `log` that can be used to evaluate
#'   the approximate density function (or its `log`) at the points `x`.
#' @export
make_d_fn <- function(ps, qs,
                      interior_method = "spline_cdf",
                      interior_args = list(),
                      tail_dist = "norm",
                      dup_tol = 1e-6, zero_tol = 1e-12) {
  interior_method <- match.arg(interior_method)

  c(ps, qs) %<-% clean_ps_and_qs(ps, qs)

  # if tail_dist is "lnorm", we treat quantiles of 0 as indicative of a
  # discrete point mass at 0, using a set up like a hurdle model
  is_hurdle <- tail_dist == "lnorm"

  # split ps and qs into discrete and continuous part, along with the
  # weight given to the discrete part
  c(disc_weight, disc_ps, disc_qs, cont_ps, cont_qs, disc_ps_range) %<-%
    split_disc_cont_ps_qs(ps, qs, dup_tol = dup_tol, zero_tol = zero_tol,
                          is_hurdle = is_hurdle)

  # throw an error if there is any discontinuous part of the distribution:
  # we cannot estimate a density function
  if (disc_weight > 0) {
    stop("make_d_fn requires the distribution to be continuous, but a discrete component was detected")
  }

  # approximate the pdf on the interior by interpolating quantiles
  interior_args <- c(
    list(ps = cont_ps, qs = cont_qs, tail_dist = tail_dist, fn_type = "d"),
    interior_args
  )
  interior_pdf <- do.call(interior_method, args = interior_args)

  # approximate the pdf in the lower tail by extrapolating from the two
  # lowest quantiles within a location-scale family
  lower_pdf <- d_ext_factory(head(cont_ps, 2), head(cont_qs, 2), tail_dist)

  # approximate the pdf in the upper tail by extrapolating from the two
  # largest quantiles within a location-scale family
  upper_pdf <- d_ext_factory(tail(cont_ps, 2), tail(cont_qs, 2), tail_dist)

  d_fn <- function(x, log = FALSE) {
    checkmate::assert_numeric(x)

    # instantiate result
    result <- rep(NA_real_, length(x))

    # interior points
    interior_idx <- (x >= cont_qs[1]) & (x <= tail(cont_qs, 1))
    if (any(interior_idx)) {
      result[interior_idx] <- interior_pdf(x[interior_idx], log = log)
    }

    # lower points
    lower_idx <- (x < cont_qs[1])
    if (any(lower_idx)) {
      result[lower_idx] <- lower_pdf(x[lower_idx], log = log)
    }

    # upper points
    upper_idx <- (x > tail(cont_qs, 1))
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
#'   interior of the provided `qs`. This package provides one method for this,
#'   `"spline_cdf"`. The user may also provide a custom function; see the
#'   details for more.
#' @param interior_args an optional named list of arguments that are passed
#'   on to the `interior_method`
#' @param tail_dist name of parametric distribution for the tails
#' @param dup_tol numeric tolerance for identifying duplicated values indicating
#'   a discrete component of the distribution. If there is a run of values where
#'   each consecutive pair is closer together than the tolerance, all are
#'   labeled as duplicates even if not all values in the run are within the
#'   tolerance.
#' @param zero_tol numeric tolerance for identifying values in `qs` that are
#'   (approximately) zero.
#'
#' @details The default `interior_method`, `"spline_cdf"`, represents the
#'   distribution as a sum of a discrete component at any points where there
#'   are duplicated `qs` for multiple different `ps` and a continuous component
#'   that is estimated by using a monotonic cubic spline that interpolates the
#'   provided `(q, p)` pairs as an estimate of the cdf.
#'
#'   Optionally, the user may provide another function that accepts arguments
#'   `ps`, `qs`, `tail_dist`, and `fn_type` (which will be either `"d"`, `"p"`,
#'   or `"q"`), and optionally additional named arguments to be specified via
#'   `interior_args`. This function should return a function with arguments
#'   `x`, `log` that evaluates the pdf or its logarithm.
#'
#' @return a function with arguments `q` and `log.p` that can be used to
#'   evaluate the approximate cumulative distribution function (or its `log`)
#'   at the points `q`.
#' @export
make_p_fn <- function(ps, qs,
                      interior_method = "spline_cdf",
                      interior_args = list(),
                      tail_dist = "norm",
                      dup_tol = 1e-6, zero_tol = 1e-12) {
  interior_method <- match.arg(interior_method)

  c(ps, qs) %<-% clean_ps_and_qs(ps, qs)

  # if tail_dist is "lnorm", we treat quantiles of 0 as indicative of a
  # discrete point mass at 0, using a set up like a hurdle model
  is_hurdle <- tail_dist == "lnorm"

  # split ps and qs into discrete and continuous part, along with the
  # weight given to the discrete part
  c(disc_weight, disc_ps, disc_qs, cont_ps, cont_qs, disc_ps_range) %<-%
    split_disc_cont_ps_qs(ps, qs, dup_tol = dup_tol, zero_tol = zero_tol,
                          is_hurdle = is_hurdle)

  if (disc_weight < 1.0) {
    # approximate the cdf on the interior by interpolating quantiles
    interior_args <- c(
      list(ps = cont_ps, qs = cont_qs, tail_dist = tail_dist, fn_type = "p"),
      interior_args
    )
    interior_cdf <- do.call(interior_method, args = interior_args)

    # approximate the cdf in the lower tail by extrapolating from the two
    # lowest quantiles within a location-scale family
    if (min(cont_ps) > 0.0) {
      lower_cdf <- p_ext_factory(head(cont_ps, 2), head(cont_qs, 2), tail_dist)
    }

    # approximate the pdf in the upper tail by extrapolating from the two
    # largest quantiles within a location-scale family
    if (max(cont_ps) < 1.0) {
      upper_cdf <- p_ext_factory(tail(cont_ps, 2), tail(cont_qs, 2), tail_dist)
    }
  }

  p_fn <- function(q, log.p = FALSE) {
    checkmate::assert_numeric(q)

    # instantiate result
    result <- rep(0.0, length(q))
    log.p_direct <- FALSE

    if (disc_weight < 1.0) {
      # directly use log.p in interior and exterior cdf methods?
      # if possible, probably saves a little numerical precision in tails
      # but we can't (easily) do this is if disc_weight > 0, in which case
      # below we will be adding probabilities from a discrete component
      log.p_direct <- (log.p && disc_weight == 0.0)

      # interior points
      interior_idx <- (q >= cont_qs[1]) & (q <= tail(cont_qs, 1))
      if (any(interior_idx)) {
        result[interior_idx] <- interior_cdf(q[interior_idx],
                                             log.p = log.p_direct)
      }

      # lower points
      lower_idx <- (q < cont_qs[1])
      if (any(lower_idx)) {
        if (min(cont_ps) > 0.0) {
          result[lower_idx] <- lower_cdf(q[lower_idx], log.p = log.p_direct)
        } # else set to 0, which is how we initialized it
      }

      # upper points
      upper_idx <- (q > tail(cont_qs, 1))
      if (any(upper_idx)) {
        if (max(cont_ps) < 1.0) {
          result[upper_idx] <- upper_cdf(q[upper_idx], log.p = log.p_direct)
        } else {
          result[upper_idx] <- 1.0
        }
      }

      # adjust by weight for continuous component
      result <- result * (1 - disc_weight)
    }

    # add discrete probabilities
    for (i in seq_along(disc_ps)) {
      inds <- (q >= disc_qs[i])
      result[inds] <- result[inds] + disc_ps[i] * disc_weight
    }

    # ensure result is within interval [0, 1]
    result <- pmin(pmax(result, 0), 1)

    # return, handling log.p argument
    if (log.p_direct) {
      # result is already on log scale
      return(result)
    } else if (log.p) {
      return(log(result))
    } else {
      return(result)
    }
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
#'   interior of the provided `qs`. This package provides one method for this,
#'   `"spline_cdf"`. The user may also provide a custom function; see the
#'   details for more.
#' @param interior_args an optional named list of arguments that are passed
#'   on to the `interior_method`
#' @param tail_dist name of parametric distribution for the tails
#' @param dup_tol numeric tolerance for identifying duplicated values indicating
#'   a discrete component of the distribution. If there is a run of values where
#'   each consecutive pair is closer together than the tolerance, all are
#'   labeled as duplicates even if not all values in the run are within the
#'   tolerance.
#' @param zero_tol numeric tolerance for identifying values in `qs` that are
#'   (approximately) zero.
#'
#' @details The default `interior_method`, `"spline_cdf"`, represents the
#'   distribution as a sum of a discrete component at any points where there
#'   are duplicated `qs` for multiple different `ps` and a continuous component
#'   that is estimated by using a monotonic cubic spline that interpolates the
#'   provided `(q, p)` pairs as an estimate of the cdf. The quantile function
#'   is then obtained by inverting this estimate of the cdf.
#'
#'   Optionally, the user may provide another function that accepts arguments
#'   `ps`, `qs`, `tail_dist`, and `fn_type` (which will be either `"d"`, `"p"`,
#'   or `"q"`), and optionally additional named arguments to be specified via
#'   `interior_args`. This function should return a function with argument `p`
#'   that evaluates the quantile function.
#'
#' @return a function with argument `p` that can be used to calculate quantiles
#'   of the approximated distribution at the probability levels `p`.
#'
#' @export
make_q_fn <- function(ps, qs,
                      interior_method = "spline_cdf",
                      interior_args = list(),
                      tail_dist = "norm",
                      dup_tol = 1e-6, zero_tol = 1e-12) {
  interior_method <- match.arg(interior_method)

  c(ps, qs) %<-% clean_ps_and_qs(ps, qs)

  # if tail_dist is "lnorm", we treat quantiles of 0 as indicative of a
  # discrete point mass at 0, using a set up like a hurdle model
  is_hurdle <- tail_dist == "lnorm"

  # split ps and qs into discrete and continuous part, along with the
  # weight given to the discrete part
  c(disc_weight, disc_ps, disc_qs, cont_ps, cont_qs, disc_ps_range) %<-%
    split_disc_cont_ps_qs(ps, qs, dup_tol = dup_tol, zero_tol = zero_tol,
                          is_hurdle = is_hurdle)

  if (disc_weight < 1.0) {
    # approximate the pdf on the interior by interpolating quantiles
    interior_args <- c(
      list(ps = cont_ps, qs = cont_qs, tail_dist = tail_dist, fn_type = "q"),
      interior_args
    )
    interior_qf <- do.call(interior_method, args = interior_args)

    # approximate the quantile function in the lower tail by extrapolating
    # from the two lowest quantiles within a location-scale family
    if (min(cont_ps) > 0.0) {
      lower_qf <- q_ext_factory(head(cont_ps, 2), head(cont_qs, 2),
                                tail_dist)
    }

    # approximate the quantile function in the upper tail by extrapolating
    # from the two largest quantiles within a location-scale family
    if (max(cont_ps) < 1.0) {
      upper_qf <- q_ext_factory(tail(cont_ps, 2), tail(cont_qs, 2), tail_dist)
    }
  }

  q_fn <- function(p) {
    checkmate::assert_numeric(p, lower = 0, upper = 1)

    # instantiate result
    result <- rep(NA_real_, length(p))

    # discrete part of distribution
    # each discrete component has an associated ((min_p, max_p), q) tuple
    # where for any p in [min_p, max_p], the quantile is q
    # along the way, we update the continuous probabilities by subtracting
    # the discrete point mass probability
    cont_inds <- rep(TRUE, length(p))
    cont_p <- p
    for (i in seq_along(disc_ps_range)) {
      dpr <- disc_ps_range[[i]]
      inds <- (p >= dpr[1]) & (p <= dpr[2])
      result[inds] <- disc_qs[i]
      cont_inds[inds] <- FALSE

      inds <- (p > dpr[2])
      cont_p[inds] <- cont_p[inds] - disc_ps[i] * disc_weight
    }

    # continuous part of the distribution
    if (disc_weight < 1.0) {
      cont_p <- cont_p / (1 - disc_weight)
      # address potential floating point errors: for example,
      # dividing by (1 - disc_weight) could take an input p = 1 to a value
      # slightly greater than 1. here, we ensure that any input p's that
      # were in [0, 1] are still in [0, 1] after adjustment
      cont_p[p >= 0] <- pmax(cont_p[p >= 0], 0)
      cont_p[p <= 1] <- pmin(cont_p[p <= 1], 1)

      # interior points
      interior_idx <- cont_inds & (cont_p >= cont_ps[1]) &
        (cont_p <= tail(cont_ps, 1))
      if (any(interior_idx)) {
        result[interior_idx] <- interior_qf(cont_p[interior_idx])
      }

      # lower points
      # no action required if min(cont_ps) == 0.0 because in that case,
      # we will never have p < cont_ps[1].  This means the outer check
      # is redundant, but may be helpful to have consistency
      if (min(cont_ps) > 0.0) {
        lower_idx <- cont_inds & (cont_p < cont_ps[1])
        if (any(lower_idx)) {
          result[lower_idx] <- lower_qf(cont_p[lower_idx])
        }
      }

      # upper points
      # no action required if max(cont_ps) == 1.0 because in that case,
      # we will never have p > tail(cont_ps, 1).  This means the outer
      # check is redundant, but may be helpful to have consistency
      if (max(cont_ps) < 1.0) {
        upper_idx <- cont_inds & (cont_p > tail(cont_ps, 1))
        if (any(upper_idx)) {
          result[upper_idx] <- upper_qf(cont_p[upper_idx])
        }
      }
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
#'   interior of the provided `qs`. This package provides one method for this,
#'   `"spline_cdf"`. The user may also provide a custom function; see the
#'   details for more.
#' @param interior_args an optional named list of arguments that are passed
#'   on to the `interior_method`
#' @param tail_dist name of parametric distribution for the tails
#' @param dup_tol numeric tolerance for identifying duplicated values indicating
#'   a discrete component of the distribution. If there is a run of values where
#'   each consecutive pair is closer together than the tolerance, all are
#'   labeled as duplicates even if not all values in the run are within the
#'   tolerance.
#' @param zero_tol numeric tolerance for identifying values in `qs` that are
#'   (approximately) zero.
#'
#' @details The default `interior_method`, `"spline_cdf"`, represents the
#'   distribution as a sum of a discrete component at any points where there
#'   are duplicated `qs` for multiple different `ps` and a continuous component
#'   that is estimated by using a monotonic cubic spline that interpolates the
#'   provided `(q, p)` pairs as an estimate of the cdf. The quantile function
#'   is then obtained by inverting this estimate of the cdf.
#'
#'   Optionally, the user may provide another function that accepts arguments
#'   `ps`, `qs`, `tail_dist`, and `fn_type` (which will be either `"d"`, `"p"`,
#'   or `"q"`), and optionally additional named arguments to be specified via
#'   `interior_args`. This function should return a function with argument `p`
#'   that evaluates the quantile function.
#'
#' @return a function with argument `n` that can be used to generate a sample of
#'   size `n` from the approximated distribution.
#' @export
make_r_fn <- function(ps, qs,
                      interior_method = "spline_cdf",
                      interior_args = list(),
                      tail_dist = "norm",
                      dup_tol = 1e-6, zero_tol = 1e-12) {
  interior_method <- match.arg(interior_method)
  q_fn <- make_q_fn(ps, qs, interior_method, interior_args, tail_dist)

  r_fn <- function(n) {
    checkmate::assert_integerish(n, lower = 0, any.missing = FALSE, len = 1)
    u <- runif(n)
    return(q_fn(u))
  }

  return(r_fn)
}

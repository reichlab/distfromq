dup_run_starts <- dup_run_ends <- NULL

#' Identify duplicated values in a sorted numeric vector, where comparison is
#' up to a specified numeric tolerance. If there is a run of values where each
#' consecutive pair is closer together than the tolerance, all are labeled as
#' duplicates even if not all values in the run are within the tolerance.
#'
#' @param x a numeric vector in which to identify duplicates
#' @param tol numeric tolerance for identifying duplicates
#' @param incl_first boolean indicator of whether or not the first entry in a
#'   run of duplicates should be indicated as a duplicate. `FALSE` mirrors the
#'   behavior of the base R function `duplicated`.
#'
#' @return a boolean vector of the same length as `x`
#' @export
duplicated_tol <- function(x, tol = 1e-6, incl_first = FALSE) {
  if (!incl_first) {
    return(c(FALSE, diff(x) < tol))
  } else {
    diff_lt_tol <- diff(x) < tol
    return(c(FALSE, diff_lt_tol) | c(diff_lt_tol, FALSE))
  }
}


#' Get unique values in a sorted numeric vector, where comparison is up to a
#' specified numeric tolerance. If there is a run of values where each
#' consecutive pair is closer together than the tolerance, all are labeled as
#' corresponding to a single unique value even if not all values in the run are
#' within the tolerance.
#'
#' @param x a numeric vector in which to identify duplicates
#' @param tol numeric tolerance for identifying duplicates
#' @param ties a function that is used to summarize groups of values that fall
#'   within the tolerance
#'
#' @return a numeric vector of the unique values in `x`
#' @export
unique_tol <- function(x, tol = 1e-6, ties = mean) {
  dups_F <- duplicated_tol(x, tol = tol, incl_first = FALSE)
  dups_T <- duplicated_tol(x, tol = tol, incl_first = TRUE)
  if (any(dups_F)) {
    c(dup_run_starts, dup_run_ends) %<-% get_dup_run_inds(dups_F)
    dup_summaries <- purrr::map2_dbl(
      dup_run_starts, dup_run_ends,
      function(start, end) {
        ties(x[start:end])
      }
    )
    ux <- sort(c(x[!dups_T], dup_summaries))
    return(ux)
  } else {
    return(x)
  }
}


#' Get indices of starts and ends of runs of duplicate values
#'
#' @param dups a boolean vector that would result from calling
#'   `duplicated_tol(..., incl_first = FALSE)`
#'
#' @return named list with entries `starts` giving indices of the first element
#'   in each sequence of runs of duplicate values and `ends` giving indices of
#'   the last element in each sequence of runs of duplicate values.
get_dup_run_inds <- function(dups) {
  if (any(dups)) {
    dups_rle <- rle(dups)
    num_runs <- length(dups_rle$lengths)
    dup_run_starts <- cumsum(c(0, dups_rle$lengths[seq_len(num_runs - 1)]))[dups_rle$values]
    dup_run_ends <- cumsum(dups_rle$lengths)[dups_rle$values]
    return(list(starts = dup_run_starts, ends = dup_run_ends))
  } else {
    return(list(starts = integer(), ends = integer()))
  }
}


#' Split ps and qs into those corresponding to discrete and continuous
#' parts of a distribution.
#'
#' @param ps vector of probability levels
#' @param qs vector of quantile values correponding to ps
#' @param dup_tol numeric tolerance for identifying duplicated values indicating
#'   a discrete component of the distribution. If there is a run of values where
#'   each consecutive pair is closer together than the tolerance, all are
#'   labeled as duplicates even if not all values in the run are within the
#'   tolerance.
#' @param zero_tol numeric tolerance for identifying values in `qs` that are
#'   (approximately) zero.
#' @param is_hurdle boolean indicating whether or not this is a hurdle model.
#'   If so, qs of zero always indicate the presence of a point mass at 0.
#'   In this case, 0 is not included among the returned `cont_qs`. Setting this
#'   argument to `TRUE` is primarily appropriate when we are working with a
#'   distributional family that is bounded above 0 (and may have density 0 at 0)
#'   such as a lognormal.
#'
#' @return named list with the following entries:
#'   - `disc_weight`: estimated numeric weight of the discrete part of the
#'     distribution.
#'   - `disc_ps`: estimated probabilities of discrete components. May be
#'     `numeric(0)` if there are no estimated discrete components.
#'   - `disc_qs`: locations of discrete components, corresponding to duplicated
#'     values in the input `qs`. May be `numeric(0)` if there are no discrete
#'     components.
#'   - `cont_ps`: probability levels for the continous part of the distribution
#'   - `cont_qs`: quantile values for the continous part of the distribution
#'   - `disc_ps_range`: a list of length equal to the number of point masses in
#'     the discrete distribution. Each entry is a numeric vector of length two
#'     with the value of the cdf approaching the point mass from the left and
#'     from the right.
#'
#' @export
split_disc_cont_ps_qs <- function(ps, qs, dup_tol = 1e-6, zero_tol = 1e-12,
                                  is_hurdle = FALSE) {
  # if zero counts as discrete and any qs are (approximately) zero, augment
  # with an additional zero so that logic below based on duplicate qs works
  if (is_hurdle && any(abs(qs) <= zero_tol)) {
    zero_ind <- min(which(abs(qs) <= zero_tol))
    if (zero_ind == 1) {
      ps <- c(0.0, ps)
    } else {
      ps <- sort(c(mean(ps[(zero_ind - 1):zero_ind]), ps))
    }
    qs <- sort(c(0.0, qs))
  }

  # Isolate the discrete portion of the distribution:
  # duplicated quantiles and the associated point mass probabilities
  # on the interior, point mass probability is the range of ps with given q
  # on the lower (upper) tail, point mass probability extends to 0 (1)
  dup_q_inds_t <- duplicated_tol(qs, tol = dup_tol, incl_first = TRUE)
  dup_q_inds_f <- duplicated_tol(qs, tol = dup_tol, incl_first = FALSE)

  # if first or last q is duplicated, ensure we get to p = 0 and p = 1
  recalc_dup_inds <- FALSE
  if (dup_q_inds_t[1]) {
    ps <- c(0.0, ps)
    qs <- c(qs[1], qs)
    recalc_dup_inds <- TRUE
  }
  if (tail(dup_q_inds_t, 1)) {
    ps <- c(ps, 1.0)
    qs <- c(qs, tail(qs, 1))
    recalc_dup_inds <- TRUE
  }
  if (recalc_dup_inds) {
    dup_q_inds_t <- duplicated_tol(qs, tol = dup_tol, incl_first = TRUE)
    dup_q_inds_f <- duplicated_tol(qs, tol = dup_tol, incl_first = FALSE)
  }

  disc_qs <- unique_tol(qs[dup_q_inds_t], tol = dup_tol)
  c(dup_run_starts, dup_run_ends) %<-% get_dup_run_inds(dup_q_inds_f)
  disc_ps_range <- purrr::map2(
    dup_run_starts, dup_run_ends,
    function(start, end) {
      range_i <- range(ps[start:end])
      return(range_i)
    }
  )
  disc_ps_mass <- purrr::map_dbl(disc_ps_range, function(range_i) diff(range_i))
  disc_cum_ps <- cumsum(disc_ps_mass)

  # Short-circuit if there is insufficient information from which to estimate
  # a continuous part of the distribution.
  uq <- unique_tol(qs, tol = dup_tol)
  if (length(uq) == 1L) {
    # all qs are duplicates
    # single point mass at that q value
    return(list(
      disc_weight = 1.0,
      disc_ps = 1.0,
      disc_qs = uq,
      cont_ps = numeric(),
      cont_qs = numeric(),
      disc_ps_range = list(c(0.0, 1.0))
    ))
  } else if (length(uq) == 2) {
    # there are only two distinct qs and at least one is at a point mass
    # we declare that both are at point masses
    # assign probabilities based on the tails (to the left of the lower q
    # and to the right of the upper q), and then allocate the mass in the
    # "gap" between the point masses proportionally
    if (length(disc_ps_mass) == 2) {
      disc_ps_unnormalized <- disc_ps_mass
    } else if (length(ps) == 2) {
      disc_ps_unnormalized <- c(ps[1], 1 - ps[2])
    } else if (dup_q_inds_t[1]) {
      disc_ps_unnormalized <- c(disc_ps_mass, 1 - tail(ps, 1))
    } else {
      disc_ps_unnormalized <- c(ps[1], disc_ps_mass)
    }
    disc_ps <- disc_ps_unnormalized / sum(disc_ps_unnormalized)
    return(list(
      disc_weight = 1.0,
      disc_ps = disc_ps,
      disc_qs = uq,
      cont_ps = numeric(),
      cont_qs = numeric(),
      disc_ps_range = list(c(0.0, disc_ps[1]), c(disc_ps[1], 1.0))
    ))
  }

  # remaining quantiles correspond to a continuous portion of the
  # distribution; extract those ps and qs, and adjust the ps,
  # removing any jumps due to point masses
  # Note that we do keep the first instance of a duplicated q.
  # That means that fits for the continuous portion of a distribution
  # will see one (q, p) pair for the duplicated q
  #
  # An exception is qs of 0 when is_hurdle == TRUE, in which case we drop
  # zero from the cont_qs; e.g. including (p=0, q=0) for a lognormal
  # is not informative.
  cont_ps_unadj <- ps[!dup_q_inds_f]
  cont_qs <- uq
  if (is_hurdle) {
    nonzero_inds <- (abs(cont_qs) >= zero_tol)
    cont_ps_unadj <- cont_ps_unadj[nonzero_inds]
    cont_qs <- cont_qs[nonzero_inds]
  }
  cont_ps <- cont_ps_unadj
  for (i in seq_along(disc_qs)) {
    adj_inds <- (cont_qs > disc_qs[i])
    cont_ps[adj_inds] <- cont_ps[adj_inds] - disc_ps_mass[i]
  }

  # adjust for the weight going to the discrete portion of the distribution
  if (length(disc_cum_ps) > 0) {
    disc_weight <- tail(disc_cum_ps, 1)
    disc_ps_mass <- disc_ps_mass / disc_weight
    cont_ps <- cont_ps / (1 - disc_weight)

    # address potential floating point errors: for example,
    # dividing by (1 - disc_weight) could take an input p = 1 to a value
    # slightly greater than 1. here, we ensure that any input p's that
    # were in [0, 1] are still in [0, 1] after adjustment
    inds_nonneg <- cont_ps_unadj >= 0
    cont_ps[inds_nonneg] <- pmax(cont_ps[inds_nonneg], 0)
    inds_lte_1 <- cont_ps_unadj <= 1
    cont_ps[inds_lte_1] <- pmin(cont_ps[inds_lte_1], 1)
  } else {
    disc_weight <- 0.0
  }

  return(list(
    disc_weight = disc_weight,
    disc_ps = disc_ps_mass,
    disc_qs = disc_qs,
    cont_ps = cont_ps,
    cont_qs = cont_qs,
    disc_ps_range = disc_ps_range
  ))
}

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

#' Modify qs by replacing a leading zero with
#' `min(lnorm_zero_buffer, qs[2] / 2)`
#'
#' @param qs numeric vector of quantile values
#' @param lnorm_zero_buffer boolean or numeric specifying how to handle zero
#'   quantiles when `dist = "lnorm"`. If `FALSE`, an error is raised if the
#'   first element of `qs` is zero. Otherwise, must be a positive numeric value,
#'   and the first element of `qs` is replaced by
#'   `min(lnorm_zero_buffer, qs[2] / 2)`.
#'
#' @return updated qs
apply_lnorm_zero_buffer <- function(qs, lnorm_zero_buffer) {
    if (!lnorm_zero_buffer) {
        stop("For dist = 'lnorm', require two positive quantiles")
    } else {
        lnorm_zero_buffer <- as.numeric(lnorm_zero_buffer)
        if (lnorm_zero_buffer <= 0.0) {
            stop("lnorm_zero_buffer must be a boolean or a non-negative numeric value")
        }
        qs[1] <- min(lnorm_zero_buffer, qs[2] / 2)
    }

    return(qs)
}

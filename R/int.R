#' Split ps and qs into those corresponding to discrete and continuous
#' parts of a distribution.
#' 
#' @param ps vector of probability levels
#' @param qs vector of quantile values correponding to ps
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
split_disc_cont_ps_qs <- function(ps, qs) {
    # Short-circuit if all qs are duplicates
    # we have no information from which to estimate a continuous
    # part of the distribution.
    if (length(unique(qs)) == 1L) {
        return(list(
            disc_weight = 1.0,
            disc_ps = 1.0,
            disc_qs = qs[1L],
            cont_ps = numeric(),
            cont_qs = numeric()
        ))
    }

    # Isolate the discrete portion of the distribution:
    # duplicated quantiles and the associated point mass probabilities
    dup_q_inds <- duplicated(qs)
    dup_qs <- qs[dup_q_inds]
    disc_qs <- sort(unique(dup_qs))
    disc_ps <- purrr::map_dbl(disc_qs, function(q) diff(range(ps[qs == q])))
    disc_cum_ps <- cumsum(disc_ps)

    # remaining quantiles correspond to a continuous portion of the
    # distribution; extract those ps and qs, and adjust the ps,
    # removing any jumps due to point masses
    # note that we do keep the first instance of a duplicated q
    # that means that fits for the continuous portion of a distribution
    # will see one (q, p) pair for the duplicated q
    cont_ps <- ps[!dup_q_inds]
    cont_qs <- qs[!dup_q_inds]
    for (i in seq_along(disc_qs)) {
        adj_inds <- (cont_qs > disc_qs[i])
        cont_ps[adj_inds] <- cont_ps[adj_inds] - disc_ps[i]
    }

    # adjust for the weight going to the discrete portion of the distribution
    if (length(disc_cum_ps) > 0) {
        disc_weight <- tail(disc_cum_ps, 1)
        disc_ps <- disc_ps / disc_weight
        cont_ps <- cont_ps / (1 - disc_weight)
    } else {
        disc_weight <- 0.0
    }

    return(list(
        disc_weight = disc_weight,
        disc_ps = disc_ps,
        disc_qs = disc_qs,
        cont_ps = cont_ps,
        cont_qs = cont_qs
    ))
}


#' Approximate density function, cdf, or quantile function on the interior of
#' provided quantiles by representing the distribution as a sum of a discrete
#' part at any duplicated `qs` and a continuous part. The cdf of the continuous
#' part is estimated using a monotonic spline that interpolates the quantiles.
#' To obtain a pdf, we differentiate this spline. To obtain a quantile function,
#' we invert the approximation to the cdf, tracking any discontinuities.
#' 
#' @param ps vector of probability levels
#' @param qs vector of quantile values correponding to ps
#' @param fn_type the type of function that is requested: `"d"` for a pdf,
#'   `"p"` for a cdf, or `"q"` for a quantile function.
#' @param spline_type the type of monotonic spline to fit: "hyman" or "monH.FC"
#' 
#' @return a function to evaluate the pdf, cdf, or quantile function.
spline_cdf <- function(ps, qs, interior_method = c("hyman", "monoH.FC")) {
    
    

    # fit a monotonic spline to the qs and ps to approximate the distribution
    # on the interior
    if (interior_method %in% c("hyman", "monoH.FC")) {
        interior_cdf_spline <- stats::splinefun(qs, ps, method = interior_method)
        interior_pdf <- function(x, log = FALSE) {
            result <- interior_cdf_spline(x, deriv = 1)
        }
    } else {
        stop("Unsupported `interior_method`")
    }

}

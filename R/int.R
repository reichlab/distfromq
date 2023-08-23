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
#' @param zero_discrete boolean indicating whether qs of zero should always
#'   indicate the presence of a point mass at 0. If so, 0 is not included among
#'   the returned `cont_qs`. Primarily appropriate when we
#'   are working with a distributional family that is bounded above 0 (and may
#'   have density 0 at 0, as with a lognormal).
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
                                  zero_discrete = FALSE) {
    # if zero counts as discrete and any qs are (approximately) zero, augment
    # with an additional zero so that logic below based on duplicate qs works
    if (zero_discrete && any(abs(qs) <= zero_tol)) {
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
    dup_q_inds_t <- duplicated_tol(qs, tol = dup_tol, incl_first = TRUE)
    dup_q_inds_f <- duplicated_tol(qs, tol = dup_tol, incl_first = FALSE)
    disc_qs <- unique_tol(qs[dup_q_inds_t], tol = dup_tol)
    c(dup_run_starts, dup_run_ends) %<-% get_dup_run_inds(dup_q_inds_f)
    disc_ps_range <- purrr::map2(
        dup_run_starts, dup_run_ends,
        function(start, end) range(ps[start:end]))
    disc_ps_mass <- purrr::map_dbl(
        disc_ps_range,
        function(range_i) diff(range_i))
    disc_cum_ps <- cumsum(disc_ps_mass)

    # Short-circuit if we have no information from which to estimate a
    # continuous part of the distribution.
    uq <- unique_tol(qs, tol = dup_tol)
    if (length(uq) == 1L) {
        # all qs are duplicates
        return(list(
            disc_weight = 1.0,
            disc_ps = 1.0,
            disc_qs = uq,
            cont_ps = numeric(),
            cont_qs = numeric(),
            disc_ps_range = list(c(0.0, 1.0))
        ))
    } else if (zero_discrete && length(uq) == 2 && any(abs(uq) < zero_tol)) {
        # zero is discrete and there is only one non-zero value
        # allocate the mass in the "gap" between the point masses proportionally
        disc_ps_unnormalized <- disc_ps_mass + c(min(ps), 1 - max(ps))
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
    cont_ps <- ps[!dup_q_inds_f]
    cont_qs <- uq
    if (zero_discrete) {
        nonzero_inds <- (abs(cont_qs) >= zero_tol)
        cont_ps <- cont_ps[nonzero_inds]
        cont_qs <- cont_qs[nonzero_inds]
    }
    for (i in seq_along(disc_qs)) {
        adj_inds <- (cont_qs > disc_qs[i])
        cont_ps[adj_inds] <- cont_ps[adj_inds] - disc_ps_mass[i]
    }

    # adjust for the weight going to the discrete portion of the distribution
    if (length(disc_cum_ps) > 0) {
        disc_weight <- tail(disc_cum_ps, 1)
        disc_ps_mass <- disc_ps_mass / disc_weight
        cont_ps <- cont_ps / (1 - disc_weight)
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

#' Create a polySpline object representing a monotonic Hermite spline
#' interpolating a given set of points.
#'
#' @param x, y: vectors giving the coordinates of the points to be
#'       interpolated.  Alternatively a single plotting structure can
#'       be specified: see ‘xy.coords’.
#'
#'       ‘y’ must be increasing or decreasing for ‘method = "hyman"’.
#'
#'    m: (for ‘splinefunH()’): vector of _slopes_ m[i] at the points
#'       (x[i],y[i]); these together determine the *H*ermite “spline”
#'       which is piecewise cubic, (only) _once_ differentiable
#'       continuously.
#'
#' @details This function essentially reproduces `stats::splinefunH`, but it
#'   returns a polynomial spline object as used in the `splines` package rather
#'   than a function that evaluates the spline, and potentially makes
#'   adjustments to the input slopes `m` to enforce monotonicity.
#'
#' @return An object of class `polySpline` with the spline object, suitable for
#'   use with other functionality from the `splines` package.
mono_Hermite_spline <- function(x, y, m) {
    # For a degree 3 Hermite polynomial, let d be an interpolation parameter
    # between 0 and delta[i] where delta[i] = x[i+1] - x[i] representing the
    # position of x between x[i] and x[i+1], and t = d/delta[i].
    # g(d) = (2 t^3 - 3 t^2 + 1) y[i] + (t^3 - 2 t^2 + t) m[i]
    #         + (-2 t^3 + 3 t^2) y[i+1] + (t^3 - t^2) m[i+1]
    # Note that g(0) = y[i], g(delta[i]) = y[i+1],
    # g'(0) = m[i] / delta_i, and g'(delta[i]) = d/dt ... * dt/dd = m[i+1] / delta_i.
    # We therefore need to adjust the slopes by the factor delta_i.
    #
    # Collecting like terms, we arrive at
    # g(d) = y[i] + m[i] t + (-3 y[i] - 2 m[i] + 3 y[i+1] - m[i+1]) t^2
    #        + (2 y[i] + m[i] - 2 y[i+1] + m[i+1]) t^3
    #      = y[i] + m[i] d / delta[i] 
    #        + (-3 y[i] - 2 m[i] + 3 y[i+1] - m[i+1]) (d / delta[i])^2
    #        + (2 y[i] + m[i] - 2 y[i+1] + m[i+1]) (d / delta[i])^3
    #      = c0 + c1 * d + c2 * d^2 + c3 * d^3, where
    #
    # c0 = y[i],
    # c1 = m[i] / delta[i],
    # c2 = (-3 y[i] - 2 m[i] + 3 y[i+1] - m[i+1]) / delta[i]^2, and
    # c3 = (2 y[i] + m[i] - 2 y[i+1] + m[i+1]) / delta[i]^3
    n <- length(y)
    x_i <- x[-n]
    x_ip1 <- x[-1]
    y_i <- y[-n]
    y_ip1 <- y[-1]
    delta <- x_ip1 - x_i
    delta2 <- delta^2
    delta3 <- delta^3
    
    # adjustments to m to ensure monotonicity; see steps 3 through 5 at
    # https://en.wikipedia.org/wiki/Monotone_cubic_interpolation#Monotone_cubic_Hermite_interpolation
    for (i in seq_along(x_i)) {
        if (y_i[i] == y_ip1[i]) {
            m[i] <- m[i+1] <- 0.0
        }
    }
    for (i in seq_along(x_i)) {
        if (y_i[i] != y_ip1[i]) {
            d <- (y_ip1[i] - y_i[i]) / (x_ip1[i] - x_i[i])
            alpha <- m[i] / d
            beta <- m[i+1] / d
            ab_sq <- alpha^2 + beta^2
            if (ab_sq > 9) {
                tau <- 3 / sqrt(ab_sq)
                m[i] <- tau * m[i]
                m[i+1] <- tau * m[i+1]
            }
        }
    }

    m_i <- m[-n]
    m_ip1 <- m[-1]
    m_i <- m_i * delta
    m_ip1 <- m_ip1 * delta

    ct_0 <- y_i
    ct_1 <- m_i
    ct_2 <- (-3 * y_i - 2 * m_i + 3 * y_ip1 - m_ip1)
    ct_3 <- (2 * y_i + m_i - 2 * y_ip1 + m_ip1)
    spl <- structure(
        list(
            knots = x,
            coefficients = rbind(
                cbind(ct_0, ct_1 / delta, ct_2 / delta2, ct_3 / delta3),
                matrix(c(tail(y, 1), tail(m, 1), 0, 0), nrow = 1)
            )
        ),
        class = c("npolySpline", "polySpline", "spline"))
    return(spl)
}


#' A factory that returns a function that performs linear interpolation,
#' allowing for "steps" or discontinuities.
#'
#' @param x: numeric vector with the "horizontal axis" coordinates of the points
#'   to interpolate.
#' @param y: numeric vector with the "vertical axis" coordinates of the points
#'   to interpolate.
#' @param cont_dir: at steps or discontinuities, the direction from which the
#'   function is continuous. This will be "right" for a cdf or "left" for a qf.
#' @param increasing: boolean indicating whether the function is increasing or
#'   decreasing. Only used in the degenerate case where there is only one unique
#'   value of `x`.
#'
#' @return a function with argument `x` that performs linear approximation of
#'   the input data points.
step_interp_factory <- function(x, y, cont_dir = c("right", "left"),
                                increasing = TRUE) {
    if ((length(x) != length(y)) || length(x) == 0) {
        stop("x and y must have the same non-zero length")
    }
    cont_dir <- match.arg(cont_dir)

    dup_inds <- duplicated(x)
    if (sum(dup_inds) == 0) {
        dup_x <- NULL
    } else {
        dup_x <- sort(x[dup_inds])
    }
    interval_endpoints <- unique(c(min(x), dup_x, max(x)))
    n_intervals <- length(interval_endpoints) - 1

    interp_funs <- purrr::map(
        seq_len(n_intervals),
        function(i) {
            l <- interval_endpoints[i]
            u <- interval_endpoints[i + 1]

            l_ind <- which(x == l)
            l_ind <- l_ind[which.max(y[l_ind])]
            u_ind <- which(x == u)
            u_ind <- u_ind[which.min(y[u_ind])]
            int_inds <- which((x > l) & (x < u))
            inds <- c(l_ind, int_inds, u_ind)

            return(approxfun(x[inds], y[inds]))
        }
    )

    # overall left/right endpoints and comparator functions for determining
    # whether x values are out of bounds on each side, and the prediction
    # at the endpoint that is not covered by half-open intervals
    l <- interval_endpoints[1]
    u <- tail(interval_endpoints, 1)
    if (cont_dir == "right") {
        l_comp <- `<`
        r_comp <- `>=`
        # approaching the right endpoint from the right in an increasing
        # function, we should predict the max of the corresponding y's;
        # in a decreasing function, predict the min of the y's
        extr_endpoint <- u
        extr_fn <- if (increasing) max else min
        extr_y <- extr_fn(y[x == extr_endpoint])
    } else {
        l_comp <- `<=`
        r_comp <- `>`
        # approaching the left endpoint from the left in an increasing
        # function, we should predict the min of the corresponding y's;
        # in a decreasing function, predict the max of the y's
        extr_endpoint <- l
        extr_fn <- if (increasing) min else max
        extr_y <- extr_fn(y[x == extr_endpoint])
    }

    # function of x that performs linear interpolation with steps
    step_interp <- function(x) {
        result <- rep(NA_real_, length(x))

        # prediction at endpoint not covered by half-open intervals
        extr_inds <- (x == extr_endpoint)
        result[extr_inds] <- extr_y

        # predictions by linear interpolation between jump points
        if (n_intervals > 0) {
            # identify which interval each x observation falls into (if any)
            interval_ind <- rep(-1L, length(x))
            ib_inds <- !(l_comp(x, l) | r_comp(x, u))
            interval_ind[ib_inds] <- cut(x[ib_inds],
                                         breaks = interval_endpoints,
                                         labels = FALSE,
                                         right = (cont_dir != "right"))
            for (i in unique(interval_ind)) {
                if (i > 0) {
                    x_inds <- (interval_ind == i)
                    result[x_inds] <- interp_funs[[i]](x[x_inds])
                }
            }
        }

        return(result)
    }

    return(step_interp)
}


#' Approximate density function, cdf, or quantile function on the interior of
#' provided quantiles by representing the distribution as a sum of a discrete
#' part at any duplicated `qs` and a continuous part for which the cdf is
#' estimated using a monotonic Hermite spline. See details for more.
#'
#' @param ps vector of probability levels
#' @param qs vector of quantile values correponding to ps
#' @param lower_tail_dist name of parametric distribution for the lower tail
#' @param upper_tail_dist name of parametric distribution for the upper tail
#' @param fn_type the type of function that is requested: `"d"` for a pdf,
#'   `"p"` for a cdf, or `"q"` for a quantile function.
#' @param n_grid grid size to use when augmenting the input `qs` to obtain a
#'   finer grid of points along which we form a piecewise linear approximation
#'   to the spline. `n_grid` evenly-spaced points are inserted between each
#'   pair of consecutive values in `qs`. The default value is 20. This can
#'   be set to `NULL`, in which case the piecewise linear approximation is not
#'   used. This is not recommended if the `fn_type` is `"q"`.
#' @param dup_tol numeric tolerance for identifying duplicated quantiles
#'   indicating a discrete component of the distribution. If there is a run of
#'   values where each consecutive pair is closer together than the tolerance,
#'   all are labeled as duplicates even if not all values in the run are within
#'   the tolerance.
#' @param zero_tol numeric tolerance for identifying whether elements of `qs`
#'   are equal to zero. This is used only for identifying point masses at zero
#'   if `lower_tail_dist` or `upper_tail_dist` is "lnorm".
#'
#' @details The cdf of the continuous part of the distribution is estimated
#' using a monotonic degree 3 Hermite spline that interpolates the quantiles
#' after subtracting the discrete distribution and renormalizing. In theory,
#' an estimate of the quantile function could be obtained by directly inverting
#' this spline. However, in practice, we have observed that this can suffer from
#' numerical problems. Therefore, the default behavior of this function is to
#' evaluate the "stage 1" cdf estimate corresponding to discrete point masses
#' plus monotonic spline at a fine grid of points, and use the "stage 2" cdf
#' estimate that linearly interpolates these points with steps at any duplicated
#' q values. The quantile function estimate is obtained by inverting this
#' "stage 2" cdf estimate. When the distribution is continuous, we can obtain an
#' estimate of the pdf by differentiating the cdf estimate, resulting in a
#' discontinuous "histogram density". The size of the grid can be specified with
#' the `n_grid` argument. In settings where it is desirable to obtain a
#' continuous density function, the "stage 1" cdf estimate can be used by
#' setting `n_grid = NULL`.
#'
#' @return a function to evaluate the pdf, cdf, or quantile function.
spline_cdf <- function(ps, qs, lower_tail_dist, upper_tail_dist,
                       fn_type = c("d", "p", "q"),
                       n_grid = 20,
                       dup_tol = 1e-6,
                       zero_tol = 1e-12) {
    fn_type <- match.arg(fn_type)

    # TODO: use duplicated_tol and unique_tol throughout spline_cdf
    if ((any(duplicated_tol(qs, tol = dup_tol)) || length(qs) == 1L) &&
            fn_type == "d") {
        stop("Distribution has a discrete component;",
             " cannot create a density function.")
    }

    # if either lower_tail_dist or upper_tail_dist is "lnorm", we treat
    # quantiles of 0 as indicative of a discrete point mass at 0
    zero_discrete <- lower_tail_dist == "lnorm" || upper_tail_dist == "lnorm"

    if (!is.null(n_grid)) {
        uq <- unique_tol(qs, tol = dup_tol)
        if (length(uq) > 1L) {
            # augment originally provided ps and qs with new grid of qs

            n_grid <- as.integer(n_grid)
            if (n_grid < 1) {
                stop("If provided, `n_grid` must be a positive integer.")
            }

            # by recursing but with n_grid = NULL, get a function that evaluates
            # the approximated cdf based on a spline + discrete steps
            p_fn <- spline_cdf(ps = ps, qs = qs,
                               lower_tail_dist = lower_tail_dist,
                               upper_tail_dist = upper_tail_dist,
                               fn_type = "p",
                               n_grid = NULL,
                               dup_tol = dup_tol,
                               zero_tol = zero_tol)

            # obtain a grid of intermediate points (between the provided qs) at
            # which we will evaluate the spline approximation
            q_grid <- purrr::map(
                seq_len(length(uq) - 1),
                function(i) {
                    new_q <- seq(from = uq[i], to = uq[i + 1],
                                 length = n_grid + 2L)
                    return(new_q[2:(n_grid + 1)])
                }) %>%
                unlist()

            # evaluate spline-based cdf approximation at new q values, and
            # append to the inputs.
            p_grid <- p_fn(q_grid)

            ps <- sort(c(ps, p_grid))
            qs <- sort(c(qs, q_grid))
        }

        # Create resulting function and return
        if (fn_type == "d") {
            # note that we have already validated that qs are all distinct and
            # there are at least two qs
            # check that we don't have a point mass at 0
            if (zero_discrete && any(abs(qs) < zero_tol)) {
                stop("Distribution has a discrete component;",
                     " cannot create a density function.")
            }

            slopes <- diff(ps) / diff(qs)
            slopes <- c(slopes, tail(slopes, 1))
            int_d_fn <- function(x, log = FALSE) {
                result <- approx(x = qs, y = slopes, xout = x,
                                 method = "constant", f = 0.0)$y
                if (log) return(log(result)) else return(result)
            }
            return(int_d_fn)
        } else if (fn_type == "p") {
            step_interp <- step_interp_factory(x = qs, y = ps,
                                               cont_dir = "right")
            int_p_fn <- function(q, log.p = FALSE) {
                result <- step_interp(q)
                if (log.p) return(log(result)) else return(result)
            }
            return(int_p_fn)
        } else if (fn_type == "q") {
            step_interp <- step_interp_factory(x = ps, y = qs,
                                                cont_dir = "left")
            int_q_fn <- function(p) {
                return(step_interp(p))
            }
            return(int_q_fn)
        }
    } else {
        # split ps and qs into discrete and continuous part, along with the
        # weight given to the discrete part
        c(disc_weight, disc_ps, disc_qs, cont_ps, cont_qs, disc_ps_range) %<-%
            split_disc_cont_ps_qs(ps, qs, dup_tol = dup_tol,
                                  zero_tol = zero_tol,
                                  zero_discrete = zero_discrete)

        # fit a monotonic spline to the qs and ps for the continuous part of the
        # distribution to approximate the cdf on the interior
        # on ends, slope of cdf approximation should match tail distribution pdf
        # on interior, slope is the mean of the slopes of the adjacent segments
        if (disc_weight < 1.0) {
            d_lower <- d_ext_factory(ps = head(cont_ps, 2),
                                     qs = head(cont_qs, 2),
                                     dist = lower_tail_dist,
                                     dup_tol, zero_tol, zero_discrete, TRUE)
            m_lower <- d_lower(cont_qs[1])
            d_upper <- d_ext_factory(ps = tail(cont_ps, 2),
                                     qs = tail(cont_qs, 2),
                                     dist = upper_tail_dist,
                                     dup_tol, zero_tol, zero_discrete, FALSE)
            m_upper <- d_upper(tail(cont_qs, 1))

            m_segments <- diff(cont_ps) / diff(cont_qs)
            n <- length(m_segments)
            m_interior <- apply(cbind(m_segments[-1], m_segments[-n]), 1, mean)
            m <- c(m_lower, m_interior, m_upper)

            interior_cdf_spline <- mono_Hermite_spline(x = cont_qs, y = cont_ps,
                                                       m = m)
        }

        # get a function that calculates the pdf, cdf, or quantile function
        if (fn_type == "d") {
            if (disc_weight > 0) {
                stop("Distribution has a discrete component;",
                    " cannot create a density function.")
            }
            int_d_fn <- function(x, log = FALSE) {
                result <- predict(interior_cdf_spline, x, deriv = 1)$y
                if (log) {
                    return(log(result))
                } else {
                    return(result)
                }
            }
            return(int_d_fn)
        } else if (fn_type == "p") {
            int_p_fn <- function(q, log.p = FALSE) {
                if (disc_weight < 1.0) {
                    result <- predict(interior_cdf_spline, q, deriv = 0)$y
                    result <- result * (1 - disc_weight)
                } else {
                    result <- rep(0.0, length(q))
                }

                for (i in seq_along(disc_ps)) {
                    inds <- (q >= disc_qs[i])
                    result[inds] <- result[inds] + disc_ps[i] * disc_weight
                }

                if (log.p) {
                    return(log(result))
                } else {
                    return(result)
                }
            }
            return(int_p_fn)
        } else if (fn_type == "q") {
            if (disc_weight < 1.0) {
                interior_qf_spline <- backSpline(interior_cdf_spline)
            }
            int_q_fn <- function(p) {
                result <- rep(NA_real_, length(p))
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

                if (disc_weight < 1.0) {
                    cont_p <- cont_p / (1 - disc_weight)
                    result[cont_inds] <- predict(interior_qf_spline,
                                                cont_p[cont_inds])$y
                }
                return(result)
            }
            return(int_q_fn)
        }
    }
}

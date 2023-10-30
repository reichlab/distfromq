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
#' @param tail_dist name of parametric distribution for the tails
#' @param fn_type the type of function that is requested: `"d"` for a pdf,
#'   `"p"` for a cdf, or `"q"` for a quantile function.
#' @param n_grid grid size to use when augmenting the input `qs` to obtain a
#'   finer grid of points along which we form a piecewise linear approximation
#'   to the spline. `n_grid` evenly-spaced points are inserted between each
#'   pair of consecutive values in `qs`. The default value is 20. This can
#'   be set to `NULL`, in which case the piecewise linear approximation is not
#'   used. This is not recommended if the `fn_type` is `"q"`.
#' @param tol numeric tolerance for identifying duplicated quantiles indicating
#'   a discrete component of the distribution. If there is a run of values where
#'   each consecutive pair is closer together than the tolerance, all are
#'   labeled as duplicates even if not all values in the run are within the
#'   tolerance.
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
spline_cdf <- function(ps, qs, tail_dist,
                       fn_type = c("d", "p", "q"),
                       n_grid = 20) {
    fn_type <- match.arg(fn_type)

    if (!is.null(n_grid)) {
        return(spline_cdf_grid_interp(ps, qs, tail_dist, fn_type, n_grid))
    } else {
        return(spline_cdf_direct(ps, qs, tail_dist, fn_type))
    }
}

spline_cdf_grid_interp <- function(ps, qs, tail_dist,
                                   fn_type = c("d", "p", "q"),
                                   n_grid = 20) {
    # augment originally provided ps and qs with new grid of qs
    c(ps, qs) %<-% grid_augment_ps_qs(ps, qs, tail_dist, n_grid)

    # Create resulting function and return
    if (fn_type == "d") {
        # note that upstream, we have already validated that qs are all
        # distinct and there are at least two qs
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
}

#' Internal function that augments ps and qs by filling in a grid of
#' intermediate values for qs. n_grid new points are inserted between each
#' pair of consecutive values in qs, and the corresponding ps are filled in
#' by evaluating the spline output from `spline_cdf_direct` at the qs grid.
grid_augment_ps_qs <- function(ps, qs, tail_dist, n_grid) {
    n_grid <- as.integer(n_grid)
    if (n_grid < 0) {
        stop("If provided, `n_grid` must be a non-negative integer.")
    } else if (n_grid == 0) {
        # if the user just wants to linearly interpolate the given qs, ps
        return(list(ps = ps, qs = qs))
    }

    # get a function that evaluates the approximated cdf based on a spline
    p_fn <- spline_cdf_direct(ps = ps, qs = qs, tail_dist = tail_dist,
                              fn_type = "p")

    # obtain a grid of intermediate points (between the provided qs) at
    # which we will evaluate the spline approximation
    q_grid <- purrr::map(
        seq_len(length(qs) - 1),
        function(i) {
            new_q <- seq(from = qs[i], to = qs[i + 1], length = n_grid + 2L)
            return(new_q[2:(n_grid + 1)])
        })
    q_grid <- unlist(q_grid)

    # evaluate spline-based cdf approximation at new q values, and
    # append to the inputs.
    p_grid <- p_fn(q_grid)

    ps <- sort(c(ps, p_grid))
    qs <- sort(c(qs, q_grid))

    return(list(ps = ps, qs = qs))
}

#' Internal function that constructs a monotonic Hermite spline interpolating
#' ps and qs.
spline_cdf_direct <- function(ps, qs, tail_dist,
                              fn_type = c("d", "p", "q")) {
    # fit a monotonic spline to the qs and ps for the continuous part of the
    # distribution to approximate the cdf on the interior
    # the vector m that we assemble here has the target slopes of the spline
    # at each data point (qs[i], ps[i])
    # on ends, slope of cdf approximation should match tail distribution pdf
    # on interior, slope is the mean of the slopes of the adjacent segments
    # note: there is probably room for improvement for this interior behavior,
    # but this seems to be the standard strategy for monotonic splines
    d_lower <- d_ext_factory(ps = head(ps, 2), qs = head(qs, 2),
                             dist = tail_dist)
    m_lower <- d_lower(qs[1])
    d_upper <- d_ext_factory(ps = tail(ps, 2), qs = tail(qs, 2),
                             dist = tail_dist)
    m_upper <- d_upper(tail(qs, 1))

    m_segments <- diff(ps) / diff(qs)
    n <- length(m_segments)
    m_interior <- apply(cbind(m_segments[-1], m_segments[-n]), 1, mean)

    # in case of errors in evaluating density function at endpoints,
    # repeat interior values
    if (is.nan(m_lower) || !is.finite(m_lower)) {
        if (n == 1) {
            m_lower <- m_segments
        } else {
            m_lower <- m_interior[1]
        }
    }
    if (is.nan(m_upper) || !is.finite(m_upper)) {
        if (n == 1) {
            m_upper <- m_segments
        } else {
            m_upper <- tail(m_interior, 1)
        }
    }

    m <- c(m_lower, m_interior, m_upper)

    interior_cdf_spline <- mono_Hermite_spline(x = qs, y = ps, m = m)

    # get a function that calculates the pdf, cdf, or quantile function
    if (fn_type == "d") {
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
            result <- predict(interior_cdf_spline, q, deriv = 0)$y

            if (log.p) {
                return(log(result))
            } else {
                return(result)
            }
        }
        return(int_p_fn)
    } else if (fn_type == "q") {
        interior_qf_spline <- backSpline(interior_cdf_spline)

        int_q_fn <- function(p) {
            return(predict(interior_qf_spline, p)$y)
        }
        return(int_q_fn)
    }
}

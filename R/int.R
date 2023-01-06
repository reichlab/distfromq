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
            cont_qs = numeric(),
            disc_ps_range = list(range(ps))
        ))
    }

    # Isolate the discrete portion of the distribution:
    # duplicated quantiles and the associated point mass probabilities
    dup_q_inds <- duplicated(qs)
    dup_qs <- qs[dup_q_inds]
    disc_qs <- sort(unique(dup_qs))
    disc_ps_range <- purrr::map(
        disc_qs,
        function(q) range(ps[qs == q]))
    disc_ps_mass <- purrr::map_dbl(
        disc_qs,
        function(q) diff(range(ps[qs == q])))
    disc_cum_ps <- cumsum(disc_ps_mass)

    # remaining quantiles correspond to a continuous portion of the
    # distribution; extract those ps and qs, and adjust the ps,
    # removing any jumps due to point masses
    # Note that we do keep the first instance of a duplicated q.
    # That means that fits for the continuous portion of a distribution
    # will see one (q, p) pair for the duplicated q
    cont_ps <- ps[!dup_q_inds]
    cont_qs <- qs[!dup_q_inds]
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
#'   than a function that evaluates the spline.
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
    # Label these coefficients as c0 = y[i], c1 = m[i] / delta[i],
    # c2 = (-3 y[i] - 2 m[i] + 3 y[i+1] - m[i+1]) / delta[i]^2, and
    # c3 = (2 y[i] + m[i] - 2 y[i+1] + m[i+1]) / delta[i]^3
    n <- length(y)
    x_i <- x[-n]
    x_ip1 <- x[-1]
    y_i <- y[-n]
    y_ip1 <- y[-1]
    m_i <- m[-n]
    m_ip1 <- m[-1]
    delta <- x_ip1 - x_i
    delta2 <- delta^2
    delta3 <- delta^3

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
#' @param spline_method the type of monotonic spline to fit: "hyman" or "monH.FC"
#' 
#' @return a function to evaluate the pdf, cdf, or quantile function.
spline_cdf <- function(ps, qs, fn_type = c("d", "p", "q"),
                       lower_tail_dist,
                       upper_tail_dist) {
    fn_type <- match.arg(fn_type)

    # split ps and qs into discrete and continuous part, along with the weight
    # given to the discrete part
    c(disc_weight, disc_ps, disc_qs, cont_ps, cont_qs, disc_ps_range) %<-%
        split_disc_cont_ps_qs(ps, qs)

    # fit a monotonic spline to the qs and ps for the continuous part of the
    # distribution to approximate the cdf on the interior
    # on ends, slope of cdf approximation should match tail distribution pdf
    # on interior, slope is the mean of the slopes of the adjacent line segments
    if (disc_weight < 1.0) {
        d_lower <- d_ext_factory(ps = head(cont_ps, 2), qs = head(cont_qs, 2),
                                dist = lower_tail_dist)
        m_lower <- d_lower(cont_qs[1])
        d_upper <- d_ext_factory(ps = tail(cont_ps, 2), qs = tail(cont_qs, 2),
                                dist = lower_tail_dist)
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
        int_p_fn <- function(x, log.p = FALSE) {
            if (disc_weight < 1.0) {
                result <- predict(interior_cdf_spline, x, deriv = 0)$y
                result <- result * (1 - disc_weight)
            } else {
                result <- rep(0.0, length(x))
            }

            for (i in seq_along(disc_ps)) {
                inds <- (x >= disc_qs[i])
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

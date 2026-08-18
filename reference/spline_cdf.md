# Approximate density function, CDF, or quantile function on the interior of provided quantiles by representing the distribution as a sum of a discrete part at any duplicated `qs` and a continuous part for which the CDF is estimated using a monotonic Hermite spline. See details for more.

Approximate density function, CDF, or quantile function on the interior
of provided quantiles by representing the distribution as a sum of a
discrete part at any duplicated `qs` and a continuous part for which the
CDF is estimated using a monotonic Hermite spline. See details for more.

## Usage

``` r
spline_cdf(ps, qs, tail_dist, fn_type = c("d", "p", "q"), n_grid = 20)
```

## Arguments

- ps:

  vector of probability levels

- qs:

  vector of quantile values corresponding to ps

- tail_dist:

  name of parametric distribution for the tails

- fn_type:

  the type of function that is requested: `"d"` for a PDF, `"p"` for a
  CDF, or `"q"` for a quantile function.

- n_grid:

  grid size to use when augmenting the input `qs` to obtain a finer grid
  of points along which we form a piecewise linear approximation to the
  spline. `n_grid` evenly-spaced points are inserted between each pair
  of consecutive values in `qs`. The default value is 20. This can be
  set to `NULL`, in which case the piecewise linear approximation is not
  used. This is not recommended if the `fn_type` is `"q"`.

## Value

a function to evaluate the PDF, CDF, or quantile function.

## Details

The CDF of the continuous part of the distribution is estimated using a
monotonic degree 3 Hermite spline that interpolates the quantiles after
subtracting the discrete distribution and renormalizing. In theory, an
estimate of the quantile function could be obtained by directly
inverting this spline. However, in practice, we have observed that this
can suffer from numerical problems. Therefore, the default behavior of
this function is to evaluate the "stage 1" CDF estimate corresponding to
discrete point masses plus monotonic spline at a fine grid of points,
and use the "stage 2" CDF estimate that linearly interpolates these
points with steps at any duplicated q values. The quantile function
estimate is obtained by inverting this "stage 2" CDF estimate. When the
distribution is continuous, we can obtain an estimate of the PDF by
differentiating the CDF estimate, resulting in a discontinuous
"histogram density". The size of the grid can be specified with the
`n_grid` argument. In settings where it is desirable to obtain a
continuous density function, the "stage 1" CDF estimate can be used by
setting `n_grid = NULL`.

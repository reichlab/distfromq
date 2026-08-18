# Creates a function that evaluates the probability density function of an approximation to a distribution obtained by interpolating and extrapolating from a set of quantiles of the distribution.

Creates a function that evaluates the probability density function of an
approximation to a distribution obtained by interpolating and
extrapolating from a set of quantiles of the distribution.

## Usage

``` r
make_d_fn(
  ps,
  qs,
  interior_method = "spline_cdf",
  interior_args = list(),
  tail_dist = "norm",
  dup_tol = 1e-06,
  zero_tol = 1e-12
)
```

## Arguments

- ps:

  vector of probability levels

- qs:

  vector of quantile values corresponding to ps

- interior_method:

  method for interpolating the distribution on the interior of the
  provided `qs`. This package provides one method for this,
  `"spline_cdf"`. The user may also provide a custom function; see the
  details for more.

- interior_args:

  an optional named list of arguments that are passed on to the
  `interior_method`

- tail_dist:

  name of parametric distribution for the tails

- dup_tol:

  numeric tolerance for identifying duplicated values indicating a
  discrete component of the distribution. If there is a run of values
  where each consecutive pair is closer together than the tolerance, all
  are labeled as duplicates even if not all values in the run are within
  the tolerance.

- zero_tol:

  numeric tolerance for identifying values in `qs` that are
  (approximately) zero.

## Value

a function with arguments `x` and `log` that can be used to evaluate the
approximate density function (or its `log`) at the points `x`.

## Details

The default `interior_method`, `"spline_cdf"`, represents the
distribution as a sum of a discrete component at any points where there
are duplicated `qs` for multiple different `ps` and a continuous
component that is estimated by using a monotonic cubic spline that
interpolates the provided `(q, p)` pairs as an estimate of the CDF. The
density function is then obtained by differentiating this estimate of
the CDF.

Optionally, the user may provide another function that accepts arguments
`ps`, `qs`, `tail_dist`, and `fn_type` (which will be either `"d"`,
`"p"`, or `"q"`), and optionally additional named arguments to be
specified via `interior_args`. This function should return a function
with arguments `x`, `log` that evaluates the pdf or its logarithm.

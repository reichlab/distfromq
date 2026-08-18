# A factory that returns a function that performs linear interpolation, allowing for "steps" or discontinuities.

A factory that returns a function that performs linear interpolation,
allowing for "steps" or discontinuities.

## Usage

``` r
step_interp_factory(x, y, cont_dir = c("right", "left"), increasing = TRUE)
```

## Arguments

- x:

  numeric vector with the "horizontal axis" coordinates of the points to
  interpolate.

- y:

  numeric vector with the "vertical axis" coordinates of the points to
  interpolate.

- cont_dir:

  at steps or discontinuities, the direction from which the function is
  continuous. This will be "right" for a CDF or "left" for a QF.

- increasing:

  boolean indicating whether the function is increasing or decreasing.
  Only used in the degenerate case where there is only one unique value
  of `x`.

## Value

a function with argument `x` that performs linear approximation of the
input data points.

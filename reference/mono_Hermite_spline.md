# Create a polySpline object representing a monotonic Hermite spline interpolating a given set of points.

Create a polySpline object representing a monotonic Hermite spline
interpolating a given set of points.

## Usage

``` r
mono_Hermite_spline(x, y, m)
```

## Arguments

- x:

  vector giving the x coordinates of the points to be interpolated.

- y:

  vector giving the y coordinates of the points to be interpolated. Must
  be increasing or decreasing for 'method = "hyman"'.

- m:

  (for 'splinefunH()') vector of *slopes* \\m_i\\ at the points
  \\(x_i,y_i)\\; these together determine the *H*ermite “spline” which
  is piecewise cubic, (only) *once* differentiable continuously.

## Value

An object of class `polySpline` with the spline object, suitable for use
with other functionality from the `splines` package.

## Details

This function essentially reproduces
[`stats::splinefunH`](https://rdrr.io/r/stats/splinefun.html), but it
returns a polynomial spline object as used in the `splines` package
rather than a function that evaluates the spline, and potentially makes
adjustments to the input slopes `m` to enforce monotonicity.

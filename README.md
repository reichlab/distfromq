# distfromq

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/reichlab/distfromq/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/reichlab/distfromq/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

Given a set of predictive quantiles from a distribution, estimate the distribution and create `d`, `p`, `q`, and `r` functions to evaluate its density function, distribution function, and quantile function, and generate random samples. On the interior of the provided quantiles, an interpolation method such as a monotonic cubic spline is used; the tails are approximated by a location-scale family.

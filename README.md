# distfromq

Given a set of predictive quantiles from a distribution, estimate the distribution and create `d`, `p`, `q`, and `r` functions to evaluate its density function, distribution function, and quantile function, and generate random samples. On the interior of the provided quantiles, an interpolation method such as a monotonic cubic spline is used; the tails are approximated by a location-scale family.

## News

### Version 0.1.2

Correct bugs in handling of situations where only one unique `q` is provided as input.

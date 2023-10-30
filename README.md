# distfromq

Given a set of predictive quantiles from a distribution, estimate the distribution and create `d`, `p`, `q`, and `r` functions to evaluate its density function, distribution function, and quantile function, and generate random samples. On the interior of the provided quantiles, an interpolation method such as a monotonic cubic spline is used; the tails are approximated by a location-scale family.

## News

### Version 1.0.2

Changes:
 - Fix a bug that came up with multiple discrete point masses, where floating point issues could cause adjusted probability levels to fall outside the interval [0, 1].
 - Add argument checking.

### Version 1.0.1

Changes:
 - Fix a bug when all quantiles occur at one of three duplicated values, one of which is zero, when the tail distribution is `"lnorm"`.

### Version 1.0.0

Changes:
 - Refactored handling of discrete point masses in the distribution, with a goal of better supporting use of a lognormal distribution in the tails.
 - Removed separate arguments for `lower_tail_dist` and `upper_tail_dist`, replaced by a single `tail_dist` that is used in both tails.
 - Introduced parameters `dup_tol`, which is used for detecting point masses when quantile values are approximate duplicates (i.e., the distance between them is less than `dup_tol`), and `zero_tol`, which is used for detecting point masses at 0 when a lognormal family is used.
 - Updated the vignette
 - Added a more comprehensive suite of unit tests.

We now recommend the use of `tail_dist = "lnorm"` when the target variable is bounded above zero. This choice will automatically handle the non-negativity constraint correctly, and will have a longer right tail than the default choice of `tail_dist = "norm"`.

### Version 0.1.2

Correct bugs in handling of situations where only one unique `q` is provided as input.

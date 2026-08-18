# Split ps and qs into those corresponding to discrete and continuous parts of a distribution.

Split ps and qs into those corresponding to discrete and continuous
parts of a distribution.

## Usage

``` r
split_disc_cont_ps_qs(
  ps,
  qs,
  dup_tol = 1e-06,
  zero_tol = 1e-12,
  is_hurdle = FALSE
)
```

## Arguments

- ps:

  vector of probability levels

- qs:

  vector of quantile values corresponding to ps

- dup_tol:

  numeric tolerance for identifying duplicated values indicating a
  discrete component of the distribution. If there is a run of values
  where each consecutive pair is closer together than the tolerance, all
  are labeled as duplicates even if not all values in the run are within
  the tolerance.

- zero_tol:

  numeric tolerance for identifying values in `qs` that are
  (approximately) zero.

- is_hurdle:

  boolean indicating whether or not this is a hurdle model. If so, qs of
  zero always indicate the presence of a point mass at 0. In this case,
  0 is not included among the returned `cont_qs`. Setting this argument
  to `TRUE` is primarily appropriate when we are working with a
  distributional family that is bounded above 0 (and may have density 0
  at 0) such as a lognormal.

## Value

named list with the following entries:

- `disc_weight`: estimated numeric weight of the discrete part of the
  distribution.

- `disc_ps`: estimated probabilities of discrete components. May be
  `numeric(0)` if there are no estimated discrete components.

- `disc_qs`: locations of discrete components, corresponding to
  duplicated values in the input `qs`. May be `numeric(0)` if there are
  no discrete components.

- `cont_ps`: probability levels for the continuous part of the
  distribution

- `cont_qs`: quantile values for the continuous part of the distribution

- `disc_ps_range`: a list of length equal to the number of point masses
  in the discrete distribution. Each entry is a numeric vector of length
  two with the value of the CDF approaching the point mass from the left
  and from the right.

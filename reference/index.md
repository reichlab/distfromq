# Package index

## All functions

- [`duplicated_tol()`](http://reichlab.io/distfromq/reference/duplicated_tol.md)
  : Identify duplicated values in a sorted numeric vector, where
  comparison is up to a specified numeric tolerance. If there is a run
  of values where each consecutive pair is closer together than the
  tolerance, all are labeled as duplicates even if not all values in the
  run are within the tolerance.

- [`get_dup_run_inds()`](http://reichlab.io/distfromq/reference/get_dup_run_inds.md)
  : Get indices of starts and ends of runs of duplicate values

- [`make_d_fn()`](http://reichlab.io/distfromq/reference/make_d_fn.md) :
  Creates a function that evaluates the probability density function of
  an approximation to a distribution obtained by interpolating and
  extrapolating from a set of quantiles of the distribution.

- [`make_p_fn()`](http://reichlab.io/distfromq/reference/make_p_fn.md) :
  Creates a function that evaluates the cumulative distribution function
  of an approximation to a distribution obtained by interpolating and
  extrapolating from a set of quantiles of the distribution.

- [`make_q_fn()`](http://reichlab.io/distfromq/reference/make_q_fn.md) :
  Creates a function that evaluates the quantile function of an
  approximation to a distribution obtained by interpolating and
  extrapolating from a set of quantiles of the distribution.

- [`make_r_fn()`](http://reichlab.io/distfromq/reference/make_r_fn.md) :
  Creates a function that generates random deviates from an
  approximation to a distribution obtained by interpolating and
  extrapolating from a set of quantiles of the distribution.

- [`mono_Hermite_spline()`](http://reichlab.io/distfromq/reference/mono_Hermite_spline.md)
  : Create a polySpline object representing a monotonic Hermite spline
  interpolating a given set of points.

- [`spline_cdf()`](http://reichlab.io/distfromq/reference/spline_cdf.md)
  :

  Approximate density function, CDF, or quantile function on the
  interior of provided quantiles by representing the distribution as a
  sum of a discrete part at any duplicated `qs` and a continuous part
  for which the CDF is estimated using a monotonic Hermite spline. See
  details for more.

- [`split_disc_cont_ps_qs()`](http://reichlab.io/distfromq/reference/split_disc_cont_ps_qs.md)
  : Split ps and qs into those corresponding to discrete and continuous
  parts of a distribution.

- [`step_interp_factory()`](http://reichlab.io/distfromq/reference/step_interp_factory.md)
  : A factory that returns a function that performs linear
  interpolation, allowing for "steps" or discontinuities.

- [`unique_tol()`](http://reichlab.io/distfromq/reference/unique_tol.md)
  : Get unique values in a sorted numeric vector, where comparison is up
  to a specified numeric tolerance. If there is a run of values where
  each consecutive pair is closer together than the tolerance, all are
  labeled as corresponding to a single unique value even if not all
  values in the run are within the tolerance.

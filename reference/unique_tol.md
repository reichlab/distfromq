# Get unique values in a sorted numeric vector, where comparison is up to a specified numeric tolerance. If there is a run of values where each consecutive pair is closer together than the tolerance, all are labeled as corresponding to a single unique value even if not all values in the run are within the tolerance.

Get unique values in a sorted numeric vector, where comparison is up to
a specified numeric tolerance. If there is a run of values where each
consecutive pair is closer together than the tolerance, all are labeled
as corresponding to a single unique value even if not all values in the
run are within the tolerance.

## Usage

``` r
unique_tol(x, tol = 1e-06, ties = mean)
```

## Arguments

- x:

  a numeric vector in which to identify duplicates

- tol:

  numeric tolerance for identifying duplicates

- ties:

  a function that is used to summarize groups of values that fall within
  the tolerance

## Value

a numeric vector of the unique values in `x`

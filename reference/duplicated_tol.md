# Identify duplicated values in a sorted numeric vector, where comparison is up to a specified numeric tolerance. If there is a run of values where each consecutive pair is closer together than the tolerance, all are labeled as duplicates even if not all values in the run are within the tolerance.

Identify duplicated values in a sorted numeric vector, where comparison
is up to a specified numeric tolerance. If there is a run of values
where each consecutive pair is closer together than the tolerance, all
are labeled as duplicates even if not all values in the run are within
the tolerance.

## Usage

``` r
duplicated_tol(x, tol = 1e-06, incl_first = FALSE)
```

## Arguments

- x:

  a numeric vector in which to identify duplicates

- tol:

  numeric tolerance for identifying duplicates

- incl_first:

  boolean indicator of whether or not the first entry in a run of
  duplicates should be indicated as a duplicate. `FALSE` mirrors the
  behavior of the base R function `duplicated`.

## Value

a boolean vector of the same length as `x`

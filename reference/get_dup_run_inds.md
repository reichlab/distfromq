# Get indices of starts and ends of runs of duplicate values

Get indices of starts and ends of runs of duplicate values

## Usage

``` r
get_dup_run_inds(dups)
```

## Arguments

- dups:

  a boolean vector that would result from calling
  `duplicated_tol(..., incl_first = FALSE)`

## Value

named list with entries `starts` giving indices of the first element in
each sequence of runs of duplicate values and `ends` giving indices of
the last element in each sequence of runs of duplicate values.

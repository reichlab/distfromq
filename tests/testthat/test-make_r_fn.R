test_that("make_r_fn errors with out-of-bounds or incorrectly typed ps, qs", {
  testthat::expect_no_error(make_r_fn(ps = c(0.0, 0.5, 1.0), qs = 1:3))
  testthat::expect_error(make_r_fn(ps = c(-1, 0.5, 1.0), qs = 1:3),
                         "Assertion on 'ps' failed: Element 1 is not >= 0.")
  testthat::expect_error(make_r_fn(ps = c(0.0, 0.5, 2.0), qs = 1:3),
                         "Assertion on 'ps' failed: Element 3 is not <= 1.")
  testthat::expect_error(make_r_fn(ps = c(0.0, "a", 1.0), qs = 1:3),
                         "Assertion on 'ps' failed: Must be of type 'numeric', not 'character'.")
  testthat::expect_error(make_r_fn(ps = c(0.0, 0.5, 1.0), qs = c(1, "a", 3)),
                         "Assertion on 'qs' failed: Must be of type 'numeric', not 'character'.")
  testthat::expect_error(make_r_fn(ps = c(0.0, 0.5, 1.0), qs = 1:4),
                         "'ps' and 'qs' must have the same length.")
})

test_that("make_r_fn result errors with out-of-bounds or incorrectly typed argument n", {
  r_fn <- make_r_fn(ps = c(0.0, 0.5, 1.0), qs = 1:3)
  testthat::expect_no_error(r_fn(5))
  testthat::expect_no_error(r_fn(1))
  testthat::expect_no_error(r_fn(0))
  testthat::expect_error(r_fn("a"),
                         "Must be of type 'integerish', not 'character'.")
  testthat::expect_error(r_fn(-1),
                         "Assertion on 'n' failed: Element 1 is not >= 0.")
  testthat::expect_error(r_fn(c(1L, 5L)),
                         "Assertion on 'n' failed: Must have length 1, but has length 2.")
})

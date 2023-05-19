test_that("make_d_fn throws error, all qs duplicated", {
    expect_error(make_d_fn(ps = c(.25, .5, .75), qs = c(0, 0, 0)))
})

test_that("make_d_fn throws error, only one q provided", {
    expect_error(make_d_fn(ps = .5, qs = 0))
})

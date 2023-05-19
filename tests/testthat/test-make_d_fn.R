test_that("make_q_fn works, all qs duplicated", {
    Q <- make_q_fn(ps = c(.25, .5, .75), qs = c(0, 0, 0))
    
    expect_equal(Q(0.5), 0.0)
    expect_equal(Q(c(0.01, 0.25, 0.29, 0.5)), rep(0.0, 4))
})

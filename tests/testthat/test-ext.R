test_that("q_ext works, repeated upper quantiles", {
  ps <- c(0.975, 0.99)
  qs <- c(28, 28)
  qext_fun <- q_ext_factory(ps = ps, qs = qs, dist = "norm")
  expect_true(
    isTRUE(all.equal(qext_fun(seq(from = 0.99, to = 1.0, length.out = 10)),
                     rep(28, 10)))
  )
})

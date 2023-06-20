test_that("q_ext works, repeated upper quantiles", {
    ps <- c(0.975, 0.99)
    qs <- c(28, 28)
    qext_fun <- q_ext_factory(ps = ps, qs = qs, dist = "norm")
    expect_true(
        isTRUE(all.equal(qext_fun(seq(from = 0.99, to = 1.0, length.out = 10)),
                  rep(28, 10)))
    )
})

test_that("q_ext works, lognormal family, non-zero quantiles", {
    ps <- c(0.975, 0.99)
    qs <- qlnorm(ps, meanlog = 0.5, sdlog = 2)
    qext_fun <- q_ext_factory(ps = ps, qs = qs, dist = "lnorm")

    test_ps <- seq(from = 0.99, to = 0.9999, length.out = 10)
    expect_true(
        isTRUE(all.equal(qext_fun(test_ps),
                         qlnorm(test_ps, meanlog = 0.5, sdlog = 2)))
    )
})

test_that("q_ext works, lognormal family, zero quantile(s) lower tail", {
    ps <- c(0.975, 0.99)
    qs <- c(0, 28)
    qext_fun <- q_ext_factory(ps = ps, qs = qs, dist = "lnorm", is_lower = TRUE)

    test_ps <- seq(from = 0.01, to = 0.97, length.out = 10)
    expect_true(
        isTRUE(all.equal(qext_fun(test_ps),
                         rep(0, 10)))
    )
})

test_that("q_ext generates error, lognormal family, negative quantile", {
    ps <- c(0.975, 0.99)
    qs <- c(-1, 28)

    expect_error(
        qext_fun <- q_ext_factory(ps = ps, qs = qs, dist = "lnorm")
    )
})

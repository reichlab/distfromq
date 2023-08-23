# Tests in this file are for 5 scenarios:
#  1. qs distinct, norm and lnorm
#  2. qs approx equal to non-zero value, norm and lnorm
#  3. qs approx equal to zero, norm and lnorm
#  4. q[1] approx 0, q[2] > 0, norm and lnorm
#  5. q[1] < 0, q[2] approx 0, norm and lnorm
#
# We evaluate each of d_ext_factory, p_ext_factory, and q_ext_factory
# under each scenario.
# d_ext_factory should generate an error under scenarios 2-5; the other methods
# should generate errors under scenario 5 when using an lnorm distribution, but
# work otherwise

qs3 <- c(0, 0 + 1e-13)
qs4 <- c(0 + 1e-13, 2)
qs5 <- c(-2, 0 + 1e-13)

test_that("ext factories work, qs distinct, norm", {
    ps <- c(0.1, 0.2)
    qs <- qnorm(ps, mean = 2, sd = 2)

    ps_test <- seq(from = 0.01, to = 0.99, by = 0.01)
    qs_test <- seq(from = -5, to = 5, length.out = 101)

    dext_fun <- d_ext_factory(ps = ps, qs = qs, dist = "norm",
                              dup_tol = 1e-6, zero_tol = 1e-12,
                              zero_discrete = TRUE,
                              is_lower = FALSE)
    pext_fun <- p_ext_factory(ps = ps, qs = qs, dist = "norm",
                              dup_tol = 1e-6, zero_tol = 1e-12,
                              zero_discrete = TRUE,
                              is_lower = FALSE)
    qext_fun <- q_ext_factory(ps = ps, qs = qs, dist = "norm",
                              dup_tol = 1e-6, zero_tol = 1e-12,
                              zero_discrete = FALSE,
                              is_lower = FALSE)

    expect_true(
        isTRUE(all.equal(dext_fun(qs_test), dnorm(qs_test, mean = 2, sd = 2)))
    )
    expect_true(
        isTRUE(all.equal(pext_fun(qs_test), pnorm(qs_test, mean = 2, sd = 2)))
    )
    expect_true(
        isTRUE(all.equal(qext_fun(ps_test), qnorm(ps_test, mean = 2, sd = 2)))
    )
})


test_that("ext factories work, qs distinct, lnorm", {
    ps <- c(0.1, 0.2)
    qs <- qlnorm(ps, meanlog = 2, sdlog = 2)

    ps_test <- seq(from = 0.01, to = 0.99, by = 0.01)
    qs_test <- seq(from = 0.01, to = 5, length.out = 101)

    dext_fun <- d_ext_factory(ps = ps, qs = qs, dist = "lnorm",
                              dup_tol = 1e-6, zero_tol = 1e-12,
                              zero_discrete = TRUE,
                              is_lower = FALSE)
    pext_fun <- p_ext_factory(ps = ps, qs = qs, dist = "lnorm",
                              dup_tol = 1e-6, zero_tol = 1e-12,
                              zero_discrete = TRUE,
                              is_lower = FALSE)
    qext_fun <- q_ext_factory(ps = ps, qs = qs, dist = "lnorm",
                              dup_tol = 1e-6, zero_tol = 1e-12,
                              zero_discrete = FALSE,
                              is_lower = FALSE)

    expect_true(
        isTRUE(all.equal(dext_fun(qs_test), dlnorm(qs_test, mean = 2, sd = 2)))
    )
    expect_true(
        isTRUE(all.equal(pext_fun(qs_test), plnorm(qs_test, mean = 2, sd = 2)))
    )
    expect_true(
        isTRUE(all.equal(qext_fun(ps_test), qlnorm(ps_test, mean = 2, sd = 2)))
    )
})


test_that("ext factories work, qs non-distinct non-zero, norm and lnorm", {
    ps <- c(0.1, 0.2)
    qs <- c(2, 2 + 1e-8)

    ps_test <- seq(from = 0.01, to = 0.99, by = 0.01)
    qs_test <- seq(from = -5, to = 5, length.out = 101)

    for (dist in c("norm", "lnorm")) {
        expect_error(
            dext_fun <- d_ext_factory(ps = ps, qs = qs, dist = dist,
                                    dup_tol = 1e-6, zero_tol = 1e-12,
                                    zero_discrete = TRUE,
                                    is_lower = FALSE)
        )
        pext_fun <- p_ext_factory(ps = ps, qs = qs, dist = dist,
                                dup_tol = 1e-6, zero_tol = 1e-12,
                                zero_discrete = TRUE,
                                is_lower = FALSE)
        qext_fun <- q_ext_factory(ps = ps, qs = qs, dist = dist,
                                dup_tol = 1e-6, zero_tol = 1e-12,
                                zero_discrete = FALSE,
                                is_lower = FALSE)

        expect_true(
            isTRUE(all.equal(pext_fun(qs_test), as.numeric(qs_test >= 2.0)))
        )
        expect_true(
            isTRUE(all.equal(qext_fun(ps_test), rep(2.0, length(ps_test))))
        )
    }
})


test_that("ext factories work, qs non-distinct at zero, norm and lnorm", {
    ps <- c(0.1, 0.2)
    qs <- c(0, 1e-13)

    ps_test <- seq(from = 0.01, to = 0.99, by = 0.01)
    qs_test <- seq(from = -5, to = 5, length.out = 101)

    for (dist in c("norm", "lnorm")) {
        for (is_lower in c(TRUE, FALSE)) {
            expect_error(
                dext_fun <- d_ext_factory(ps = ps, qs = qs, dist = dist,
                                        dup_tol = 1e-6, zero_tol = 1e-12,
                                        zero_discrete = TRUE,
                                        is_lower = is_lower)
            )
            pext_fun <- p_ext_factory(ps = ps, qs = qs, dist = dist,
                                    dup_tol = 1e-6, zero_tol = 1e-12,
                                    zero_discrete = (dist == "lnorm"),
                                    is_lower = is_lower)
            qext_fun <- q_ext_factory(ps = ps, qs = qs, dist = dist,
                                    dup_tol = 1e-6, zero_tol = 1e-12,
                                    zero_discrete = (dist == "lnorm"),
                                    is_lower = is_lower)

            expect_true(
                isTRUE(all.equal(pext_fun(qs_test), as.numeric(qs_test >= 0.0)))
            )
            expect_true(
                isTRUE(all.equal(qext_fun(ps_test), rep(0.0, length(ps_test))))
            )
        }
    }
})


test_that("ext factories work, q[1] approx zero, q[2] > 0, norm and lnorm", {
    qs <- c(1e-13, 2)
    ps <- pnorm(qs, mean = 2, sd = 2)

    ps_test <- seq(from = 0.01, to = 0.99, by = 0.01)
    qs_test <- seq(from = -5, to = 5, length.out = 101)

    for (dist in c("norm", "lnorm")) {
        for (is_lower in c(TRUE, FALSE)) {
            if (dist == "lnorm") {
                expect_error(
                    dext_fun <- d_ext_factory(ps = ps, qs = qs, dist = dist,
                                            dup_tol = 1e-6, zero_tol = 1e-12,
                                            zero_discrete = (dist == "lnorm"),
                                            is_lower = is_lower)
                )
            } else {
                dext_fun <- d_ext_factory(ps = ps, qs = qs, dist = dist,
                                          dup_tol = 1e-6, zero_tol = 1e-12,
                                          zero_discrete = (dist == "lnorm"),
                                          is_lower = is_lower)
            }

            pext_fun <- p_ext_factory(ps = ps, qs = qs, dist = dist,
                                      dup_tol = 1e-6, zero_tol = 1e-12,
                                      zero_discrete = (dist == "lnorm"),
                                      is_lower = is_lower)
            qext_fun <- q_ext_factory(ps = ps, qs = qs, dist = dist,
                                      dup_tol = 1e-6, zero_tol = 1e-12,
                                      zero_discrete = (dist == "lnorm"),
                                      is_lower = is_lower)

            if (dist == "lnorm") {
                if (is_lower) {
                    expect_true(
                        isTRUE(all.equal(pext_fun(qs_test),
                                        as.numeric(qs_test >= 0.0)))
                    )
                    expect_true(
                        isTRUE(all.equal(qext_fun(ps_test),
                                        rep(0.0, length(ps_test))))
                    )
                } else {
                    expect_true(
                        isTRUE(all.equal(pext_fun(qs_test),
                                        as.numeric(qs_test >= 2.0)))
                    )
                    expect_true(
                        isTRUE(all.equal(qext_fun(ps_test),
                                        rep(2.0, length(ps_test))))
                    )
                }
            } else {
                # effectively duplicating checks done in first scenario
                expect_true(
                    isTRUE(all.equal(dext_fun(qs_test),
                                     dnorm(qs_test, mean = 2, sd = 2)))
                )
                expect_true(
                    isTRUE(all.equal(pext_fun(qs_test),
                                     pnorm(qs_test, mean = 2, sd = 2)))
                )
                expect_true(
                    isTRUE(all.equal(qext_fun(ps_test),
                                     qnorm(ps_test, mean = 2, sd = 2)))
                )
            }
        }
    }
})


test_that("ext factories work, q[1] < 0, q[2] approx 0, norm and lnorm", {
    qs <- c(-1, 1e-13)
    ps <- pnorm(qs, mean = 2, sd = 2)

    ps_test <- seq(from = 0.01, to = 0.99, by = 0.01)
    qs_test <- seq(from = -5, to = 5, length.out = 101)

    for (dist in c("norm", "lnorm")) {
        for (is_lower in c(TRUE, FALSE)) {
            if (dist == "lnorm") {
                expect_error(
                    dext_fun <- d_ext_factory(ps = ps, qs = qs, dist = dist,
                                            dup_tol = 1e-6, zero_tol = 1e-12,
                                            zero_discrete = (dist == "lnorm"),
                                            is_lower = is_lower)
                )
                expect_error(
                    pext_fun <- p_ext_factory(ps = ps, qs = qs, dist = dist,
                                            dup_tol = 1e-6, zero_tol = 1e-12,
                                            zero_discrete = (dist == "lnorm"),
                                            is_lower = is_lower)
                )
                expect_error(
                    qext_fun <- q_ext_factory(ps = ps, qs = qs, dist = dist,
                                            dup_tol = 1e-6, zero_tol = 1e-12,
                                            zero_discrete = (dist == "lnorm"),
                                            is_lower = is_lower)
                )
            } else {
                dext_fun <- d_ext_factory(ps = ps, qs = qs, dist = dist,
                                          dup_tol = 1e-6, zero_tol = 1e-12,
                                          zero_discrete = (dist == "lnorm"),
                                          is_lower = is_lower)
                pext_fun <- p_ext_factory(ps = ps, qs = qs, dist = dist,
                                        dup_tol = 1e-6, zero_tol = 1e-12,
                                        zero_discrete = (dist == "lnorm"),
                                        is_lower = is_lower)
                qext_fun <- q_ext_factory(ps = ps, qs = qs, dist = dist,
                                        dup_tol = 1e-6, zero_tol = 1e-12,
                                        zero_discrete = (dist == "lnorm"),
                                        is_lower = is_lower)

                # effectively duplicating checks done in first scenario
                expect_true(
                    isTRUE(all.equal(dext_fun(qs_test),
                                     dnorm(qs_test, mean = 2, sd = 2)))
                )
                expect_true(
                    isTRUE(all.equal(pext_fun(qs_test),
                                     pnorm(qs_test, mean = 2, sd = 2)))
                )
                expect_true(
                    isTRUE(all.equal(qext_fun(ps_test),
                                     qnorm(ps_test, mean = 2, sd = 2)))
                )
            }
        }
    }
})

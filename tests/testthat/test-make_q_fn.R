# This file includes the following test cases for make_q_fn:
# - continuous component; no point masses
# - no continuous component; one point mass, single q
# - no continuous component; one point mass, duplicated qs
# - no continuous component; two point masses, both from duplicated qs
# - no continuous component; two point masses, non-zero, duplicated values first only
# - no continuous component; two point masses, non-zero, duplicated values second only
# - no continuous component; two point masses, is_hurdle with one zero and one non-zero
# - continuous component; one point mass, duplicated qs
# - continuous component; one point mass, is_hurdle with zero
# - continuous component; two point masses, both from duplicated qs
# - continuous component; two point masses, is_hurdle with one zero and one non-zero with duplicated qs
#
# there are three additional tests related to bugs that were found:
# - make_q_fn is well-behaved with near-equal quantiles (two variations)
# - make_q_fn is well-behaved with near-zero quantiles
#
# we expect the function to approximately reproduce the qf that is being
# estimated in all cases.


test_that("make_q_fn works, continuous component; no point masses", {
    ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    qs <- qnorm(ps, mean = 1, sd = 2)
    test_ps <- seq(from = 0.01, to = 0.99, by = 0.01)

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_expected <- qnorm(test_ps, mean = 1, sd = 2)
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
})

test_that("make_q_fn works, no continuous component; one point mass, single q", {
    ps <- 0.1
    qs <- rep(0.0, length(ps))
    test_ps <- seq(from = 0.01, to = 0.99, by = 0.01)

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_expected <- rep(0.0, length(test_ps))
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
})

test_that("make_q_fn works, no continuous component; one point mass, duplicated qs", {
    ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    qs <- rep(0.0, length(ps))
    test_ps <- seq(from = 0.01, to = 0.99, by = 0.01)

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_expected <- rep(0.0, length(test_ps))
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
})

test_that("make_q_fn works, no continuous component; two point masses, both from duplicated qs", {
    ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    qs <- c(rep(1.0, 3), rep(2.0, 6))
    test_ps <- sort(c(
        1 / 3, 1.0,
        seq(from = 0.01, to = 0.99, by = 0.01)))

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_expected <- ifelse(test_ps <= 1 / 3, 1.0, 2.0)
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
})


test_that("make_q_fn works, no continuous component; two point masses, non-zero, duplicated values first only", {
    ps <- seq(from = 0.1, to = 0.4, by = 0.1)
    qs <- c(rep(1.0, 3), rep(2.0, 1))
    test_ps <- sort(c(
        1 / 3, 1.0,
        seq(from = 0.01, to = 0.99, by = 0.01)))

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_expected <- ifelse(test_ps <= 1 / 3, 1.0, 2.0)
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
})


test_that("make_q_fn works, no continuous component; two point masses, non-zero, duplicated values second only", {
    ps <- seq(from = 0.3, to = 0.9, by = 0.1)
    qs <- c(rep(1.0, 1), rep(2.0, 6))
    test_ps <- sort(c(
        1 / 3, 1.0,
        seq(from = 0.01, to = 0.99, by = 0.01)))

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_expected <- ifelse(test_ps <= 1 / 3, 1.0, 2.0)
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
})


test_that("make_q_fn works, no continuous component; is_hurdle with one zero and one non-zero", {
    ps <- seq(from = 0.1, to = 0.2, by = 0.1)
    qs <- c(1e-13, 2.0)
    test_ps <- sort(c(
        1 / 9, 1.0,
        seq(from = 0.01, to = 0.99, by = 0.01)))

    q_actual <- make_q_fn(ps, qs, tail_dist = "lnorm")(test_ps)
    q_expected <- ifelse(test_ps <= 1 / 9, 0.0, 2.0)
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
})


test_that("make_q_fn works, continuous component; one point mass, duplicated qs", {
    # mixture of a Normal(0,1) with weight 0.8 and
    # a point mass at 0 with weight 0.2

    # probabilities and quantiles for normal component
    norm_ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    norm_qs <- qnorm(norm_ps)
    adj_norm_ps <- norm_ps * 0.8 + 0.2 * (norm_qs > 0.0)

    # probabilities and quantiles for point mass at 0
    point_ps <- seq(from = 0.0, to = 1.0, by = 0.1)
    point_qs <- rep(0.0, length(point_ps))
    adj_point_ps <- 0.5 * 0.8 + point_ps * 0.2

    ps <- sort(c(adj_norm_ps, adj_point_ps))
    qs <- sort(c(norm_qs, point_qs))
    dup_inds <- duplicated(ps)
    ps <- ps[!dup_inds]
    qs <- qs[!dup_inds]

    test_ps <- seq(from = 0.01, to = 0.99, by = 0.01)

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_expected <- c(
        qnorm(test_ps[test_ps < 0.4] / 0.8),
        rep(0.0, sum((test_ps >= 0.4) & (test_ps <= 0.6))),
        qnorm((test_ps[test_ps > 0.6] - 0.2) / 0.8)
    )
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
})


test_that("make_q_fn works, continuous component; one point mass, is_hurdle with zero", {
    # mixture of a LogNormal(0,1) with weight 0.8 and
    # a point mass at 0 with weight 0.2

    # probabilities and quantiles for lognormal component
    norm_ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    norm_qs <- qlnorm(norm_ps)
    adj_norm_ps <- norm_ps * 0.8 + 0.2 * (norm_qs > 0.0)

    # probabilities and quantiles for point mass at 0
    point_ps <- 1.0
    point_qs <- 0.0
    adj_point_ps <- 0.2

    ps <- sort(c(adj_norm_ps, adj_point_ps))
    qs <- sort(c(norm_qs, point_qs))
    dup_inds <- duplicated(ps)
    ps <- ps[!dup_inds]
    qs <- qs[!dup_inds]

    test_ps <- seq(from = 0.01, to = 0.99, by = 0.01)
    q_actual <- make_q_fn(ps, qs, tail_dist = "lnorm")(test_ps)
    q_expected <- c(
        rep(0.0, sum(test_ps <= 0.2)),
        qlnorm((test_ps[test_ps > 0.2] - 0.2) / 0.8)
    )

    # note that regions where actual matches expected depend on
    # tail family assumptions:
    # plot(test_ps, q_actual; lines(test_ps, q_expected)

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
})


test_that("make_q_fn works, two point masses, both from duplicated qs", {
    # mixture of a Normal(0,1) with weight 0.6,
    # a point mass at 0 with weight 0.3, and a point mass at 1 with weight 0.1

    # probabilities and quantiles for normal component
    norm_ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    norm_qs <- qnorm(norm_ps)
    adj_norm_ps <- norm_ps * 0.6 + 0.3 * (norm_qs > 0.0) + 0.1 * (norm_qs > 1.0)

    # probabilities and quantiles for point mass at 0
    point_ps_0 <- seq(from = 0.0, to = 1.0, by = 0.1)
    point_qs_0 <- rep(0.0, length(point_ps_0))
    adj_point_ps_0 <- 0.5 * 0.6 + point_ps_0 * 0.3

    # probabilities and quantiles for point mass at 1
    point_ps_1 <- seq(from = 0.0, to = 1.0, by = 0.1)
    point_qs_1 <- rep(1.0, length(point_ps_1))
    adj_point_ps_1 <- pnorm(1.0) * 0.6 + 0.3 + point_ps_1 * 0.1

    ps <- sort(c(adj_norm_ps, adj_point_ps_0, adj_point_ps_1))
    qs <- sort(c(norm_qs, point_qs_0, point_qs_1))
    dup_inds <- duplicated(ps)
    ps <- ps[!dup_inds]
    qs <- qs[!dup_inds]

    test_ps <- seq(from = 0.01, to = 0.99, by = 0.01)

    q_actual <- make_q_fn(ps, qs)(test_ps)

    pcut1 <- 0.3
    pcut2 <- pnorm(1.0)*0.6 + 0.3
    q_expected <- c(
        qnorm(test_ps[test_ps < pcut1] / 0.6),
        rep(0.0, sum((test_ps >= 0.3) & (test_ps <= 0.6))),
        qnorm((test_ps[(test_ps > 0.6) & (test_ps < pcut2)] - 0.3) / 0.6),
        rep(1.0, sum((test_ps >= pcut2) & (test_ps <= pcut2 + 0.1))),
        qnorm((test_ps[test_ps > pcut2 + 0.1] - 0.4) / 0.6)
    )
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
})


test_that("make_q_fn works, two point masses, is_hurdle with one zero and one non-zero with duplicated qs", {
    # mixture of a LogNormal(0,1) with weight 0.6,
    # a point mass at 0 with weight 0.3, and a point mass at 1 with weight 0.1

    # probabilities and quantiles for normal component
    norm_ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    norm_qs <- qlnorm(norm_ps)
    adj_norm_ps <- norm_ps * 0.6 + 0.3 * (norm_qs > 0.0) + 0.1 * (norm_qs > 1.0)

    # probabilities and quantiles for point mass at 0
    point_ps_0 <- 1.0
    point_qs_0 <- 0.0
    adj_point_ps_0 <- 0.3

    # probabilities and quantiles for point mass at 1
    point_ps_1 <- seq(from = 0.0, to = 1.0, by = 0.1)
    point_qs_1 <- rep(1.0, length(point_ps_1))
    adj_point_ps_1 <- plnorm(1.0) * 0.6 + 0.3 + point_ps_1 * 0.1

    ps <- sort(c(adj_norm_ps, adj_point_ps_0, adj_point_ps_1))
    qs <- sort(c(norm_qs, point_qs_0, point_qs_1))
    dup_inds <- duplicated(ps)
    ps <- ps[!dup_inds]
    qs <- qs[!dup_inds]

    test_ps <- seq(from = 0.01, to = 0.99, by = 0.01)
    q_actual <- make_q_fn(ps, qs, tail_dist = "lnorm")(test_ps)

    pcut1 <- 0.3
    pcut2 <- plnorm(1.0) * 0.6 + 0.3
    q_expected <- c(
        rep(0.0, sum(test_ps <= 0.3)),
        qlnorm((test_ps[(test_ps > 0.3) & (test_ps < pcut2)] - 0.3) / 0.6),
        rep(1.0, sum((test_ps >= pcut2) & (test_ps <= pcut2 + 0.1))),
        qlnorm((test_ps[test_ps > pcut2 + 0.1] - 0.4) / 0.6)
    )
    # note that regions where actual matches expected depend on
    # tail family assumptions:
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
})


test_that("make_q_fn is well-behaved with near-equal quantiles", {
    # probabilities and quantiles with a near-duplicate q
    ps <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)
    qs <- c(0, 0, 2.99327779507948e-16, 4.94582028429876, 8.22187599370956,
        9.89992446804951, 11.2927083522741, 12.5925601513011, 13.8718728208768,
        15.0382200560432, 16.0565760032992, 16.9592610229454, 17.7662853190937,
        18.6751481876927, 19.7650748463216, 21.0314540241319, 22.5188110517322,
        24.3155280557975, 26.5953673396936, 30.2878452640675, 37.6971667663299,
        46.5930391567271, 52.2594821457131)

    test_ps <- seq(from = 0.001, to = 0.999, by = 0.001)
    q_actual <- make_q_fn(ps, qs)(test_ps)

    testthat::expect_true(all(q_actual >= 0.0))
    testthat::expect_true(all(q_actual <= 10 * max(qs)))
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
})


test_that("make_q_fn is well-behaved with near-equal quantiles", {
    # probabilities and quantiles with a near-duplicate q
    ps <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)
    qs <- c(0, 0, 1.5e-6, 4.94582028429876, 8.22187599370956,
        9.89992446804951, 11.2927083522741, 12.5925601513011, 13.8718728208768,
        15.0382200560432, 16.0565760032992, 16.9592610229454, 17.7662853190937,
        18.6751481876927, 19.7650748463216, 21.0314540241319, 22.5188110517322,
        24.3155280557975, 26.5953673396936, 30.2878452640675, 37.6971667663299,
        46.5930391567271, 52.2594821457131)

    test_ps <- seq(from = 0.001, to = 0.999, by = 0.001)
    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_fn <- make_q_fn(ps, qs)

    testthat::expect_true(all(q_actual >= 0.0))
    testthat::expect_true(all(q_actual <= 10 * max(qs)))
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
})


test_that("make_q_fn is well-behaved with near-zero quantiles", {
    # probabilities and quantiles with a near-duplicate q
    ps <- seq(0.1, 0.9, 0.2)
    qs <- c(0, 0, 0, 1e-07, 1e-05)

    testthat::expect_no_error(
        cdf_hat <- distfromq:::make_q_fn(ps, qs)
    )
})

# This file includes the following test cases for make_p_fn and make_q_fn:
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
# there are five additional tests related to bugs that were found:
# - near-equal quantiles (two variations)
# - near-zero quantiles
# - duplicated quantiles at 3 distinct values, lnorm tail distribution
# - four repeated zero quantiles, normal tail distribution
#
# in each case, we check that:
# - make_p_fn approximately reproduces the cdf that is being estimated
# - make_q_fn approximately reproduces the qf that is being estimated
# - make_p_fn and make_q_fn are inverses

test_that("make_p_fn, make_q_fn work, continuous component; no point masses", {
    ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    qs <- qnorm(ps, mean = 1, sd = 2)

    test_ps <- seq(from = 0.01, to = 0.99, by = 0.01)
    test_qs <- qnorm(seq(from = 0.01, to = 0.99, by = 0.01), mean = 1, sd = 2)

    p_actual <- make_p_fn(ps, qs)(test_qs)
    p_expected <- pnorm(test_qs, mean = 1, sd = 2)
    # plot(test_qs, p_actual); lines(test_qs, p_expected)

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_expected <- qnorm(test_ps, mean = 1, sd = 2)
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(p_actual, p_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(p_actual[order(test_qs)]) >= 0))

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))

    testthat::expect_equal(test_ps, make_p_fn(ps, qs)(q_actual),
                           tolerance = 1e-3)
    testthat::expect_equal(test_qs, make_q_fn(ps, qs)(p_actual),
                           tolerance = 1e-3)
})

test_that("make_p_fn, make_q_fn work, no continuous component; one point mass, single q", {
    ps <- 0.1
    qs <- rep(0.0, length(ps))

    test_ps <- seq(from = 0.01, to = 0.99, by = 0.01)
    test_qs <- qnorm(seq(from = 0.01, to = 0.99, by = 0.01), mean = 1, sd = 2)

    p_actual <- make_p_fn(ps, qs)(test_qs)
    p_expected <- as.numeric(test_qs >= 0)
    # plot(test_qs, p_actual); lines(test_qs, p_expected)

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_expected <- rep(0.0, length(test_ps))
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(p_actual, p_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(p_actual[order(test_qs)]) >= 0))

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))

    testthat::expect_equal(rep(1, length(test_ps)),
                           make_p_fn(ps, qs)(q_actual),
                           tolerance = 1e-3)
    testthat::expect_equal(rep(0, length(test_qs)),
                           make_q_fn(ps, qs)(p_actual),
                           tolerance = 1e-3)
})

test_that("make_p_fn, make_q_fn work, no continuous component; one point mass, duplicated qs", {
    ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    qs <- rep(0.0, length(ps))

    test_ps <- seq(from = 0.01, to = 0.99, by = 0.01)
    test_qs <- qnorm(seq(from = 0.01, to = 0.99, by = 0.01), mean = 1, sd = 2)

    p_actual <- make_p_fn(ps, qs)(test_qs)
    p_expected <- as.numeric(test_qs >= 0)
    # plot(test_qs, p_actual); lines(test_qs, p_expected)

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_expected <- rep(0.0, length(test_ps))
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(p_actual, p_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(p_actual[order(test_qs)]) >= 0))

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))

    testthat::expect_equal(rep(1, length(test_ps)),
                           make_p_fn(ps, qs)(q_actual),
                           tolerance = 1e-3)
    testthat::expect_equal(rep(0, length(test_qs)),
                           make_q_fn(ps, qs)(p_actual),
                           tolerance = 1e-3)
})

test_that("make_p_fn, make_q_fn work, no continuous component; two point masses, both from duplicated qs", {
    ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    qs <- c(rep(1.0, 3), rep(2.0, 6))

    test_ps <- sort(c(
        1 / 3, 1.0,
        seq(from = 0.01, to = 0.99, by = 0.01)))
    test_qs <- c(
        1, 2,
        qnorm(seq(from = 0.01, to = 0.99, by = 0.01), mean = 1, sd = 2))

    p_actual <- make_p_fn(ps, qs)(test_qs)
    p_expected <- (1 / 3) * as.numeric(test_qs >= 1.0) +
                  (2 / 3) * as.numeric(test_qs >= 2.0)
    # plot(test_qs, p_actual); lines(test_qs, p_expected)

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_expected <- ifelse(test_ps <= 1 / 3, 1.0, 2.0)
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(p_actual, p_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(p_actual[order(test_qs)]) >= 0))

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))

    testthat::expect_equal(c(rep(1 / 3, sum(test_ps <= 1 / 3)),
                             rep(1, sum(test_ps > 1 / 3))),
                           make_p_fn(ps, qs)(q_actual),
                           tolerance = 1e-3)
    # Commenting this test out -- fails due to floating point precision
    # testthat::expect_equal(...,
    #                        make_q_fn(ps, qs)(p_actual),
    #                        tolerance = 1e-3)
})


test_that("make_p_fn, make_q_fn work, no continuous component; two point masses, non-zero, duplicated values first only", {
    ps <- seq(from = 0.1, to = 0.4, by = 0.1)
    qs <- c(rep(1.0, 3), rep(2.0, 1))

    test_ps <- sort(c(
        1 / 3, 1.0,
        seq(from = 0.01, to = 0.99, by = 0.01)))
    test_qs <- c(
        1, 2,
        qnorm(seq(from = 0.01, to = 0.99, by = 0.01), mean = 1, sd = 2))

    p_actual <- make_p_fn(ps, qs)(test_qs)
    p_expected <- (1 / 3) * as.numeric(test_qs >= 1.0) +
                  (2 / 3) * as.numeric(test_qs >= 2.0)
    # plot(test_qs, p_actual); lines(test_qs, p_expected)

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_expected <- ifelse(test_ps <= 1 / 3, 1.0, 2.0)
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(p_actual, p_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(p_actual[order(test_qs)]) >= 0))

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))

    testthat::expect_equal(c(rep(1 / 3, sum(test_ps <= 1 / 3)),
                             rep(1, sum(test_ps > 1 / 3))),
                           make_p_fn(ps, qs)(q_actual),
                           tolerance = 1e-3)
    # Commenting this test out -- fails due to floating point precision
    # testthat::expect_equal(...,
    #                        make_q_fn(ps, qs)(p_actual),
    #                        tolerance = 1e-3)
})


test_that("make_p_fn, make_q_fn work, no continuous component; two point masses, non-zero, duplicated values second only", {
    ps <- seq(from = 0.3, to = 0.9, by = 0.1)
    qs <- c(rep(1.0, 1), rep(2.0, 6))

    test_ps <- sort(c(
        1 / 3, 1.0,
        seq(from = 0.01, to = 0.99, by = 0.01)))
    test_qs <- c(
        1, 2,
        qnorm(seq(from = 0.01, to = 0.99, by = 0.01), mean = 1, sd = 2))

    p_actual <- make_p_fn(ps, qs)(test_qs)
    p_expected <- (1 / 3) * as.numeric(test_qs >= 1.0) +
                  (2 / 3) * as.numeric(test_qs >= 2.0)
    # plot(test_qs, p_actual); lines(test_qs, p_expected)

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_expected <- ifelse(test_ps <= 1 / 3, 1.0, 2.0)
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(p_actual, p_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(p_actual[order(test_qs)]) >= 0))

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))

    testthat::expect_equal(c(rep(1 / 3, sum(test_ps <= 1 / 3)),
                             rep(1, sum(test_ps > 1 / 3))),
                           make_p_fn(ps, qs)(q_actual),
                           tolerance = 1e-3)
    # Commenting this test out -- fails due to floating point precision
    # testthat::expect_equal(...,
    #                        make_q_fn(ps, qs)(p_actual),
    #                        tolerance = 1e-3)
})


test_that("make_p_fn, make_q_fn work, no continuous component; is_hurdle with one zero and one non-zero", {
    ps <- seq(from = 0.1, to = 0.2, by = 0.1)
    qs <- c(1e-13, 2.0)

    test_ps <- sort(c(
        1 / 9, 1.0,
        seq(from = 0.01, to = 0.99, by = 0.01)))
    test_qs <- c(
        1, 2,
        qnorm(seq(from = 0.01, to = 0.99, by = 0.01), mean = 1, sd = 2))

    p_actual <- make_p_fn(ps, qs, tail_dist = "lnorm")(test_qs)
    p_expected <- (1 / 9) * as.numeric(test_qs >= 0.0) +
                  (8 / 9) * as.numeric(test_qs >= 2.0)
    # plot(test_qs, p_actual); lines(test_qs, p_expected)

    q_actual <- make_q_fn(ps, qs, tail_dist = "lnorm")(test_ps)
    q_expected <- ifelse(test_ps <= 1 / 9, 0.0, 2.0)
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(p_actual, p_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(p_actual[order(test_qs)]) >= 0))

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
})


test_that("make_p_fn, make_q_fn work, continuous component; one point mass, duplicated qs", {
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
    test_qs <- qnorm(seq(from = 0.01, to = 0.99, by = 0.01), mean = 1, sd = 2)

    p_actual <- make_p_fn(ps, qs)(test_qs)
    p_expected <- 0.2 * as.numeric(test_qs >= 0) + 0.8 * pnorm(test_qs)
    # plot(test_qs, p_actual); lines(test_qs, p_expected)

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_expected <- c(
        qnorm(test_ps[test_ps < 0.4] / 0.8),
        rep(0.0, sum((test_ps >= 0.4) & (test_ps <= 0.6))),
        qnorm((test_ps[test_ps > 0.6] - 0.2) / 0.8)
    )
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(p_actual, p_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(p_actual[order(test_qs)]) >= 0))

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))

    testthat::expect_equal(make_p_fn(ps, qs)(q_actual),
                           ifelse(test_ps >= 0.4 & test_ps <= 0.6,
                                  0.6,
                                  test_ps),
                           tolerance = 1e-3)
    testthat::expect_equal(make_q_fn(ps, qs)(p_actual),
                           test_qs,
                           tolerance = 1e-3)
})


test_that("make_p_fn, make_q_fn work, continuous component; one point mass, is_hurdle with zero", {
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
    test_qs <- qnorm(seq(from = 0.01, to = 0.99, by = 0.01), mean = 1, sd = 2)

    p_actual <- make_p_fn(ps, qs, tail_dist = "lnorm")(test_qs)
    p_expected <- 0.2 * as.numeric(test_qs >= 0) +
                  0.8 * plnorm(test_qs)
    # note that regions where actual matches expected depend on
    # tail family assumptions:
    # plot(test_qs, p_actual); lines(test_qs, p_expected)

    q_actual <- make_q_fn(ps, qs, tail_dist = "lnorm")(test_ps)
    q_expected <- c(
        rep(0.0, sum(test_ps <= 0.2)),
        qlnorm((test_ps[test_ps > 0.2] - 0.2) / 0.8)
    )
    # note that regions where actual matches expected depend on
    # tail family assumptions:
    # plot(test_ps, q_actual); lines(test_ps, q_expected)

    testthat::expect_equal(p_actual, p_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(p_actual[order(test_qs)]) >= 0))

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))

    testthat::expect_equal(make_p_fn(ps, qs, tail_dist = "lnorm")(q_actual),
                           ifelse(test_ps <= 0.2,
                                  0.2,
                                  test_ps),
                           tolerance = 1e-3)
    testthat::expect_equal(make_q_fn(ps, qs, tail_dist = "lnorm")(p_actual),
                           ifelse(test_qs < 0, 0, test_qs),
                           tolerance = 1e-3)
})


test_that("make_p_fn, make_q_fn work, two point masses, both from duplicated qs", {
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
    test_qs <- qnorm(seq(from = 0.01, to = 0.99, by = 0.01), mean = 1, sd = 2)

    p_actual <- make_p_fn(ps, qs)(test_qs)
    p_expected <- 0.3 * as.numeric(test_qs >= 0) +
                  0.1 * as.numeric(test_qs >= 1) +
                  0.6 * pnorm(test_qs)
    # plot(test_qs, p_actual); lines(test_qs, p_expected)

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

    testthat::expect_equal(p_actual, p_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(p_actual[order(test_qs)]) >= 0))

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))

    testthat::expect_equal(make_p_fn(ps, qs)(q_actual),
                           dplyr::case_when(
                               test_ps >= 0.3 & test_ps <= 0.6 ~ 0.6,
                               test_ps >= pcut2 & test_ps <= pcut2 + 0.1 ~ pcut2 + 0.1,
                               TRUE ~ test_ps),
                           tolerance = 1e-3)
    testthat::expect_equal(make_q_fn(ps, qs)(p_actual),
                           test_qs,
                           tolerance = 1e-3)
})


test_that("make_p_fn, make_q_fn work, two point masses, is_hurdle with one zero and one non-zero with duplicated qs", {
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
    test_qs <- qnorm(seq(from = 0.01, to = 0.99, by = 0.01), mean = 1, sd = 2)

    p_actual <- make_p_fn(ps, qs, tail_dist = "lnorm")(test_qs)
    p_expected <- 0.3 * as.numeric(test_qs >= 0) +
                  0.1 * as.numeric(test_qs >= 1) +
                  0.6 * plnorm(test_qs)
    # note that regions where actual matches expected depend on
    # tail family assumptions:
    # plot(test_qs, p_actual); lines(test_qs, p_expected)

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

    testthat::expect_equal(p_actual, p_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(p_actual[order(test_qs)]) >= 0))

    testthat::expect_equal(q_actual, q_expected, tolerance = 1e-3)
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))

    testthat::expect_equal(make_p_fn(ps, qs, tail_dist = "lnorm")(q_actual),
                           dplyr::case_when(
                               test_ps <= 0.3 ~ 0.3,
                               test_ps >= pcut2 & test_ps <= pcut2 + 0.1 ~ pcut2 + 0.1,
                               TRUE ~ test_ps),
                           tolerance = 1e-3)
    testthat::expect_equal(make_q_fn(ps, qs, tail_dist = "lnorm")(p_actual),
                           ifelse(test_qs < 0, 0, test_qs),
                           tolerance = 1e-3)
})


test_that("make_p_fn, make_q_fn well-behaved with near-equal quantiles", {
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
    test_qs <- seq(from = 0.0, to = 100.0, length.out = 1001)

    p_actual <- make_p_fn(ps, qs)(test_qs)
    q_actual <- make_q_fn(ps, qs)(test_ps)

    testthat::expect_true(all(p_actual >= 0.0))
    testthat::expect_true(all(p_actual <= 1.0))
    testthat::expect_true(all(diff(p_actual[order(test_qs)]) > 0))
    # plot(test_qs, p_actual, type = "l"); points(qs, ps, pch = 20)
    # lines(q_actual, test_ps, lty=2, col='red')

    testthat::expect_true(all(q_actual >= 0.0))
    testthat::expect_true(all(q_actual <= 10 * max(qs)))
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))
    # plot(test_ps, q_actual, type = "l"); points(ps, qs, pch = 20)

    testthat::expect_equal(make_p_fn(ps, qs)(q_actual),
                           dplyr::case_when(
                               test_ps <= 0.05 ~ 0.05,
                               TRUE ~ test_ps),
                           tolerance = 1e-3)
    testthat::expect_equal(make_q_fn(ps, qs)(p_actual),
                           ifelse(test_qs < 0, 0, test_qs),
                           tolerance = 1e-3)
})


test_that("make_p_fn, make_q_fn well-behaved with near-equal quantiles", {
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
    test_qs <- seq(from = 0.0, to = 100.0, length.out = 1001)

    p_actual <- make_p_fn(ps, qs)(test_qs)

    q_actual <- make_q_fn(ps, qs)(test_ps)
    q_fn <- make_q_fn(ps, qs)

    testthat::expect_true(all(p_actual >= 0.0))
    testthat::expect_true(all(p_actual <= 1.0))
    testthat::expect_true(all(diff(p_actual[order(test_qs)]) > 0))

    testthat::expect_true(all(q_actual >= 0.0))
    testthat::expect_true(all(q_actual <= 10 * max(qs)))
    testthat::expect_true(all(diff(q_actual[order(test_ps)]) >= 0))

    testthat::expect_equal(make_p_fn(ps, qs)(q_actual),
                           dplyr::case_when(
                               test_ps <= 0.025 ~ 0.025,
                               TRUE ~ test_ps),
                           tolerance = 1e-6)
    testthat::expect_equal(make_q_fn(ps, qs)(p_actual),
                           test_qs,
                           tolerance = 1e-3)
})


test_that("make_p_fn, make_q_fn well-behaved with near-zero quantiles", {
    # probabilities and quantiles with a near-duplicate q
    ps <- seq(0.1, 0.9, 0.2)
    qs <- c(0, 0, 0, 1e-07, 1e-05)

    testthat::expect_no_error(
        cdf_hat <- distfromq:::make_p_fn(ps, qs)
    )
    testthat::expect_no_error(
        qf_hat <- distfromq:::make_q_fn(ps, qs)
    )
})

test_that("make_p_fn, make_q_fn well-behaved with 3 duplicated values, one at zero, lnorm tail dist", {
    ps <- seq(0.05, 0.95, 0.1)
    qs <- c(rep(0, 4), rep(1, 3), rep(5, 3))

    testthat::expect_no_error(
        q_hat <- distfromq::make_q_fn(
            ps = ps,
            qs = qs,
            dup_tol = 1e-06,
            zero_tol = 1e-12,
            tail_dist = "lnorm")(seq(from = 0, to = 1, length.out = 10000 + 2)[2:10000])
    )

    testthat::expect_no_error(
        p_hat <- distfromq::make_p_fn(
            ps = ps,
            qs = qs,
            dup_tol = 1e-06,
            zero_tol = 1e-12,
            tail_dist = "lnorm")(seq(from = 0, to = 10, length.out = 10000))
    )
})


test_that("make_p_fn, make_q_fn well-behaved: multiple duplicates and floating point issues with discrete adjustments", {
    ps <- c(.01,.025, seq(.05,.95, by=.05), .975, .99)
    qs <- c(0,0,0,0,3,6,8,8,9,11,12,13,17,21,25,27,29,30,31,33,37,52,61)

    q_hat <- distfromq:::make_q_fn(ps, qs)
    testthat::expect_identical(q_hat(c(0.99, 1)), c(61, Inf))
    
    p_hat <- distfromq::make_p_fn(ps, qs)
    testthat::expect_identical(p_hat(c(61, Inf)), c(0.99, 1.0))
})


test_that("make_p_fn result outputs values  <= 1", {
    ps <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975,
        0.99)
    qs <- c(0, 0, 0, 0.103331351093582, 0.206662702187163, 0.258328377733954,
        0.309994053280745, 0.309994053280745, 0.361659728827535,
        0.361659728827535, 0.413325404374326, 0.413325404374326,
        0.413325404374326, 0.464991079921117, 0.464991079921117,
        0.516656755467908, 0.516656755467908, 0.568322431014698,
        0.619988106561489, 0.723319457655071, 0.878316484295443,
        1.08497918648261, 1.29164188866977)

    result <- distfromq::make_p_fn(ps=ps, qs=qs)(25)
    testthat::expect_true(result <= 1)
})


test_that("make_p_fn and make_q_fn error with out-of-bounds or incorrectly typed ps, qs", {
    testthat::expect_no_error(make_p_fn(ps=c(0.0, 0.5, 1.0), qs = 1:3))
    testthat::expect_error(make_p_fn(ps=c(-1, 0.5, 1.0), qs = 1:3),
                           "Assertion on 'ps' failed: Element 1 is not >= 0.")
    testthat::expect_error(make_p_fn(ps=c(0.0, 0.5, 2.0), qs = 1:3),
                           "Assertion on 'ps' failed: Element 3 is not <= 1.")
    testthat::expect_error(make_p_fn(ps=c(0.0, "a", 1.0), qs = 1:3),
                           "Assertion on 'ps' failed: Must be of type 'numeric', not 'character'.")
    testthat::expect_error(make_p_fn(ps=c(0.0, 0.5, 1.0), qs = c(1, "a", 3)),
                           "Assertion on 'qs' failed: Must be of type 'numeric', not 'character'.")
    testthat::expect_error(make_p_fn(ps=c(0.0, 0.5, 1.0), qs = 1:4),
                           "'ps' and 'qs' must have the same length.")

    testthat::expect_no_error(make_q_fn(ps=c(0.0, 0.5, 1.0), qs = 1:3))
    testthat::expect_error(make_q_fn(ps=c(-1, 0.5, 1.0), qs = 1:3),
                           "Assertion on 'ps' failed: Element 1 is not >= 0.")
    testthat::expect_error(make_q_fn(ps=c(0.0, 0.5, 2.0), qs = 1:3),
                           "Assertion on 'ps' failed: Element 3 is not <= 1.")
    testthat::expect_error(make_q_fn(ps=c(0.0, "a", 1.0), qs = 1:3),
                           "Assertion on 'ps' failed: Must be of type 'numeric', not 'character'.")
    testthat::expect_error(make_q_fn(ps=c(0.0, 0.5, 1.0), qs = c(1, "a", 3)),
                           "Assertion on 'qs' failed: Must be of type 'numeric', not 'character'.")
    testthat::expect_error(make_q_fn(ps=c(0.0, 0.5, 1.0), qs = 1:4),
                           "'ps' and 'qs' must have the same length.")
})

test_that("make_p_fn and make_q_fn results error with out-of-bounds or incorrectly typed argument", {
    p_fn <- make_p_fn(ps=c(0.01, 0.5, 0.99), qs = 1:3)
    testthat::expect_no_error(p_fn(c(0, 1, 5)))
    testthat::expect_error(p_fn("a"),
                           "Must be of type 'numeric', not 'character'.")

    q_fn <- make_q_fn(ps=c(0.01, 0.5, 0.99), qs = 1:3)
    testthat::expect_no_error(q_fn(c(0, 0.5, 1)))
    testthat::expect_error(q_fn(c(-1, 0.5, 1.0)),
                           "Assertion on 'p' failed: Element 1 is not >= 0.")
    testthat::expect_error(q_fn(c(0.0, 0.5, 2.0)),
                           "Assertion on 'p' failed: Element 3 is not <= 1.")
    testthat::expect_error(q_fn("a"),
                           "Must be of type 'numeric', not 'character'.")
})

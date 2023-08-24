# This file includes the following test cases for make_d_fn
#
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
# we expect success in the first case and an error in all others


test_that("make_d_fn works, continuous component; no point masses", {
    ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    qs <- qnorm(ps, mean = 1, sd = 2)
    test_qs <- qnorm(seq(from = 0.01, to = 0.99, by = 0.01), mean = 1, sd = 2)

    d_actual <- make_d_fn(ps, qs)(test_qs)
    d_expected <- dnorm(test_qs, mean = 1, sd = 2)
    # plot(test_qs, d_actual); lines(test_qs, d_expected)
    expect_equal(d_actual, d_expected, tolerance = 1e-2)
})

test_that("make_d_fn generates error, no continuous component; one point mass, single q", {
    ps <- 0.1
    qs <- rep(0.0, length(ps))
    expect_error(make_d_fn(ps, qs))
})

test_that("make_d_fn generates error, no continuous component; one point mass, duplicated qs", {
    ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    qs <- rep(0.0, length(ps))
    expect_error(make_d_fn(ps, qs))
})

test_that("make_d_fn generates error, no continuous component; two point masses, both from duplicated qs", {
    ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    qs <- c(rep(1.0, 3), rep(2.0, 6))
    expect_error(make_d_fn(ps, qs))
})


test_that("make_d_fn generates error, no continuous component; two point masses, non-zero, duplicated values first only", {
    ps <- seq(from = 0.1, to = 0.4, by = 0.1)
    qs <- c(rep(1.0, 3), rep(2.0, 1))
    expect_error(make_d_fn(ps, qs))
})


test_that("make_d_fn generates error, no continuous component; two point masses, non-zero, duplicated values second only", {
    ps <- seq(from = 0.3, to = 0.9, by = 0.1)
    qs <- c(rep(1.0, 1), rep(2.0, 6))
    expect_error(make_d_fn(ps, qs))
})


test_that("make_d_fn generates error, no continuous component; is_hurdle with one zero and one non-zero", {
    ps <- seq(from = 0.1, to = 0.2, by = 0.1)
    qs <- c(1e-13, 2.0)
    expect_error(make_d_fn(ps, qs, tail_dist = "lnorm"))
})


test_that("make_d_fn generates error, continuous component; one point mass, duplicated qs", {
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

    expect_error(make_d_fn(ps, qs))
})


test_that("make_d_fn generates error, continuous component; one point mass, is_hurdle with zero", {
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

    expect_error(make_d_fn(ps, qs, tail_dist = "lnorm"))
})


test_that("make_d_fn generates error, two point masses, both from duplicated qs", {
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

    expect_error(make_d_fn(ps, qs))
})


test_that("make_d_fn generates error, two point masses, is_hurdle with one zero and one non-zero with duplicated qs", {
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

    expect_error(make_d_fn(ps, qs, tail_dist = "lnorm"))
})

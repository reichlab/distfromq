test_that("split_disc_cont_ps_qs works, no discrete component", {
    ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    qs <- qnorm(ps)
    expect_equal(
        split_disc_cont_ps_qs(ps, qs),
        list(
            disc_weight = 0.0,
            disc_ps = numeric(), disc_qs = numeric(),
            cont_ps = ps, cont_qs = qs,
            disc_ps_range = list()
        )
    )
})

test_that("split_disc_cont_ps_qs works, no continuous component", {
    ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    qs <- rep(0.0, length(ps))
    expect_equal(
        split_disc_cont_ps_qs(ps, qs),
        list(
            disc_weight = 1.0,
            disc_ps = 1.0, disc_qs = 0.0,
            cont_ps = numeric(), cont_qs = numeric(),
            disc_ps_range = list(c(0.1, 0.9))
        )
    )
})

test_that("split_disc_cont_ps_qs works, one discrete component", {
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

    expect_equal(
      split_disc_cont_ps_qs(ps, qs),
      list(
        disc_weight = 0.2,
        disc_ps = 1.0,
        disc_qs = 0.0,
        cont_ps = norm_ps,
        cont_qs = norm_qs,
        disc_ps_range = list(c(0.4, 0.6))
      )
    )
})

test_that("split_disc_cont_ps_qs works, two discrete components", {
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

    expect_equal(
      split_disc_cont_ps_qs(ps, qs),
      list(
        disc_weight = 0.4,
        disc_ps = c(0.75, 0.25),
        disc_qs = c(0.0, 1.0),
        cont_ps = sort(c(norm_ps, pnorm(1.0))),
        cont_qs = sort(c(norm_qs, 1.0)),
        disc_ps_range = list(
          range(adj_point_ps_0),
          range(adj_point_ps_1)
        )
      )
    )
})


test_that("split_disc_cont_ps_qs fails, one discrete component mismatched ps", {
    # mixture of a Normal(0,1) with weight 0.8 and
    # a point mass at 0 with weight 0.2

    # probabilities and quantiles for normal component
    norm_ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    norm_ps <- norm_ps[norm_ps != 0.5]
    norm_qs <- qnorm(norm_ps)
    adj_norm_ps <- norm_ps * 0.8 + 0.2 * (norm_qs > 0.0)

    # probabilities and quantiles for point mass at 0
    point_ps <- seq(from = 0.2, to = 1.0, by = 0.1)
    point_qs <- rep(0.0, length(point_ps))
    adj_point_ps <- 0.5 * 0.8 + point_ps * 0.2

    ps <- sort(c(adj_norm_ps, adj_point_ps))
    qs <- sort(c(norm_qs, point_qs))
    dup_inds <- duplicated(ps)
    ps <- ps[!dup_inds]
    qs <- qs[!dup_inds]

    # NOTE: ps does not include 0.4, the value at which we switch
    # between the continuous and discrete distributions. This
    # results in an underestimate of the weight of the discrete
    # component.

    # expect_equal(
    #   split_disc_cont_ps_qs(ps, qs),
    #   list(
    #     disc_weight = 0.2,
    #     disc_ps = 1.0,
    #     disc_qs = 0.0,
    #     cont_ps = norm_ps,
    #     cont_qs = norm_qs
    #   )
    # )
})


test_that("spline_cdf recovers cdf, no discrete component", {
    ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    qs <- qnorm(ps)

    cdf_hat <- spline_cdf(ps = ps, qs = qs, fn_type = "p",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm")
    cdf_hat_lin <- spline_cdf(ps = ps, qs = qs, fn_type = "p",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm",
                          n_grid = 2L)

    test_qs <- seq(from = min(qs), to = max(qs), length.out = 101)
    test_ps <- pnorm(test_qs)
    test_p_hats <- cdf_hat(test_qs)
    test_p_hats_lin <- cdf_hat_lin(test_qs)

    expect_equal(test_ps, test_p_hats, tolerance = 1e-3)
    expect_equal(test_p_hats_lin, test_p_hats, tolerance = 1e-3)
    expect_equal(mean(test_ps - test_p_hats), 0.0)
})


test_that("spline_cdf recovers cdf, one discrete component", {
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

    cdf_hat <- spline_cdf(ps = ps, qs = qs, fn_type = "p",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm")
    cdf_hat_lin <- spline_cdf(ps = ps, qs = qs, fn_type = "p",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm",
                          n_grid = 2L)

    test_qs <- seq(from = min(qs), to = max(qs), length.out = 101)
    test_ps <- pnorm(test_qs) * 0.8 + 0.2 * (test_qs >= 0.0)
    test_p_hats <- cdf_hat(test_qs)
    test_p_hats_lin <- cdf_hat_lin(test_qs)

    expect_equal(test_ps, test_p_hats, tolerance = 1e-3)
    expect_equal(mean(test_ps - test_p_hats), 0.0)
    expect_equal(test_p_hats_lin, test_p_hats, tolerance = 1e-3)
})

test_that("spline_cdf recovers cdf, two discrete components", {
    # mixture of a Normal(0,1) with weight 0.6,
    # a point mass at 0 with weight 0.3, and a point mass at 1 with weight 0.1

    # probabilities and quantiles for normal component
    norm_ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    norm_qs <- qnorm(norm_ps)
    adj_norm_ps <- norm_ps * 0.6 + 0.3 * (norm_qs >= 0.0) + 0.1 * (norm_qs >= 1.0)

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

    cdf_hat <- spline_cdf(ps = ps, qs = qs, fn_type = "p",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm")
    cdf_hat_lin <- spline_cdf(ps = ps, qs = qs, fn_type = "p",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm",
                          n_grid = 2L)

    test_qs <- seq(from = min(qs), to = max(qs), length.out = 101)
    test_ps <- pnorm(test_qs) * 0.6 + 0.3 * (test_qs >= 0.0) +
                0.1 * (test_qs >= 1.0)
    test_p_hats <- cdf_hat(test_qs)
    test_p_hats_lin <- cdf_hat_lin(test_qs)

    expect_equal(test_ps, test_p_hats, tolerance = 1e-3)
    expect_equal(test_p_hats_lin, test_p_hats, tolerance = 1e-3)
    expect_equal(mean(test_ps - test_p_hats), 0.0, tolerance = 1e-5)
})



test_that("spline_cdf recovers pdf, no discrete component", {
    ps <- seq(from = 0.1, to = 0.9, by = 0.01)
    qs <- qnorm(ps)

    pdf_hat <- spline_cdf(ps = ps, qs = qs, fn_type = "d",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm")
    pdf_hat_lin <- spline_cdf(ps = ps, qs = qs, fn_type = "d",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm",
                          n_grid = 20L)

    test_qs <- seq(from = min(qs), to = max(qs), length.out = 101)
    test_ps <- dnorm(test_qs)
    test_p_hats <- pdf_hat(test_qs)
    test_p_hats_lin <- pdf_hat_lin(test_qs)

    expect_equal(test_ps, test_p_hats, tolerance = 1e-3)
    expect_equal(test_p_hats_lin, test_p_hats, tolerance = 1e-3)
    expect_equal(mean(test_ps - test_p_hats), 0.0, tolerance = 1e-5)
})

test_that("spline_cdf errors when recovering pdf, one discrete component", {
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

    expect_error(spline_cdf(ps = ps, qs = qs, fn_type = "d",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm"))
    expect_error(spline_cdf(ps = ps, qs = qs, fn_type = "d",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm",
                          n_grid = 2L))
})




test_that("spline_cdf recovers qf, no discrete component", {
    ps <- c(0.01, 0.025,
            seq(from = 0.05, to = 0.95, by = 0.05), 0.975, 0.99)
    qs <- qnorm(ps)

    qf_hat <- spline_cdf(ps = ps, qs = qs, fn_type = "q",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm")
    qf_hat_lin <- spline_cdf(ps = ps, qs = qs, fn_type = "q",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm",
                          n_grid = 10L)

    cdf_hat <- spline_cdf(ps = ps, qs = qs, fn_type = "p",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm")
    cdf_hat_lin <- spline_cdf(ps = ps, qs = qs, fn_type = "p",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm",
                          n_grid = 10L)

    test_ps <- seq(from = min(ps), to = max(ps), length.out = 101)
    test_qs <- qnorm(test_ps)
    test_q_hats <- qf_hat(test_ps)
    test_q_hats_lin <- qf_hat_lin(test_ps)

    expect_equal(test_qs, test_q_hats, tolerance = 1e-3)
    expect_equal(test_q_hats_lin, test_q_hats, tolerance = 1e-3)
    expect_equal(mean(test_qs - test_q_hats), 0.0)
    expect_equal(cdf_hat(qf_hat(test_ps)), test_ps, tolerance = 1e-4)
    expect_equal(cdf_hat_lin(qf_hat_lin(test_ps)), test_ps, tolerance = 1e-12)
})


test_that("spline_cdf recovers qf, one discrete component", {
    # mixture of a Normal(0,1) with weight 0.8 and
    # a point mass at 0 with weight 0.2

    # probabilities and quantiles for normal component
    # norm_ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    norm_ps <- c(0.01, 0.025,
                 seq(from = 0.05, to = 0.95, by = 0.05), 0.975, 0.99)
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

    qf_hat <- spline_cdf(ps = ps, qs = qs, fn_type = "q",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm")
    qf_hat_lin <- spline_cdf(ps = ps, qs = qs, fn_type = "q",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm",
                          n_grid = 10L)

    cdf_hat <- spline_cdf(ps = ps, qs = qs, fn_type = "p",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm")
    cdf_hat_lin <- spline_cdf(ps = ps, qs = qs, fn_type = "p",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm",
                          n_grid = 10L)

    test_ps <- seq(from = min(ps), to = max(ps), length.out = 101)
    test_qs <- c(
      qnorm(test_ps[test_ps < 0.4] / 0.8),
      rep(0.0, sum((test_ps >= 0.4) & (test_ps <= 0.6))),
      qnorm((test_ps[test_ps > 0.6] - 0.2) / 0.8)
    )
    test_q_hats <- qf_hat(test_ps)
    test_q_hats_lin <- qf_hat_lin(test_ps)

    expect_equal(test_q_hats, test_qs, tolerance = 1e-3)
    expect_equal(test_q_hats, test_q_hats_lin, tolerance = 1e-3)
    expect_equal(mean(test_qs - test_q_hats), 0.0)
    expected_test_ps <- test_ps
    expected_test_ps[(test_ps >= 0.4) & (test_ps <= 0.6)] <- 0.6
    expect_equal(cdf_hat(qf_hat(test_ps)), expected_test_ps, tolerance = 1e-3)
    expect_equal(cdf_hat_lin(qf_hat_lin(test_ps)), expected_test_ps, tolerance = 1e-12)
})

test_that("spline_cdf recovers qf, two discrete components", {
    # mixture of a Normal(0,1) with weight 0.6,
    # a point mass at 0 with weight 0.3, and a point mass at 1 with weight 0.1

    # probabilities and quantiles for normal component
    # norm_ps <- seq(from = 0.1, to = 0.9, by = 0.1)
    norm_ps <- c(0.01, 0.025,
                 seq(from = 0.05, to = 0.95, by = 0.05), 0.975, 0.99)
    norm_qs <- qnorm(norm_ps)
    adj_norm_ps <- norm_ps * 0.6 + 0.3 * (norm_qs >= 0.0) + 0.1 * (norm_qs >= 1.0)

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

    qf_hat <- spline_cdf(ps = ps, qs = qs, fn_type = "q",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm")
    qf_hat_lin <- spline_cdf(ps = ps, qs = qs, fn_type = "q",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm",
                          n_grid = 10L)

    cdf_hat <- spline_cdf(ps = ps, qs = qs, fn_type = "p",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm")
    cdf_hat_lin <- spline_cdf(ps = ps, qs = qs, fn_type = "p",
                          lower_tail_dist = "norm",
                          upper_tail_dist = "norm",
                          n_grid = 10L)

    test_ps <- seq(from = min(ps), to = max(ps), length.out = 101)
    pcut1 <- 0.3
    pcut2 <- pnorm(1.0)*0.6 + 0.3
    test_qs <- c(
      qnorm(test_ps[test_ps < pcut1] / 0.6),
      rep(0.0, sum((test_ps >= 0.3) & (test_ps <= 0.6))),
      qnorm((test_ps[(test_ps > 0.6) & (test_ps < pcut2)] - 0.3) / 0.6),
      rep(1.0, sum((test_ps >= pcut2) & (test_ps <= pcut2 + 0.1))),
      qnorm((test_ps[test_ps > pcut2 + 0.1] - 0.4) / 0.6)
    )
    test_q_hats <- qf_hat(test_ps)
    test_q_hats_lin <- qf_hat_lin(test_ps)

    expect_equal(test_q_hats, test_qs, tolerance = 1e-3)
    expect_equal(test_q_hats, test_q_hats_lin, tolerance = 1e-3)
    expect_equal(mean(test_qs - test_q_hats), 0.0, tolerance = 1e-4)
    expected_test_ps <- test_ps
    expected_test_ps[(test_ps >= 0.3) & (test_ps <= 0.6)] <- 0.6
    expected_test_ps[(test_ps >= pcut2) & (test_ps <= pcut2 + 0.1)] <- pcut2 + 0.1
    expect_equal(cdf_hat(qf_hat(test_ps)), expected_test_ps, tolerance = 1e-3)
    expect_equal(cdf_hat_lin(qf_hat_lin(test_ps)), expected_test_ps, tolerance = 1e-12)
})

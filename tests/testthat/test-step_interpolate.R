test_that("step_interp works, no discrete component", {
  ps <- seq(from = 0.1, to = 0.9, by = 0.1)
  qs <- qnorm(ps)

  step_interp_r <- step_interp_factory(x = ps, y = qs, cont_dir = "right")
  step_interp_l <- step_interp_factory(x = ps, y = qs, cont_dir = "left")

  expect_equal(step_interp_r(ps[1:8]), qs[1:8])
  expect_equal(step_interp_l(ps[2:9]), qs[2:9])

  expect_equal(step_interp_r(c(0.17, 0.32)),
               c((1 - 0.7) * qs[1] + 0.7 * qs[2],
                 (1 - 0.2) * qs[3] + 0.2 * qs[4]))
  expect_equal(step_interp_l(c(0.17, 0.32)),
               c((1 - 0.7) * qs[1] + 0.7 * qs[2],
                 (1 - 0.2) * qs[3] + 0.2 * qs[4]))
})

test_that("step_interp works, no continuous component", {
  ps <- seq(from = 0.1, to = 0.9, by = 0.1)
  qs <- rep(1.0, length(ps))

  step_interp_r <- step_interp_factory(x = qs, y = ps, cont_dir = "right")
  step_interp_l <- step_interp_factory(x = qs, y = ps, cont_dir = "left")

  test_qs <- c(0.99999, 1.0, 1.00001)
  expect_equal(
    step_interp_r(test_qs),
    c(NA_real_, 0.9, NA_real_)
  )
  expect_equal(
    step_interp_l(test_qs),
    c(NA_real_, 0.1, NA_real_)
  )
})

test_that("step_interp works, one discrete component", {
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

  step_interp_cdf_r <- step_interp_factory(x = qs, y = ps, cont_dir = "right")
  step_interp_cdf_l <- step_interp_factory(x = qs, y = ps, cont_dir = "left")

  test_input_qs <- sort(c(range(qs), seq(from = min(qs) - 0.00001,
                                         to = max(qs) + 0.00001,
                                         length.out = 101)))

  # manual calculation for step_interp_cdf_r
  expected <- rep(NA_real_, length(test_input_qs))
  for (i in seq_len(length(norm_qs) - 1)) {
    qi <- norm_qs[i]
    qip1 <- norm_qs[i + 1]
    pi <- adj_norm_ps[i]
    pip1 <- adj_norm_ps[i + 1]
    if (qi == 0.0) {
      pi <- pi + 0.2
    }
    delta <- qip1 - qi
    j <- which(test_input_qs >= qi & test_input_qs < qip1)
    expected[j] <- pi + (test_input_qs[j] - qi) * (pip1 - pi) / (delta)
  }
  expected[test_input_qs == 0.0] <- max(ps[qs == 0.0])
  expected[test_input_qs == max(qs)] <- max(ps)
  expect_equal(step_interp_cdf_r(test_input_qs), expected, tolerance = 1e-8)

  # only change in expected output from step_interp_cdf_l is at the step
  expected[test_input_qs == 0.0] <- min(ps[qs == 0.0])
  expect_equal(step_interp_cdf_l(test_input_qs), expected, tolerance = 1e-8)
})

test_that("step_interp works, two discrete components", {
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

  step_interp_cdf_r <- step_interp_factory(x = qs, y = ps, cont_dir = "right")
  step_interp_cdf_l <- step_interp_factory(x = qs, y = ps, cont_dir = "left")

  test_input_qs <- sort(c(range(qs), seq(from = min(qs) - 0.00001,
                                         to = max(qs) + 0.00001,
                                         length.out = 101)))

  # manual calculation for step_interp_cdf_r
  expected <- rep(NA_real_, length(test_input_qs))
  norm_qs <- sort(c(norm_qs, 1.0))
  adj_norm_ps <- sort(c(adj_norm_ps, pnorm(1.0) * 0.6 + 0.3))
  for (i in seq_len(length(norm_qs) - 1)) {
    qi <- norm_qs[i]
    qip1 <- norm_qs[i + 1]
    pi <- adj_norm_ps[i]
    pip1 <- adj_norm_ps[i + 1]
    if (qi == 0.0) {
      pi <- pi + 0.3
    } else if (qi == 1.0) {
      pi <- pi + 0.1
    }
    delta <- qip1 - qi
    j <- which(test_input_qs >= qi & test_input_qs < qip1)
    expected[j] <- pi + (test_input_qs[j] - qi) * (pip1 - pi) / (delta)
  }
  expected[test_input_qs == 0.0] <- max(ps[qs == 0.0])
  expected[test_input_qs == 1.0] <- max(ps[qs == 1.0])
  expected[test_input_qs == max(qs)] <- max(ps)
  expect_equal(step_interp_cdf_r(test_input_qs), expected, tolerance = 1e-8)

  # only change in expected output from step_interp_cdf_l is at the steps
  expected[test_input_qs == 0.0] <- min(ps[qs == 0.0])
  expected[test_input_qs == 1.0] <- min(ps[qs == 1.0])
  expect_equal(step_interp_cdf_l(test_input_qs), expected, tolerance = 1e-8)
})

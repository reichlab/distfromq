# This file includes the following test cases for split_disc_cont_ps_qs
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
# additionally there is a "test" with a known (and unaddressable) failure in a
# setting where the ps do not coincide with the start of a point mass; in this
# case, there is not enough information in the input ps and qs to get the point
# mass probability right

test_that("split_disc_cont_ps_qs works, continuous component; no point masses", {
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

test_that("split_disc_cont_ps_qs works, no continuous component; one point mass, single q", {
  ps <- 0.1
  qs <- rep(0.0, length(ps))
  expect_equal(
    split_disc_cont_ps_qs(ps, qs),
    list(
      disc_weight = 1.0,
      disc_ps = 1.0, disc_qs = 0.0,
      cont_ps = numeric(), cont_qs = numeric(),
      disc_ps_range = list(c(0.0, 1.0))
    )
  )
})

test_that("split_disc_cont_ps_qs works, no continuous component; one point mass, duplicated qs", {
  ps <- seq(from = 0.1, to = 0.9, by = 0.1)
  qs <- rep(0.0, length(ps))
  expect_equal(
    split_disc_cont_ps_qs(ps, qs),
    list(
      disc_weight = 1.0,
      disc_ps = 1.0, disc_qs = 0.0,
      cont_ps = numeric(), cont_qs = numeric(),
      disc_ps_range = list(c(0.0, 1.0))
    )
  )
})

test_that("split_disc_cont_ps_qs works, no continuous component; two point masses, both from duplicated qs", {
  ps <- seq(from = 0.1, to = 0.9, by = 0.1)
  qs <- c(rep(1.0, 3), rep(2.0, 6))
  expect_equal(
    split_disc_cont_ps_qs(ps, qs),
    list(
      disc_weight = 1.0,
      disc_ps = c(3 / 9, 6 / 9), disc_qs = c(1.0, 2.0),
      cont_ps = numeric(), cont_qs = numeric(),
      disc_ps_range = list(c(0.0, 3 / 9), c(3 / 9, 1.0))
    )
  )
})


test_that("split_disc_cont_ps_qs works, no continuous component;
          two point masses, non-zero, duplicated values first only", {
            ps <- seq(from = 0.1, to = 0.4, by = 0.1)
            qs <- c(rep(1.0, 3), rep(2.0, 1))
            expect_equal(
              split_disc_cont_ps_qs(ps, qs),
              list(
                disc_weight = 1.0,
                disc_ps = c(3 / 9, 6 / 9), disc_qs = c(1.0, 2.0),
                cont_ps = numeric(), cont_qs = numeric(),
                disc_ps_range = list(c(0.0, 3 / 9), c(3 / 9, 1.0))
              )
            )
          })


test_that("split_disc_cont_ps_qs works, no continuous component;
          two point masses, non-zero, duplicated values second only", {
            ps <- seq(from = 0.3, to = 0.9, by = 0.1)
            qs <- c(rep(1.0, 1), rep(2.0, 6))
            expect_equal(
              split_disc_cont_ps_qs(ps, qs),
              list(
                disc_weight = 1.0,
                disc_ps = c(3 / 9, 6 / 9), disc_qs = c(1.0, 2.0),
                cont_ps = numeric(), cont_qs = numeric(),
                disc_ps_range = list(c(0.0, 3 / 9), c(3 / 9, 1.0))
              )
            )
          })


test_that("split_disc_cont_ps_qs works, no continuous component; is_hurdle with one zero and one non-zero", {
  ps <- seq(from = 0.1, to = 0.2, by = 0.1)
  qs <- c(1e-13, 2.0)
  expect_equal(
    split_disc_cont_ps_qs(ps, qs, is_hurdle = TRUE),
    list(
      disc_weight = 1.0,
      disc_ps = c(1 / 9, 8 / 9), disc_qs = c(1e-13 / 2, 2.0),
      cont_ps = numeric(), cont_qs = numeric(),
      disc_ps_range = list(c(0.0, 1 / 9), c(1 / 9, 1.0))
    )
  )
})


test_that("split_disc_cont_ps_qs works, continuous component; one point mass, duplicated qs", {
  # mixture of a Normal(0, 1) with weight 0.8 and
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


test_that("split_disc_cont_ps_qs works, continuous component; one point mass, is_hurdle with zero", {
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

  expect_equal(
    split_disc_cont_ps_qs(ps, qs, is_hurdle = TRUE),
    list(
      disc_weight = 0.2,
      disc_ps = 1.0,
      disc_qs = 0.0,
      cont_ps = norm_ps,
      cont_qs = norm_qs,
      disc_ps_range = list(c(0.0, 0.2))
    )
  )
})


test_that("split_disc_cont_ps_qs works, two point masses, both from duplicated qs", {
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


test_that(
  "split_disc_cont_ps_qs works, two point masses, is_hurdle with one zero and one non-zero with duplicated qs",
  {
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

    expect_equal(
      split_disc_cont_ps_qs(ps, qs, is_hurdle = TRUE),
      list(
        disc_weight = 0.4,
        disc_ps = c(0.75, 0.25),
        disc_qs = c(0.0, 1.0),
        cont_ps = sort(c(norm_ps)),
        cont_qs = sort(c(norm_qs)),
        disc_ps_range = list(
          c(0.0, 0.3),
          range(adj_point_ps_1)
        )
      )
    )
  }
)


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

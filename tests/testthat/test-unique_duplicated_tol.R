test_that("get_dup_run_inds, unique_tol and duplicated_tol work, first value duplicated", {
    x <- c(1, 1 + 1e-7, 2, 3, 4, 4, 4, 5, 6,
           7 - 2e-7, 7 - 1e-7, 7 + 1e-7, 7 + 2e-7, 8, 9)

    dxf <- duplicated_tol(x, incl_first = FALSE)
    dxt <- duplicated_tol(x, incl_first = TRUE)
    ux <- unique_tol(x)
    dri <- get_dup_run_inds(dxf)
    expect_true(all.equal(
        dxf,
        c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE,
          FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)))
    expect_true(all.equal(
        dxt,
        c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE,
          TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)))
    expect_true(all.equal(ux, c(1 + 5e-8, 2:9)))
    expect_equal(
        dri,
        list(starts = c(1, 5, 10), ends = c(2, 7, 13))
    )
})

test_that("get_dup_run_inds, unique_tol and duplicated_tol work, first value not duplicated", {
    x <- c(0, 1, 2, 3, 4, 4, 4, 5, 6,
           7 - 2e-7, 7 - 1e-7, 7 + 1e-7, 7 + 2e-7, 8, 9)

    dxf <- duplicated_tol(x, incl_first = FALSE)
    dxt <- duplicated_tol(x, incl_first = TRUE)
    ux <- unique_tol(x)
    dri <- get_dup_run_inds(dxf)
    expect_true(all.equal(
        dxf,
        c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE,
          FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)))
    expect_true(all.equal(
        dxt,
        c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE,
          TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)))
    expect_true(all.equal(ux, 0:9))
    expect_equal(
        dri,
        list(starts = c(5, 10), ends = c(7, 13))
    )
})


test_that("get_dup_run_inds, unique_tol and duplicated_tol work, consecutive runs of different duplicated values", {
    x <- c(rep(0, 5), rep(1, 4))

    dxf <- duplicated_tol(x, incl_first = FALSE)
    dxt <- duplicated_tol(x, incl_first = TRUE)
    ux <- unique_tol(x)
    dri <- get_dup_run_inds(dxf)
    expect_true(all.equal(
        dxf,
        c(FALSE, TRUE, TRUE, TRUE, TRUE,
          FALSE, TRUE, TRUE, TRUE)))
    expect_true(all.equal(
        dxt,
        c(TRUE, TRUE, TRUE, TRUE, TRUE,
          TRUE, TRUE, TRUE, TRUE)))
    expect_true(all.equal(ux, 0:1))
    expect_equal(
        dri,
        list(starts = c(1, 6), ends = c(5, 9))
    )
})

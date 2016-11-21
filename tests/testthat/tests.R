context("Test")

test_that("Generated xs have right dimensions", {
    expect_equal(dim(gen_xs_default(5, 2)), c(5, 1))
    expect_equal(dim(gen_xs_corr(5, 2)),    c(5, 1))
})

test_that("Generated ys have right dimensions", {
    n <- 100
    ys <- gen_data(n, c(1, 2))$ys
    expect_equal(length(ys), n)
})

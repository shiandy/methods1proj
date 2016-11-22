context("General Unit Tests")

test_that("Generated xs have right dimensions", {
    expect_equal(dim(gen_xs_default(5, 2)), c(5, 1))
    expect_equal(dim(gen_xs_corr(5, 2)),    c(5, 1))
})

test_that("Generated ys have right dimensions", {
    n <- 100
    ys <- gen_data(n, c(1, 2))$ys
    expect_equal(length(ys), n)
})

test_that("run_sim stops with wrong parameters", {
    expect_error(run_sim(1, 10, c(1, 2, 3), split_prop = 1.1))
})

test_that("Train/test have unique rows and correct # of rows", {
    xs <- data.frame(a = 1:10, b = 1:10)
    ys <- data.frame(a = 20:30, b = 20:30)
    split_prop <- 0.7
    train_test <- train_test_split(xs, ys, split_prop)
    xs_train <- train_test$xs_train
    xs_test <- train_test$xs_test
    expect_false(any(xs_train$a %in% xs_test$a))
    expect_equal(nrow(xs_train), floor(split_prop * nrow(xs)))
})

test_that("train_test_split stops when split_prop wrong", {
    expect_error(train_test_split(1, 1, -1))
    expect_error(train_test_split(1, 1, 1))
})

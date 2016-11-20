gen_xs_default <- function(n, p) {
    xs <- matrix(rnorm(n * (p - 1)), nrow = n, ncol = p - 1)
    return(xs)
}

gen_data <- function(n, true_betas, gen_dist = gen_xs_default,
                     error_dist = rnorm) {
    p <- length(true_betas)
    xs <- gen_dist(n, p)
    errors <- error_dist(n)
    ys <- true_betas[1] + xs %*% true_betas[-1] + errors
    return(list(xs = as.data.frame(xs), ys = ys))
}

pick_fit_model_leaps <- function(xs, ys, method) {
    res <- leaps::leaps(xs, ys, method = method)
    best_mod <- which.max(get(method, res))
    which_coefs <- res$which[best_mod, ]
    best_fit <- lm(ys ~ ., data = xs[which_coefs])
    return(confint(best_fit))
}

pick_fit_model_step <- function(xs, ys, method) {

}

pick_fit_model <- function(xs, ys, method = "leaps") {
    method <- match.arg(method)
    if (method == "leaps") {
        if(dim(xs)[2] > 6) {
            stop("Dimensions too high for leaps")
        }
        else {
            pick_fit_model_leaps(xs, ys, "adjr2")
        }
    }
    # return model and coefficients
}

run_sim <- function(nreps, true_betas) {

}

plot_sim <- function(true_betas, beta_cis) {
    # make a plot
}

train_test_split <- function(xs, split_prop = 0.5) {
    # return splitted data
}

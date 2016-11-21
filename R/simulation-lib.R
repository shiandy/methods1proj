gen_xs_default <- function(n, p) {
    xs <- matrix(rnorm(n * (p - 1)), nrow = n, ncol = p - 1)
    return(xs)
}

gen_xs_corr <- function(n, p, covar = 0.5) {
    mus <- rep(0, p - 1)
    covar_mat <- matrix(covar, nrow = p - 1, ncol = p - 1)
    diag(covar_mat) <- 1
    xs <- MASS::mvrnorm(n = n, mus, covar_mat)
    return(xs)
}

gen_data <- function(n, true_betas, gen_dist = gen_xs_default,
                     error_dist = rnorm, ...) {
    p <- length(true_betas)
    xs <- gen_dist(n, p, ...)
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

pick_fit_model_step <- function(xs, ys, direction) {
    if (direction != "backward") {
        fit_null <- lm(ys ~ 1, data = xs)
        formula_full <- paste("~", paste(colnames(xs),
                                         collapse = " + "))
        step_res <- step(fit_null, as.formula(formula_full),
                         direction = direction,
                         trace = 0)
    }
    else {
        fit_full <- lm(ys ~ ., data = xs)
        step_res <- step(fit_full, direction = "backward",
                         trace = 0)
    }
    return(confint(step_res))
}

pick_fit_model <- function(xs, ys, method = c("step", "leaps"),
                           direction = "both") {
    method <- match.arg(method)
    if (method == "leaps") {
        if(dim(xs)[2] > 6) {
            stop("Dimensions too high for leaps")
        }
        else {
            ci <- pick_fit_model_leaps(xs, ys, "adjr2")
        }
    }
    else {
        ci <- pick_fit_model_step(xs, ys, direction)
    }
    return(ci)
}

run_sim <- function(nreps, n, true_betas,
                    select_method = c("step", "leaps"),
                    direction = "both",
                    gen_dist = gen_xs_default, error_dist = rnorm,
                    ...) {
    sim_num_lst <- list()
    ci_lst <- list()
    for (i in 1:nreps) {
        dat <- gen_data(n, true_betas, gen_dist, error_dist, ...)
        ci <- pick_fit_model(dat$xs, dat$ys, select_method, direction)
        ci_lst[[i]] <- ci
        sim_num_lst[[i]] <- rep(i, dim(ci)[1])
    }

    ci_mat <- do.call(rbind, ci_lst)
    ret_df <- data.frame(lb = ci_mat[, 1], ub = ci_mat[, 2],
                         coef_name = rownames(ci_mat),
                         sim_num = unlist(sim_num_lst))
    return(ret_df)
}


plot_sim <- function(true_betas, beta_cis) {
    # make a plot
}

train_test_split <- function(xs, split_prop = 0.5) {
    # return splitted data
}

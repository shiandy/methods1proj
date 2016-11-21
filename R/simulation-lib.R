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
    # point estimate is midpoint
    ret_df$pt_est <- ret_df$lb + (ret_df$ub - ret_df$lb) / 2
    ret_df$signif <- !(ret_df$lb < 0 & 0 < ret_df$ub)
    return(ret_df)
}


plot_sim <- function(res_df, true_betas) {
    p <- length(true_betas)
    plots_lst <- list()
    coverages <- rep(NA, p)
    for (i in 1:p) {
        true_val <- true_betas[i]
        beta_name <- names(true_betas)[i]
        cur_dat <- res_df[res_df$coef_name == beta_name, ]
        cur_dat <- cur_dat[order(cur_dat$pt_est),]
        cur_dat$index <- 1:nrow(cur_dat)
        coverage <- mean(cur_dat$lb < true_val & true_val < cur_dat$ub)
        coverages[i] <- coverage

        plt <- ggplot2::ggplot(data = cur_dat,
                               ggplot2::aes(x = index, y = pt_est)) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = lb, ymax = ub),
                                   alpha = 0.5) +
            ggplot2::geom_hline(yintercept = true_val, color = "red") +
            ggplot2::ggtitle(paste(beta_name,
                                   "Coverage: ", round(coverage, 4)))
        plots_lst[[i]] <- plt
    }
    return(list(plots = plots_lst, coverages = coverages))
}

train_test_split <- function(xs, split_prop = 0.5) {
    # return splitted data
}

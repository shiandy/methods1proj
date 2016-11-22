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
        # direction is ignored for leaps
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
                    direction = "both", split_prop = -1,
                    gen_dist = gen_xs_default, error_dist = rnorm,
                    ...) {
    sim_num_lst <- list()
    ci_lst <- list()
    for (i in 1:nreps) {
        dat <- gen_data(n, true_betas, gen_dist, error_dist, ...)
        if (split_prop > 0) {
            stopifnot(split_prop < 1)
            train_test <- train_test_split(dat$xs, dat$ys, split_prop)
            picked_ci <- pick_fit_model(train_test$xs_train,
                                        train_test$ys_train,
                                        select_method, direction)
            #picked_form <- extract_form(picked_ci)
            vars <- rownames(picked_ci)[-1]
            form_str <- paste("train_test$ys_test", "~",
                              paste(vars, collapse = " + "))
            picked_form <- as.formula(form_str)
            test_fit <- lm(picked_form, data = train_test$xs_test)
            ci <- confint(test_fit)
            ci_lst[[i]] <- ci
        }
        else {
            ci <- pick_fit_model(dat$xs, dat$ys, select_method,
                                 direction)
            ci_lst[[i]] <- ci
        }
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

get_coverages <- function(res_df, true_betas) {
    p <- length(true_betas)
    coverages <- rep(NA, p)
    names(coverages) <- names(true_betas)
    for (i in 1:p) {
        true_val <- true_betas[i]
        beta_name <- names(true_betas)[i]
        cur_dat <- res_df[res_df$coef_name == beta_name, ]
        if (nrow(cur_dat) == 0) {
            next
        }
        coverage <- mean(cur_dat$lb < true_val & true_val < cur_dat$ub)
        coverages[i] <- coverage
    }
    return(coverages)
}

# helper function to run the simulation and generate plots
run_coverage_selected <- function(nreps, n, true_betas,
                                  select_method = c("step", "leaps"),
                                  direction = "both", split_prop = -1,
                                  gen_dist = gen_xs_default,
                                  error_dist = rnorm, ...) {

    p <- length(true_betas)
    start <- Sys.time()
    sim_res <- run_sim(nreps, n, true_betas,
                       select_method, direction, split_prop, gen_dist,
                       error_dist, ...)
    elapsed <- Sys.time() - start
    print(elapsed)
    names(true_betas) <- c("(Intercept)", paste0("V", 1:(p - 1)))

    print("Number of times selected:")
    print(summary(sim_res$coef_name))

    #sim_res %>% group_by(coef_name) %>% summarize(pct = sum(signif))

    coverages <- get_coverages(sim_res, true_betas)

    selected_df <- sim_res %>% dplyr::group_by(coef_name) %>%
            dplyr::summarize(pct_selected = n() / nreps)

    select_probs <- rep(NA, p)
    names(select_probs) <- names(true_betas)
    selected_names <- selected_df$coef_name
    select_probs[selected_names] <- selected_df$pct_selected
    select_probs[is.na(select_probs)] <- 0
    #stopifnot(nrow(selected_df) == p)
    #stopifnot(selected_df$coef_name == names(coverages))
    #select_probs <- selected_df$pct_selected
    #names(select_probs) <- names(true_betas)

    retval <- list(coverages = coverages, select_probs = select_probs)
    return(retval)
}


plot_sim <- function(res_df, true_betas) {
    p <- length(true_betas)
    plots_lst <- list()
    coverages <- rep(NA, p)
    names(coverages) <- names(true_betas)
    for (i in 1:p) {
        true_val <- true_betas[i]
        beta_name <- names(true_betas)[i]
        cur_dat <- res_df[res_df$coef_name == beta_name, ]
        if (nrow(cur_dat) == 0) {
            next
        }
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
                                   "Coverage: ", round(coverage, 4))) +
            ggplot2::xlab("Index") + ggplot2::ylab("95% CI")
        plots_lst[[i]] <- plt
    }
    return(list(plots = plots_lst, coverages = coverages))
}

train_test_split <- function(xs, ys, split_prop = 0.5) {
    stopifnot(split_prop > 0, split_prop < 1)
    samp_size <- floor(split_prop * nrow(xs))
    train_ind <- sample(seq_len(nrow(xs)), samp_size)
    xs_train <- xs[train_ind, ]
    xs_test <- xs[-train_ind, ]
    ys_train <- ys[train_ind, ]
    ys_test <- ys[-train_ind, ]
    ret <- list(xs_train = xs_train, xs_test = xs_test,
                ys_train = ys_train, ys_test = ys_test)
    return(ret)
}

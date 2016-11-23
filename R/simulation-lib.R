#' Default function to generate covariates X.
#'
#' Generates the X_i independently, with each component of X_i
#' also being an independent Normal(0, 1). Exported to the user.
#'
#' @param n Number of observations to generate
#' @param p Number of dimensions for each observation (includes the
#'  intercept)
#'
#' @return An n * (p-1) matrix of generated covariates
#' @export
#'
#' @examples
#' gen_xs_default(10, 2)
gen_xs_default <- function(n, p) {
    xs <- matrix(rnorm(n * (p - 1)), nrow = n, ncol = p - 1)
    return(xs)
}

#' Generate correlated covariates X
#'
#' Generate the X_i independently, but with each component of X_i
#' distributed Multivariate Normal with mean 0 and covariance matrix
#' with 1 on the diagonal and \code{covar} on the off-diagonal. Exported
#' to the user.
#'
#' @param covar Off-diagonal element of covariance matrix.
#' @inheritParams gen_xs_default
#'
#' @return An n * (p-1) matrix of generated covariates
#' @export
#'
#' @examples
#' gen_xs_corr(10, 2)
#' gen_xs_corr(10, 2, 0.5)
#'
#' set.seed(1)
#' # these should be equivalent
#' x1 <- gen_xs_corr(10, 2, 0)
#' set.seed(1)
#' x2 <- gen_xs_default(10, 2)
#' all(x1 == x2)
gen_xs_corr <- function(n, p, covar = 0.5) {
    mus <- rep(0, p - 1)
    covar_mat <- matrix(covar, nrow = p - 1, ncol = p - 1)
    diag(covar_mat) <- 1
    xs <- MASS::mvrnorm(n = n, mus, covar_mat)
    return(xs)
}

#' Generate both the covariates and the outcomes
#'
#' Generates X using \code{gen_dist} and then generates
#' Y = beta_0 + XB + eps, #' where eps is generated according to
#' \code{error_dist}. NOT exported to the user.
#'
#' @param n Number of observations
#' @inheritParams run_sim
#'
#'
#' @return
#' A list \code{lst}, with element \code{lst$xs} corresponding to the
#' covariates and \code{lst$ys} corresponding to the responses.
gen_data <- function(n, true_betas, gen_dist = gen_xs_default,
                     error_dist = rnorm, ...) {
    p <- length(true_betas)
    xs <- gen_dist(n, p, ...)
    errors <- error_dist(n)
    ys <- true_betas[1] + xs %*% true_betas[-1] + errors
    return(list(xs = as.data.frame(xs), ys = ys))
}

#' Pick the best model using \code{leaps} and return its fitted
#' confidence intervals.
#'
#' Pick the best model using \code{leaps} and return its fitted
#' confidence intervals. NOT exported to the user.
#'
#' @param xs covariates
#' @param ys responses
#' @param method A valid argument to the \code{methods} argument of the
#' leaps command
#'
#' @return 95% confidence intervals from the selected model.
pick_fit_model_leaps <- function(xs, ys, method) {
    res <- leaps::leaps(xs, ys, method = method)
    best_mod <- which.max(get(method, res))
    which_coefs <- res$which[best_mod, ]
    best_fit <- lm(ys ~ ., data = xs[which_coefs])
    return(confint(best_fit))
}

#' Pick the best model using \code{step} and return its fitted
#' confidence intervals.
#'
#' Pick the best model using \code{step} and return its fitted
#' confidence intervals. NOT exported to the user.
#'
#' @param xs covariates
#' @param ys responses
#' @param direction A valid argument to the \code{direction} argument of
#' the \code{step} command
#'
#' @return 95% confidence intervals from the selected model.
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

#' Pick the best linear regression model and fit it.
#'
#' Calls \code{pick_fit_model_leaps} or \code{pick_fit_model_step} to
#' actually do the underlying work. NOT exported to the user.
#'
#' @param xs Covariates
#' @param ys Responses
#' @param method If leaps, run leaps. Otherwise, if step, run step-wise
#' selection
#' @param direction Only valid of \code{method = "step"}. Run step-wise
#' selection using \code{step} in the specified direction
#'
#' @return 95% confidence intervals from the selected model.
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

#' Run one replicate of the simulation.
#'
#' Run one replicate of the simulation. Exported to the user.
#'
#' @param nreps Number of repetitions to run, or how many datasets to
#' generate
#' @param n Number of observations per dataset.
#' @param select_method If leaps, run leaps. Otherwise, if step, run
#' step-wise selection
#' @param direction Only valid of \code{method = "step"}. Run step-wise
#' selection using \code{step} in the specified direction
#' @param split_prop If non-positive, run naive method, where the
#' confidence #' intervals are calculated on the same data used to
#' select the variables. Otherwise, if 0 < \code{split_prop} < 1,
#' randomly subset \code{split_prop} of the data to be in the training
#' set, and the other part to be in the test set. Run variable selection
#' on the training set, and report the final confidence interval as
#' calculated by fitting the selected model on the test set.
#' @param true_betas A vector of true values of regression coefficients.
#' Should include the intercept as the first item, beta1 as the second,
#' etc.
#' @param gen_dist A function specifiying how to generate the
#' covariates. Should take \code{n} as its first argument, corresponding
#' to the number of observations, and \code{p} as the second,
#' corresponding to the number of dimensions per observation (including
#' the intercept, so each X_i has p-1 dimensions).
#' @param error_dist A function specifying how to generate the errors.
#' Should take \code{n} as the first argument, specifying the number of
#' observations
#' @param ... Additional arguments to gen_dist
#' @return A dataframe, with each row corresponding to one coefficient
#' that has been selected in one iteration, and columns:
#'
#' \item{\code{lb}}{Lower bound of the confidence interval}
#' \item{\code{ub}}{Upper bound of the confidence interval}
#' \item{\code{coef_name}}{Name of the coefficient}
#' \item{\code{sim_num}}{Which iteration of the simulation is this
#' coefficient from?}
#'
#' @export
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
            # extract names of CI and build a formula to refit the
            # model on the test set
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
        # keep track of which simulation this was
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

#' Get the coverage of each coefficient
#'
#' @param res_df data.frame returned by \code{\link{run_sim}}
#' @inheritParams run_sim
#'
#' @return A vector of coverages, one per each coefficient. The vector
#'   is nameded, with each element corresponding to the coverage
#'   probability of one coefficient in one iteration of the simulation
#'   (is NA if that coefficient was never selected).
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

#' Obtain coverage and selection probabilities
#'
#' Calls \link{\code{run_sim}} once and processes the result to obtain
#'  coverage and selection probabilities.
#'
#' Might be a bug with \code{get_coverage} and this function. Run with
#' n = 500 and split_prob = 0.8 and they're all NA.
#'
#' @inheritParams run_sim
#' @return A list \code{lst} containing two vectors. The vector
#' \code{lst$coverages} is named, with each element corresponding to
#' the coverage probability of one coefficient in one iteration of the
#' simulation (is NA if that coefficient was never selected). The
#' vector \code{lst$select_probs} is like \code{lst$coverages}, but
#' instead contains the selection probabilities.
#' @export
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

    stopifnot(all(names(coverages) == names(select_probs)))

    retval <- list(coverages = coverages, select_probs = select_probs)
    return(retval)
}


#' Draws coverage plots for one run of the simulation.
#'
#' Not really recommended because you ideally want to run the simulation
#' multiple times to get a sense of the variation. NOT exported to the
#' user
#'
#' @param res_df A data.frame returned by \link{\code{run_sim}}.
#' @param true_betas A vector of true values of beta, including the
#' intercept
#'
#' @return A list of ggplot2 objects, one per coefficient
plot_sim <- function(res_df, true_betas) {
    p <- length(true_betas)
    plots_lst <- list()
    beta_names_vec <- unique(res_df$coef_name)
    for (i in 1:p) {
        true_val <- true_betas[i]
        beta_name <- beta_names_vec[i]
        cur_dat <- res_df[res_df$coef_name == beta_name, ]
        if (nrow(cur_dat) == 0) {
            next
        }
        cur_dat <- cur_dat[order(cur_dat$pt_est),]
        cur_dat$index <- 1:nrow(cur_dat)
        coverage <- mean(cur_dat$lb < true_val & true_val < cur_dat$ub)

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
    return(plots_lst)
}

#' Train-test split
#'
#' Split the data randomly into a training and test set.
#' NOT exported to the user
#'
#' @param xs Covariates
#' @param ys Responses
#' @param split_prop Proportion in the training set. Must be between 0
#' and 1, non-inclusive
#'
#' @return A list with elements
#' \item{\code{xs_train}}{Training set of covariates}
#' \item{\code{xs_test}}{Test set of covariates}
#' \item{\code{ys_train}}{Training set of responses}
#' \item{\code{ys_test}}{Test set of responses}
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

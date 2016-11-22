library("methods1proj")
library("ggplot2")
library("cowplot")
library("doParallel")
library("foreach")

run_sim_parallel <- function(true_betas, sim_conditions, nreps = 10) {
    all_data <- foreach(i = 1:nrow(sim_conditions)) %dopar% {
        require("methods1proj")
        require("magrittr")
        require("dplyr")
        set.seed(1262 + i)
        p <- length(true_betas)
        cur_n <- sim_conditions$n[i]
        cur_select <- sim_conditions$select_method[i]
        cur_split_prop <- sim_conditions$split_prop[i]
        cur_covar <- sim_conditions$covar[i]
        cur_trial <- sim_conditions$trial_index[i]

        if (cur_select != "leaps") {
            select_method <- "step"
        } else {
            select_method <- "leaps"
        }

        # essentially, covar == 0 but we can't test floating pt
        # equality
        if (cur_covar < 0.00001) {
            ret <- run_coverage_selected(nreps, cur_n, true_betas,
                                         select_method,
                                         direction = cur_select,
                                         split_prop = cur_split_prop)
        } else {
            ret <- run_coverage_selected(nreps, cur_n, true_betas,
                                         select_method,
                                         direction = cur_select,
                                         split_prop = cur_split_prop,
                                         gen_dist = gen_xs_corr,
                                         covar = cur_covar)
        }
        res_df <- data.frame(trial = rep(cur_trial, p),
                             select_method = rep(cur_select, p),
                             n = rep(cur_n, p),
                             covar = rep(cur_covar, p),
                             coef_name = names(ret$coverages),
                             coverage = ret$coverages,
                             select_prob = ret$select_probs)
    }
    final_df <- do.call(rbind, all_data)
    return(final_df)
}

# "Easy case"
true_betas1 <- c(1, 2, 0, 0.1)
ntrials <- 10
nreps <- 1000
ns <- c(50, 500)
select_method <- c("leaps", "forward", "backward", "both")
covar <- c(0, 0.5)
split_prop <- c(-1, 0.5, 0.8)

sim_conditions1 <- expand.grid(trial_index = 1:ntrials,
                              select_method = select_method,
                              n = ns,
                              covar = covar,
                              split_prop = split_prop,
                              stringsAsFactors = FALSE)

cl <- makeCluster(4)
registerDoParallel(cl)

start_time <- Sys.time()
beta1_df <- run_sim_parallel(true_betas1, sim_conditions1, nreps)
write.csv(beta1_df, "generated-data/beta1_df.csv")
print(Sys.time() - start_time)

# Harder case
true_betas2 <- rep(0, 20)
true_betas2[1] <- 1
true_betas2[2] <- 2
true_betas[20] <- 0.1

select_method2 <- c("forward", "backward", "both")

sim_conditions2 <- expand.grid(trial_index = 1:ntrials,
                              select_method = select_method2,
                              n = ns,
                              covar = covar,
                              split_prop = split_prop,
                              stringsAsFactors = FALSE)

start_time <- Sys.time()
beta2_df <- run_sim_parallel(true_betas2, sim_conditions2, nreps)
write.csv(beta2_df, "generated-data/beta2_df.csv")
print(Sys.time() - start_time)

stopCluster(cl)

library("methods1proj")
library("dplyr")
library("ggplot2")
library("cowplot")

set.seed(1262)

# helper function to run the simulation and generate plots
run_plot_sim <- function(nreps, n, true_betas,
                    select_method = c("step", "leaps"),
                    direction = "both", split_prop = -1,
                    gen_dist = gen_xs_default, error_dist = rnorm,
                    plot_coefs = FALSE, ...) {
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

    plot_res <- plot_sim(sim_res, true_betas)
    plots_lst <- plot_res$plots
    coverages <- plot_res$coverages

    if (plot_coefs) {
        plots_lst$labels <- LETTERS[1:4]
        print(do.call(plot_grid, plots_lst))
    }

    coverages_df <- data.frame(coef_name = names(coverages),
                               coverage = coverages)
    p_cover <- ggplot(data = coverages_df, aes(x = coef_name,
                                               y = coverage)) +
        geom_point() + ylim(min(0.949, min(coverages)), 1) +
        geom_hline(yintercept = 0.95, color = "red") +
        ggtitle("Coverage Probability") + xlab("Coefficient") +
        ylab("Coverage Probability")
    #print(p_cover)

    selected_df <- sim_res %>% group_by(coef_name) %>%
            summarize(pct_selected = n() / nreps)
    p_select <- ggplot(data = selected_df, aes(x = coef_name,
                                        y = pct_selected)) +
        geom_point() + ggtitle("Selection Probability") +
        xlab("Coefficient") + ylab("Selection Probability")
    #print(p_select)
}



# "Easy case"
true_betas1 <- c(1, 2, 0, 0.1)
nreps <- 1000
n1 <- 50
#par(mfrow = c(2, 2))
sim_types <- c("leaps", "forward", "backward", "both")

# TESTING, REMOVE LATER
for (m in sim_types) {
    select_method = "leaps"
    if (m != "leaps") {
        select_method = "step"
    }
    run_plot_sim(nreps, n1, true_betas1, plot_coefs = FALSE,
                 select_method = select_method, direction = m)
}

run_plot_sim(nreps, n1, true_betas1, plot_coefs = FALSE,
             select_method = select_method, direction = m,
             split_prop =  0.5)

# try with correlation
for (m in sim_types) {
    select_method = "leaps"
    if (m != "leaps") {
        select_method = "step"
    }
    run_plot_sim(nreps, n1, true_betas1, plot_coefs = FALSE,
                 select_method = select_method, direction = m,
                 gen_dist = gen_xs_corr, covar = 0.95)

}

# Bigger n
par(mfrow = c(2, 2))
n2 <- 500
for (m in sim_types) {
    select_method = "leaps"
    if (m != "leaps") {
        select_method = "step"
    }
    run_plot_sim(nreps, n2, true_betas1, plot_coefs = FALSE,
                 select_method = select_method, direction = m)

}

# Harder case
true_betas2 <- rep(0, 20)
true_betas2[1] <- 1
true_betas2[2] <- 2
step_dir <- c("forward", "backward", "both")
for (s in step_dir) {
    run_plot_sim(nreps, n1, true_betas2, plot_coefs = FALSE,
                 select_method = "step", direction = s)
}

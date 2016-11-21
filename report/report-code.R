library("methods1proj")
library("dplyr")
library("ggplot2")
library("cowplot")

set.seed(1262)

# helper function to run the simulation and generate plots
run_plot_sim <- function(nreps, n, true_betas,
                    select_method = c("step", "leaps"),
                    direction = "both",
                    gen_dist = gen_xs_default, error_dist = rnorm,
                    plot_coefs = FALSE, ...) {
    p <- length(true_betas)
    start <- Sys.time()
    sim_res <- run_sim(nreps, n, true_betas,
                       select_method = select_method)
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

    index <- seq(0, p - 1, 1)
    plot(index, coverages, ylim = c(min(coverages), 1))
    abline(h = 0.95, col = "red")
}

# "Easy case"
true_betas1 <- c(1, 2, 0, 0.1)
nreps <- 1000
n1 <- 50
par(mfrow = c(2, 2))
run_plot_sim(nreps, n1, true_betas1, plot_coefs = FALSE,
             select_method = "leaps")

run_plot_sim(nreps, n1, true_betas1, plot_coefs = FALSE,
             select_method = "step", direction = "forward")

run_plot_sim(nreps, n1, true_betas1, plot_coefs = FALSE,
             select_method = "step", direction = "backward")

run_plot_sim(nreps, n1, true_betas1, plot_coefs = FALSE,
             select_method = "step", direction = "both")

# Bigger n
par(mfrow = c(2, 2))
n2 <- 500
run_plot_sim(nreps, n2, true_betas1, plot_coefs = FALSE,
             select_method = "leaps")

run_plot_sim(nreps, n2, true_betas1, plot_coefs = FALSE,
             select_method = "step", direction = "forward")

run_plot_sim(nreps, n2, true_betas1, plot_coefs = FALSE,
             select_method = "step", direction = "backward")

run_plot_sim(nreps, n2, true_betas1, plot_coefs = FALSE,
             select_method = "step", direction = "both")

# Harder case
true_betas2 <- rep(0, 10)
true_betas2[1] <- 1
true_betas2[2] <- 2

run_plot_sim(nreps, n1, true_betas2, plot_coefs = FALSE,
             select_method = "step", direction = "forward")

run_plot_sim(nreps, n1, true_betas2, plot_coefs = FALSE,
             select_method = "step", direction = "backward")

run_plot_sim(nreps, n1, true_betas2, plot_coefs = FALSE,
             select_method = "step", direction = "both")

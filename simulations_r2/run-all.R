# =============================================================================
# run-all.R -- master driver for the round-2 simulation study.
# Usage:
#   Rscript run-all.R quick     (smoke test, ~2-5 minutes)
#   Rscript run-all.R full      (paper-quality runs)
# Outputs land in simulations_r2/outputs/.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else "quick"
cores <- max(1, parallel::detectCores() - 2)
t0 <- Sys.time()

if (mode == "full") {
  CFG_A <- list(B = 1000, ns = c(1000, 2000, 4000),
                alphas = c(0.10, 0.15, 0.20, 0.30, 0.40, 0.50),
                delta = 0.1, seed = 20260609, cores = cores)
  CFG_B <- list(B = 1000, n = 2000, alphas = seq(0.10, 0.90, by = 0.05),
                alphas_w = c(0.10, 0.20, 0.30, 0.40, 0.50), N_test = 50000,
                pi_hi = 0.20, delta = 0.1, seed = 20260610, cores = cores)
  CFG_C <- list(B = 300, T = 25, n_train = 2000, n_eval = 53,
                alpha_grid = seq(0.10, 0.95, by = 0.05), delta = 0.1,
                tau = 0.03, seed = 20260611, cores = cores)
  CFG_D <- list(B_C = 150, B_A = 300, T = 25, n_train = 2000, n_eval = 53,
                alpha_grid = seq(0.10, 0.95, by = 0.05),
                alphas_A = c(0.10, 0.20, 0.30, 0.50), delta = 0.1, tau = 0.03,
                seed = 20260612, cores = cores)
} else {
  CFG_A <- list(B = 200, ns = c(2000),
                alphas = c(0.10, 0.15, 0.20, 0.30, 0.40, 0.50),
                delta = 0.1, seed = 20260609, cores = cores)
  CFG_B <- list(B = 100, n = 2000, alphas = seq(0.10, 0.90, by = 0.10),
                alphas_w = c(0.10, 0.30, 0.50), N_test = 30000,
                pi_hi = 0.20, delta = 0.1, seed = 20260610, cores = cores)
  CFG_C <- list(B = 20, T = 25, n_train = 2000, n_eval = 53,
                alpha_grid = seq(0.10, 0.95, by = 0.05), delta = 0.1,
                tau = 0.03, seed = 20260611, cores = cores)
  CFG_D <- list(B_C = 10, B_A = 60, T = 25, n_train = 2000, n_eval = 53,
                alpha_grid = seq(0.10, 0.95, by = 0.05),
                alphas_A = c(0.10, 0.30, 0.50), delta = 0.1, tau = 0.03,
                seed = 20260612, cores = cores)
}

cat("Mode:", mode, "| cores:", cores, "\n\n")
source("sim-panelA.R"); cat("\n")
source("sim-panelB.R"); cat("\n")
source("sim-panelC.R"); cat("\n")
source("sim-panelD.R"); cat("\n")

writeLines(capture.output(sessionInfo()), "outputs/sessionInfo.txt")
cat(sprintf("All panels finished in %.1f minutes.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

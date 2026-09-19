# =============================================================================
# Panel D: robustness.
#   D1: time-varying interest rate r_t in the insurance model
#       (a) analyst uses the correct r_t path; (b) analyst wrongly assumes
#       the constant benchmark r = 0.05. Answers Reviewer 3, Q1.
#   D2: generation-side jitter (+/- 15%) of the insurance calibration
#       (g, lambda, gamma, delta.op); analyst recovers with the benchmark.
#   D3: heteroscedastic (QDA-type) DGP with a misspecified logistic learner;
#       the umbrella's missed-crisis control must still hold (Panel A rerun).
# =============================================================================
source("sim-helpers.R")

if (!exists("CFG_D")) {
  CFG_D <- list(B_C = 30, B_A = 100, T = 25, n_train = 2000, n_eval = 53,
                alpha_grid = seq(0.10, 0.95, by = 0.05),
                alphas_A = c(0.10, 0.20, 0.30, 0.50), delta = 0.1, tau = 0.03,
                seed = 20260612, cores = max(1, parallel::detectCores() - 2))
}

pars <- sim_params()
th <- benchmark_theta()
th$pi.bar <- pars$pi0
r_path <- seq(0.03, 0.07, length.out = CFG_D$T)

scenarios <- list(
  baseline   = list(r_gen = NULL,  r_rec = NULL,  jitter = 0),
  D1a_rt_known   = list(r_gen = r_path, r_rec = r_path, jitter = 0),
  D1b_rt_ignored = list(r_gen = r_path, r_rec = NULL,   jitter = 0),
  D2_calib_jitter = list(r_gen = NULL, r_rec = NULL,    jitter = 0.15)
)

summD <- do.call(rbind, lapply(names(scenarios), function(sc) {
  s <- scenarios[[sc]]
  res <- panelC_run(CFG_D$B_C, CFG_D$T, CFG_D$n_train, CFG_D$n_eval,
                    CFG_D$alpha_grid, CFG_D$delta, CFG_D$tau, pars, th,
                    seed = CFG_D$seed, cores = CFG_D$cores,
                    r_path_gen = s$r_gen, r_path_rec = s$r_rec,
                    th_gen_jitter = s$jitter)
  sm <- panelC_summarize(res)
  data.frame(scenario = sc, mean_abs_bias = mean(abs(sm$bias)),
             mean_rmse = mean(sm$rmse))
}))

# D3: umbrella control under a misspecified base learner
parsQ <- sim_params(sigma1_scale = 1.5)
resQ <- panelA_run(CFG_D$B_A, CFG_D$n_train, CFG_D$alphas_A, CFG_D$delta,
                   parsQ, seed = CFG_D$seed + 99, cores = CFG_D$cores)
sQ <- panelA_summarize(resQ, CFG_D$delta)
sQ <- sQ[sQ$method == "NP umbrella", c("alpha", "R1", "viol")]
sQ$scenario <- "D3_qda_misspecified"

write.csv(summD, "outputs/panelD_recovery_robustness.csv", row.names = FALSE)
write.csv(sQ, "outputs/panelD_qda_control.csv", row.names = FALSE)

cat("Panel D done.\nRecovery robustness (Panel C metric under perturbations):\n")
print(summD, row.names = FALSE, digits = 3)
cat("\nD3: umbrella violation rates under the misspecified (QDA) DGP:\n")
print(sQ, row.names = FALSE, digits = 3)

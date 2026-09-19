# =============================================================================
# Panel C: recovery of a known, time-varying risk tolerance alpha*_t.
# The policymaker generates observed reserves through the two-stage pipeline
# at the true alpha*_t; an analyst with an independent training draw recovers
# alpha via the Section 5.2 calibration (squared reserve error). Validates
# "revealed risk tolerance" as an estimator: bias, RMSE, turning points.
# =============================================================================
source("sim-helpers.R")

if (!exists("CFG_C")) {
  CFG_C <- list(B = 50, T = 25, n_train = 2000, n_eval = 53,
                alpha_grid = seq(0.10, 0.95, by = 0.05), delta = 0.1,
                tau = 0.03, seed = 20260611,
                cores = max(1, parallel::detectCores() - 2))
}

pars <- sim_params()
th <- benchmark_theta()
th$pi.bar <- pars$pi0

resC <- panelC_run(CFG_C$B, CFG_C$T, CFG_C$n_train, CFG_C$n_eval,
                   CFG_C$alpha_grid, CFG_C$delta, CFG_C$tau, pars, th,
                   seed = CFG_C$seed, cores = CFG_C$cores)
saveRDS(resC, "outputs/panelC_raw.rds")
summC <- panelC_summarize(resC)
write.csv(summC, "outputs/panelC_summary.csv", row.names = FALSE)

pdf("outputs/panelC_alpha_recovery.pdf", width = 8, height = 5)
par(mar = c(4.2, 4.4, 2.5, 1))
plot(summC$t, summC$alpha_true, type = "l", lwd = 2, ylim = c(0, 1),
     xlab = "period t", ylab = expression(alpha),
     main = "True vs recovered risk tolerance")
lines(summC$t, summC$alpha_hat_mean, lty = 2)
polygon(c(summC$t, rev(summC$t)),
        c(summC$alpha_hat_mean - summC$alpha_hat_sd,
          rev(summC$alpha_hat_mean + summC$alpha_hat_sd)),
        border = NA, col = rgb(0, 0, 0, 0.12))
legend("bottomleft", c("true alpha*_t", "recovered (MC mean)", "MC +/- 1 sd"),
       lwd = c(2, 1, 8), lty = c(1, 2, 1),
       col = c("black", "black", rgb(0, 0, 0, 0.12)), bty = "n", cex = 0.85)
dev.off()

cat("Panel C done. Overall recovery metrics:\n")
cat(sprintf("  mean |bias| = %.3f,  mean RMSE = %.3f\n",
            mean(abs(summC$bias)), mean(summC$rmse)))
print(summC[, c("t", "alpha_true", "alpha_hat_mean", "bias", "rmse")],
      row.names = FALSE, digits = 3)

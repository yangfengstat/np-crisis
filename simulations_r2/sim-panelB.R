# =============================================================================
# Panel B: welfare analysis. Sub-panels:
#   B1 (exact, no MC): Proposition 1 asymmetry. For the benchmark insurance
#       model, the welfare cost of a downward deviation of the estimated
#       crisis probability exceeds that of an upward deviation of the same
#       size by a factor of about 2.3 to 5.9 across the relevant risk range.
#       This is the welfare basis for prioritizing missed crises.
#   B1b (MC diagnostic): the Eq. (12) weights omega0/omega1 under a
#       consistent estimator are dominated by the class priors, so the
#       ratio omega1/omega0 is NOT generically above one. Kept as an honest
#       diagnostic; the paper's asymmetry argument should rest on B1 and on
#       the state-contingent decision costs, not on Eq. (12) generically.
#   B2: expected welfare cost of the aggregated two-stage pipeline along the
#       NP path under low (pi0) and elevated (pi_hi) true risk; the
#       welfare-optimal alpha falls as true risk rises, while the symmetric
#       rule is stuck at one dominated point.
#   B3: analytic NP <-> cost-sensitive crossover at alpha_bar = Phi(-Delta/2)
#       (generalized Proposition 2).
# =============================================================================
source("sim-helpers.R")

if (!exists("CFG_B")) {
  CFG_B <- list(B = 150, n = 2000, alphas = seq(0.10, 0.90, by = 0.10),
                alphas_w = c(0.10, 0.20, 0.30, 0.40, 0.50),
                N_test = 50000, pi_hi = 0.20, delta = 0.1, seed = 20260610,
                cores = max(1, parallel::detectCores() - 2))
}

pars <- sim_params()
th <- benchmark_theta()
th$pi.bar <- pars$pi0   # internal consistency: insurers see the true base rate

## ---- B1: Proposition 1 asymmetry (exact) -----------------------------------
pi_grid <- seq(0.03, 0.30, by = 0.005)
rel_devs <- c(0.25, 0.50, 0.75)
asym <- do.call(rbind, lapply(rel_devs, function(rd) {
  dn <- vapply(pi_grid, function(p) welfare_cost(p, p * (1 - rd), th), 0)
  up <- vapply(pi_grid, function(p) welfare_cost(p, p * (1 + rd), th), 0)
  data.frame(pi = pi_grid, rel_dev = rd, wc_down = dn, wc_up = up,
             ratio = dn / up)
}))
write.csv(asym, "outputs/panelB1_prop1_asymmetry.csv", row.names = FALSE)

## ---- B1b: Eq. (12) weights diagnostic (MC) ---------------------------------
resW <- panelB_weights(CFG_B$B, CFG_B$n, CFG_B$N_test, CFG_B$alphas_w,
                       pars, th, seed = CFG_B$seed + 1, cores = CFG_B$cores)
saveRDS(resW, "outputs/panelB1b_weights_raw.rds")
sW <- aggregate(cbind(omega0, omega1, ratio) ~ alpha, data = resW, FUN = mean)
write.csv(sW, "outputs/panelB1b_eq12_weights.csv", row.names = FALSE)

## ---- B2: pipeline welfare cost along the NP path ---------------------------
resB <- panelB_run(CFG_B$B, CFG_B$n, CFG_B$alphas, CFG_B$delta, pars, th,
                   seed = CFG_B$seed, cores = CFG_B$cores)
saveRDS(resB, "outputs/panelB2_raw.rds")
wc <- aggregate(cbind(share_np, share_sym) ~ alpha, data = resB, FUN = mean)
wc$wc_np_lo <- vapply(wc$share_np, function(s) welfare_cost(pars$pi0, s, th), 0)
wc$wc_np_hi <- vapply(wc$share_np, function(s) welfare_cost(CFG_B$pi_hi, s, th), 0)
wc_sym_lo <- welfare_cost(pars$pi0, wc$share_sym[1], th)
wc_sym_hi <- welfare_cost(CFG_B$pi_hi, wc$share_sym[1], th)
a_star_lo <- wc$alpha[which.min(wc$wc_np_lo)]
a_star_hi <- wc$alpha[which.min(wc$wc_np_hi)]
write.csv(wc, "outputs/panelB2_welfare.csv", row.names = FALSE)

## ---- B3: analytic crossover -------------------------------------------------
cross <- data.frame(alpha = seq(0.02, 0.98, by = 0.01))
cross$logD <- pars$Delta^2 / 2 + qnorm(cross$alpha) * pars$Delta
write.csv(cross, "outputs/panelB3_crossover.csv", row.names = FALSE)

## ---- figure -----------------------------------------------------------------
pdf("outputs/panelB_welfare.pdf", width = 12, height = 4.2)
par(mfrow = c(1, 3), mar = c(4.2, 4.4, 2.8, 1))
plot(NA, xlim = range(pi_grid), ylim = c(1, max(asym$ratio) * 1.05),
     xlab = expression("true probability " * pi),
     ylab = "welfare cost ratio: down / up",
     main = "B1: Prop. 1 asymmetry")
abline(h = 1, col = "grey60")
for (i in seq_along(rel_devs)) {
  ai <- asym[asym$rel_dev == rel_devs[i], ]
  lines(ai$pi, ai$ratio, lty = i)
}
legend("topright", sprintf("deviation = %d%% of pi", round(100 * rel_devs)),
       lty = seq_along(rel_devs), bty = "n", cex = 0.85)
yl <- range(c(wc$wc_np_lo, wc$wc_np_hi, wc_sym_lo, wc_sym_hi))
plot(wc$alpha, wc$wc_np_lo, type = "b", pch = 19, ylim = yl,
     xlab = expression(alpha), ylab = "expected welfare cost",
     main = "B2: pipeline welfare cost")
lines(wc$alpha, wc$wc_np_hi, type = "b", pch = 17, lty = 2)
abline(h = wc_sym_lo, lty = 3)
abline(h = wc_sym_hi, lty = 4)
legend("topleft", c(sprintf("NP, true risk = %.3f (opt %.1f)", pars$pi0, a_star_lo),
                    sprintf("NP, true risk = %.2f (opt %.1f)", CFG_B$pi_hi, a_star_hi),
                    "symmetric, low risk", "symmetric, high risk"),
       lty = c(1, 2, 3, 4), pch = c(19, 17, NA, NA), bty = "n", cex = 0.8)
plot(cross$alpha, cross$logD, type = "l", xlab = expression(alpha),
     ylab = expression(log * " " * D[alpha]),
     main = "B3: equal-weight crossover")
abline(h = 0, col = "grey60")
abline(v = pars$alpha_bar, lty = 3)
legend("topleft", legend = sprintf("alpha_bar = %.3f", pars$alpha_bar),
       lty = 3, bty = "n")
dev.off()

cat("Panel B done.\n")
cat(sprintf("B1: Prop-1 asymmetry ratio (down/up) ranges %.2f to %.2f over pi in [%.2f, %.2f]\n",
            min(asym$ratio), max(asym$ratio), min(pi_grid), max(pi_grid)))
cat("B1b (diagnostic): Eq-12 weight ratio omega1/omega0 by alpha:\n")
print(sW, row.names = FALSE, digits = 3)
cat(sprintf("B2: welfare-optimal alpha: %.2f at true risk %.3f; %.2f at true risk %.2f\n",
            a_star_lo, pars$pi0, a_star_hi, CFG_B$pi_hi))
cat(sprintf("    symmetric-rule welfare cost is %.1fx the NP optimum (low risk), %.1fx (high risk)\n",
            wc_sym_lo / min(wc$wc_np_lo), wc_sym_hi / min(wc$wc_np_hi)))
cat(sprintf("B3: analytic alpha_bar = %.3f\n", pars$alpha_bar))

# =============================================================================
# sim-helpers.R: engine for the revision simulation study (Panels A-D)
#
# NEW file for the second revision. Does not modify any existing analysis code.
# Design document: simulation-study-spec.md in the Overleaf project folder.
#
# Provenance of shared formulas:
#   - The NP umbrella threshold rule mirrors npc_aux.core() in
#     ../np-crisis-helper-funcs.R (Tong, Feng and Li, 2018, JMLR "nproc").
#   - rho_star(), p_of(), eta_of(), f_of(), U_real() mirror rho.star() in
#     ../np-crisis-helper-funcs.R and Section S.2 of the Supplementary Material.
#
# Dependencies: base R + stats + parallel only.
# =============================================================================

# ---------------------------------------------------------------------------
# 1. Data-generating process (homoscedastic Gaussian; QDA variant via scale)
# ---------------------------------------------------------------------------
# Class 0 (normal): X ~ N(0, I_d)
# Class 1 (SSGI):   X ~ N(mu1, S1), S1 = I except first s coords scaled
# Signal strength is set by the target Bayes AUC: Delta = sqrt(2)*qnorm(auc),
# so that the equal-weight threshold alpha_bar = pnorm(-Delta/2) is known.
sim_params <- function(auc = 0.70, d = 30, s = 5, pi0 = 0.042,
                       sigma1_scale = 1) {
  Delta <- sqrt(2) * qnorm(auc)
  mu1 <- c(rep(Delta / sqrt(s), s), rep(0, d - s))
  list(auc = auc, d = d, s = s, pi0 = pi0, Delta = Delta, mu1 = mu1,
       sigma1_scale = sigma1_scale, alpha_bar = pnorm(-Delta / 2))
}

gen_data <- function(n, pars) {
  y <- rbinom(n, 1, pars$pi0)
  X <- matrix(rnorm(n * pars$d), n, pars$d)
  i1 <- which(y == 1)
  if (length(i1) > 0) {
    if (pars$sigma1_scale != 1) {
      X[i1, seq_len(pars$s)] <- X[i1, seq_len(pars$s)] * pars$sigma1_scale
    }
    X[i1, ] <- sweep(X[i1, , drop = FALSE], 2, pars$mu1, `+`)
  }
  list(x = X, y = y)
}

# ---------------------------------------------------------------------------
# 2. Closed-form population errors for a linear rule: flag crisis iff w'x > c
# ---------------------------------------------------------------------------
pop_errors <- function(w, cthr, pars) {
  s0 <- sqrt(sum(w^2))                       # sd of w'X under class 0
  v1 <- w^2
  v1[seq_len(pars$s)] <- v1[seq_len(pars$s)] * pars$sigma1_scale^2
  s1 <- sqrt(sum(v1))                        # sd of w'X under class 1
  m1 <- sum(w * pars$mu1)
  c(R0 = 1 - pnorm(cthr / s0),               # false-alarm rate  P(flag | Y=0)
    R1 = pnorm((cthr - m1) / s1))            # missed-crisis rate P(no flag | Y=1)
}

# Population share of observations flagged as crises (the representative-
# country probability input of Section 5.2, in population form)
pop_share <- function(pe, pars, eps = 1e-3) {
  ph <- pars$pi0 * (1 - pe[["R1"]]) + (1 - pars$pi0) * pe[["R0"]]
  min(max(ph, eps), 1 - eps)
}

# ---------------------------------------------------------------------------
# 3. Scoring function (logistic regression; score = linear predictor)
# ---------------------------------------------------------------------------
fit_score <- function(x, y) {
  df <- data.frame(y = y, x)
  fit <- suppressWarnings(glm(y ~ ., data = df, family = binomial()))
  co <- coef(fit)
  co[is.na(co)] <- 0
  list(w = as.numeric(co[-1]), b = as.numeric(co[1]))
}
lp <- function(model, x) drop(x %*% model$w)   # intercept irrelevant: absorbed in threshold

# ---------------------------------------------------------------------------
# 4. Threshold rules on a common scoring function
# ---------------------------------------------------------------------------
# (a) NP umbrella order-statistic rule (mirrors npc_aux.core):
#     scores_out = scores of LEFT-OUT crisis observations; flag iff score > c.
#     Violation guarantee: P( R1 > alpha ) <= delta.
#     Feasible iff m >= log(delta)/log(1-alpha).
np_threshold <- function(scores_out, alpha, delta = 0.1) {
  m <- length(scores_out)
  ks <- seq_len(m)
  feas <- pbinom(ks - 1, m, alpha) <= delta
  if (!any(feas)) return(NA_real_)
  sort(scores_out)[max(ks[feas])]
}

# (b) Naive plug-in rule: empirical alpha-quantile of the same left-out crisis
#     scores. Targets R1 = alpha on average but has NO high-probability
#     guarantee (violation rate is about 1/2). This is the "tune to hit the
#     empirical constraint" practice shown to fail in Tong et al. (2018).
thr_plugin <- function(scores_out, alpha) {
  as.numeric(quantile(scores_out, alpha, type = 1))
}

# (c) Symmetric ("literature paradigm") rule: threshold minimizing the
#     empirical sum of the two error rates R0 + R1 on the training scores.
thr_symmetric <- function(scores0, scores1) {
  cand <- quantile(c(scores0, scores1), seq(0, 1, length.out = 512),
                   names = FALSE, type = 1)
  cand <- unique(cand)
  err <- vapply(cand, function(cc) mean(scores0 > cc) + mean(scores1 <= cc),
                numeric(1))
  cand[which.min(err)]
}

# One umbrella fit: split the crisis class (whose error is controlled),
# train the scoring function on all non-crises + half the crises,
# keep the other half of the crises out for thresholding.
np_umbrella <- function(x, y, split.ratio = 0.5) {
  i1 <- which(y == 1)
  n1 <- length(i1)
  i1.tr <- sample(i1, round(n1 * split.ratio))
  i1.out <- setdiff(i1, i1.tr)
  tr <- c(which(y == 0), i1.tr)
  mod <- fit_score(x[tr, , drop = FALSE], y[tr])
  list(model      = mod,
       scores_out = lp(mod, x[i1.out, , drop = FALSE]),
       scores1    = lp(mod, x[i1, , drop = FALSE]),
       scores0    = lp(mod, x[y == 0, , drop = FALSE]))
}

# ---------------------------------------------------------------------------
# 5. Insurance model and welfare (mirrors rho.star() and Supp. S.2)
# ---------------------------------------------------------------------------
# Benchmark calibration midpoints (Supplementary Material, Table S.1);
# pi.bar is reset to the DGP's pi0 at runtime for internal consistency.
benchmark_theta <- function() {
  list(sigma = 2, r = 0.05, g = 0.034, lambda = 0.069, gamma = 0.065,
       pi.bar = 0.042, delta.op = 0.011)
}

# p from (delta.op, pi.bar), exactly as in rho.star() in np-crisis-helper-funcs.R
p_of <- function(th) 1 - th$delta.op / ((1 - th$pi.bar) * (th$pi.bar + th$delta.op))

eta_of <- function(pi, th) {
  ((pi * (1 - th$pi.bar)) / (th$pi.bar * (1 - pi)) * p_of(th))^(1 / th$sigma)
}

rho_star <- function(pi, th, r = th$r) {
  eta <- eta_of(pi, th)
  (1 - (r - th$g) / (1 + th$g) * th$lambda -
     (1 / eta) * (1 - th$gamma - (1 + r) / (1 + th$g) * th$lambda)) /
    (1 - (1 - 1 / eta) * th$delta.op)
}

f_of <- function(pi, th, r = th$r) {
  p <- p_of(th)
  Theta <- p * (1 - th$pi.bar) * (1 - (r - th$g) / (1 + th$g) * th$lambda) +
    th$pi.bar * (1 - th$gamma - (1 + r) / (1 + th$g) * th$lambda)
  Theta / (p * (1 - th$pi.bar) + eta_of(pi, th) * th$pi.bar)
}

U_real <- function(pi_true, pi_hat, th, r = th$r) {
  s <- th$sigma
  f <- f_of(pi_hat, th, r)
  eta <- eta_of(pi_hat, th)
  (pi_true * (eta * f)^(1 - s) + (1 - pi_true) * f^(1 - s)) / (1 - s)
}

welfare_cost <- function(pi_true, pi_hat, th, r = th$r) {
  U_real(pi_true, pi_true, th, r) - U_real(pi_true, pi_hat, th, r)
}

# True posterior probability P(Y=1 | x) for the homoscedastic Gaussian DGP
true_posterior <- function(x, pars) {
  stopifnot(pars$sigma1_scale == 1)
  logLR <- drop(x %*% pars$mu1) - sum(pars$mu1^2) / 2
  plogis(qlogis(pars$pi0) + logLR)
}

# ---------------------------------------------------------------------------
# 6. Panel engines
# ---------------------------------------------------------------------------

# ---- Panel A: population-level missed-crisis control ----------------------
panelA_run <- function(B, n, alphas, delta, pars, seed, cores = 1) {
  one <- function(b) {
    set.seed(seed + b)
    dat <- gen_data(n, pars)
    if (sum(dat$y) < 8) return(NULL)               # degenerate draw guard
    um <- np_umbrella(dat$x, dat$y)
    c_sym <- thr_symmetric(um$scores0, um$scores1)
    pe_sym <- pop_errors(um$model$w, c_sym, pars)
    rows <- lapply(alphas, function(a) {
      c_np <- np_threshold(um$scores_out, a, delta)
      c_pi <- thr_plugin(um$scores_out, a)
      pe_np <- if (is.na(c_np)) c(R0 = NA, R1 = NA)
               else pop_errors(um$model$w, c_np, pars)
      pe_pi <- pop_errors(um$model$w, c_pi, pars)
      data.frame(rep = b, alpha = a,
                 method = c("NP umbrella", "plug-in quantile", "symmetric"),
                 R0 = c(pe_np[["R0"]], pe_pi[["R0"]], pe_sym[["R0"]]),
                 R1 = c(pe_np[["R1"]], pe_pi[["R1"]], pe_sym[["R1"]]))
    })
    do.call(rbind, rows)
  }
  res <- parallel::mclapply(seq_len(B), one, mc.cores = cores)
  do.call(rbind, res[!vapply(res, is.null, logical(1))])
}

panelA_summarize <- function(res, delta) {
  agg <- aggregate(cbind(R0, R1) ~ alpha + method, data = res, FUN = mean,
                   na.rm = TRUE, na.action = NULL)
  viol <- aggregate(cbind(viol = R1 > alpha) ~ alpha + method, data = res,
                    FUN = mean, na.rm = TRUE, na.action = NULL)
  infeas <- aggregate(cbind(infeasible = is.na(R1)) ~ alpha + method,
                      data = res, FUN = mean, na.action = NULL)
  out <- Reduce(function(a, b) merge(a, b, by = c("alpha", "method")),
                list(agg, viol, infeas))
  out$delta_nominal <- delta
  out[order(out$method, out$alpha), ]
}

# ---- Panel B: welfare cost along the NP path vs the symmetric rule --------
panelB_run <- function(B, n, alphas, delta, pars, th, seed, cores = 1) {
  one <- function(b) {
    set.seed(seed + b)
    dat <- gen_data(n, pars)
    if (sum(dat$y) < 8) return(NULL)
    um <- np_umbrella(dat$x, dat$y)
    c_sym <- thr_symmetric(um$scores0, um$scores1)
    sh_sym <- pop_share(pop_errors(um$model$w, c_sym, pars), pars)
    rows <- lapply(alphas, function(a) {
      c_np <- np_threshold(um$scores_out, a, delta)
      if (is.na(c_np)) return(NULL)
      sh_np <- pop_share(pop_errors(um$model$w, c_np, pars), pars)
      data.frame(rep = b, alpha = a,
                 share_np = sh_np, wc_np = welfare_cost(pars$pi0, sh_np, th),
                 share_sym = sh_sym, wc_sym = welfare_cost(pars$pi0, sh_sym, th))
    })
    do.call(rbind, rows)
  }
  res <- parallel::mclapply(seq_len(B), one, mc.cores = cores)
  do.call(rbind, res[!vapply(res, is.null, logical(1))])
}

# ---- Panel B1: welfare-implied error weights (main-text Eq. 12) -----------
# Monte Carlo integration of omega0 = E[dw | false alarm] P(Y=0) and
# omega1 = E[dw | missed crisis] P(Y=1), where dw = welfare_cost(pi(x), pihat(x))
# with pi(x) the true posterior and pihat(x) a finite-sample logistic estimate.
# Demonstrates the welfare-cost asymmetry omega1 > omega0 in a continuous-X
# model (generalizing the two-point illustration of Supp. S.3).
panelB_weights <- function(B, n, N_test, alphas, pars, th, seed, cores = 1) {
  one <- function(b) {
    set.seed(seed + b)
    dat <- gen_data(n, pars)
    if (sum(dat$y) < 8) return(NULL)
    mod <- fit_score(dat$x, dat$y)
    te <- gen_data(N_test, pars)
    sc <- lp(mod, te$x)
    ph <- pmin(pmax(plogis(mod$b + sc), 1e-3), 1 - 1e-3)
    pt <- pmin(pmax(true_posterior(te$x, pars), 1e-6), 1 - 1e-6)
    dw <- welfare_cost(pt, ph, th)
    s1 <- lp(mod, dat$x[dat$y == 1, , drop = FALSE])
    rows <- lapply(alphas, function(a) {
      c_a <- as.numeric(quantile(s1, a, type = 1))
      flag <- sc > c_a
      fa <- flag & te$y == 0
      mc <- (!flag) & te$y == 1
      if (sum(fa) < 20 || sum(mc) < 20) return(NULL)
      w0 <- mean(dw[fa]) * (1 - pars$pi0)
      w1 <- mean(dw[mc]) * pars$pi0
      data.frame(rep = b, alpha = a, omega0 = w0, omega1 = w1,
                 ratio = w1 / w0)
    })
    do.call(rbind, rows)
  }
  res <- parallel::mclapply(seq_len(B), one, mc.cores = cores)
  do.call(rbind, res[!vapply(res, is.null, logical(1))])
}

# ---- Panel C: recovery of a known, time-varying risk tolerance ------------
# Generation: in period t, the policymaker fits the NP model at the TRUE
#   alpha*_t on their own training draw, flags a fresh cross-section of
#   n_eval "countries", converts the flagged share into rho*, and observed
#   reserves = rho* x (1 + noise), noise ~ N(0, tau^2).
# Recovery: an analyst with an INDEPENDENT training draw of the same DGP
#   re-runs the pipeline over a grid of alpha and picks the alpha minimizing
#   the squared reserve error (mirrors reserves_match.mse.alpha()).
alpha_true_path <- function(T) {
  knots_t <- c(1, 6, 10, 18, 25)
  knots_a <- c(0.50, 0.90, 0.45, 0.90, 0.40)
  approx(knots_t, knots_a, xout = seq(1, 25, length.out = T))$y
}

empirical_share <- function(model, cthr, x_eval, eps = 1e-3) {
  sh <- mean(lp(model, x_eval) > cthr)
  min(max(sh, eps), 1 - eps)
}

panelC_run <- function(B, T, n_train, n_eval, alpha_grid, delta, tau,
                       pars, th, seed, cores = 1,
                       r_path_gen = NULL, r_path_rec = NULL,
                       th_gen_jitter = 0) {
  a_true <- alpha_true_path(T)
  one <- function(b) {
    set.seed(seed + b)
    out <- vector("list", T)
    for (t in seq_len(T)) {
      r_gen <- if (is.null(r_path_gen)) th$r else r_path_gen[t]
      r_rec <- if (is.null(r_path_rec)) th$r else r_path_rec[t]
      th_gen <- th
      if (th_gen_jitter > 0) {
        for (nm in c("g", "lambda", "gamma", "delta.op")) {
          th_gen[[nm]] <- th[[nm]] * (1 + runif(1, -th_gen_jitter, th_gen_jitter))
        }
      }
      x_eval <- gen_data(n_eval, pars)$x        # current-year cross-section
      # policymaker side
      dp <- gen_data(n_train, pars)
      if (sum(dp$y) < 8) next
      up <- np_umbrella(dp$x, dp$y)
      c_p <- np_threshold(up$scores_out, a_true[t], delta)
      if (is.na(c_p)) next
      sh_p <- empirical_share(up$model, c_p, x_eval)
      rho_obs <- rho_star(sh_p, th_gen, r = r_gen) * (1 + rnorm(1, 0, tau))
      # analyst side (independent draw, same DGP and information structure)
      da <- gen_data(n_train, pars)
      if (sum(da$y) < 8) next
      ua <- np_umbrella(da$x, da$y)
      rho_a <- vapply(alpha_grid, function(a) {
        c_a <- np_threshold(ua$scores_out, a, delta)
        if (is.na(c_a)) return(NA_real_)
        rho_star(empirical_share(ua$model, c_a, x_eval), th, r = r_rec)
      }, numeric(1))
      a_hat <- alpha_grid[which.min((rho_a - rho_obs)^2)]
      out[[t]] <- data.frame(rep = b, t = t, alpha_true = a_true[t],
                             alpha_hat = a_hat, rho_obs = rho_obs,
                             share_pol = sh_p)
    }
    do.call(rbind, out)
  }
  res <- parallel::mclapply(seq_len(B), one, mc.cores = cores)
  do.call(rbind, res[!vapply(res, is.null, logical(1))])
}

panelC_summarize <- function(res) {
  agg <- aggregate(alpha_hat ~ t + alpha_true, data = res,
                   FUN = function(z) c(mean = mean(z), sd = sd(z)))
  out <- data.frame(t = agg$t, alpha_true = agg$alpha_true,
                    alpha_hat_mean = agg$alpha_hat[, "mean"],
                    alpha_hat_sd = agg$alpha_hat[, "sd"])
  out$bias <- out$alpha_hat_mean - out$alpha_true
  rmse <- aggregate((alpha_hat - alpha_true)^2 ~ t, data = res, FUN = mean)
  out$rmse <- sqrt(rmse[[2]][match(out$t, rmse$t)])
  out[order(out$t), ]
}

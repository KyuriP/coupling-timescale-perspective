# =============================================================================
# Revision analyses: "Rethinking Group Differences in Psychopathology Networks:
# A Slow-Fast Perspective on Context and Symptom Activation"
#
# Addresses the Nature Mental Health reviews (decision 6 Oct 2026):
#   A  Revised Slow Risk Load (SRL): adversity coding of work variables, removal
#      of healthcare cultural distance, winsorizing, minimum-completion rule  (R1.2)
#   B  Non-overlapping SRL (drops indicators that may reflect depression)  (R1.3)
#   C  Exact moderated Ising model with continuous SRL moderating activation
#      (h) and interactions (J); LRTs, AIC/BIC, common-slope test, implied
#      probabilities                                         (R1.1, R1.4, R1.5, R1.7, R2.3, R2.4)
#   C0 Median-split multigroup models (replicates Table 1 as a check)
#   D  Domain-specific moderators                                       (R1.8)
#   E  Alternative symptom cut-off (item >= 2); ordinal mgm cross-check (R1.4)
#   F  Ethnicity: moderation of J, SRL x ethnicity on h, and multigroup
#      models within each ethnic group                                  (R1.10)
#   G  Simulated parameters for the redrawn Fig. 1                      (R2.2)
#
# The core estimator fits, by exact maximum likelihood over all 2^9 = 512
# symptom states,
#     h_i(n)  = Wh[n, ] %*% Bh[, i]          (activation, i = 1..9)
#     J_ij(n) = WJ[n, ] %*% BJ[, ij]         (interactions, 36 pairs)
# A binary group dummy in Wh and WJ reproduces the fully group-specific
# multigroup model exactly; a continuous SRL gives a moderated Ising model
# (Haslbeck, 2022, without regularization).
#
# Edit the CONFIG block, then run:   Rscript revision_analyses.R
# Run with  RUN_SELFTEST <- TRUE  first to check the estimator on simulated data
# (gradient and Hessian vs finite differences, parameter recovery, and that the
# dummy-coded model equals separate group fits).
#
# Run time: ~5 min for n = 6,000 with OpenBLAS; expect ~20-30 min for the full
# HELIUS sample. R's default reference BLAS is ~5x slower: on macOS switch R to
# Accelerate/vecLib, on Linux install libopenblas. The mgm block can take longer.
# Only base R is required for everything except the mgm cross-check.
#
# Section C0 is a built-in check: with the original SRL (ORIG_SRL_VAR) and the
# original median split, it must reproduce manuscript Table 1. If it does not,
# the data preparation here differs from the original pipeline.
# =============================================================================

suppressPackageStartupMessages({
  library(stats)
  library(utils)
})

# ----------------------------------------------------------------------------
# CONFIG  -- every variable name marked TODO must be checked against the
#           HELIUS file. Indicator names are those listed in Supplementary
#           Table S5.
# ----------------------------------------------------------------------------
RUN_SELFTEST <- FALSE           # TRUE: run the simulation checks and stop
RUN_MGM      <- TRUE            # ordinal mgm cross-check (slow; needs package mgm)
DATA_FILE    <- "helius_baseline.csv"   # TODO path (csv; or replace read step)
OUT_DIR      <- "revision_output"
SEED         <- 20261020

# PHQ-9 items, scored 0-3, in this order                                   TODO
PHQ_ITEMS <- c(anh = "H1_PHQ9_1", dep = "H1_PHQ9_2", slp = "H1_PHQ9_3",
               ene = "H1_PHQ9_4", app = "H1_PHQ9_5", glt = "H1_PHQ9_6",
               con = "H1_PHQ9_7", mot = "H1_PHQ9_8", sui = "H1_PHQ9_9")

AGE_VAR       <- "H1_lft"         # TODO age in years
SEX_VAR       <- "H1_geslacht"    # TODO
FEMALE_CODE   <- 2                # TODO code for female
ETH_VAR       <- "H1_EtnTotaal"   # TODO ethnicity
DUTCH_CODE    <- 1                # TODO code for Dutch origin
ORIG_SRL_VAR  <- NA               # optional: name of the original SRL score,
                                  # to replicate Table 1 as a check

# Missing-value codes used in the original preparation
MISSING_CODES <- c(-9, -1)

# SRL indicators. direction = +1 if higher = more adverse, -1 if protective.
# recode = function applied to the raw vector before standardization.
adverse_worksit <- function(x) ifelse(is.na(x), NA, as.numeric(x %in% c(6, 7, 8)))
  # [DECIDE] unemployed (6), unable to work (7), social benefits (8) = 1;
  # paid work (1-4), retired (5), homemaker (9), student (10) = 0
adverse_labour  <- function(x) ifelse(is.na(x), NA, as.numeric(x %in% c(3, 4)))
  # [DECIDE] unemployed/welfare (3), incapacitated (4) = 1
stress_recode   <- function(x) ifelse(x %in% 1:4, x, NA)   # 5 = does not apply
identity_recode <- function(x) x

INDICATORS <- list(
  # domain "ses": socioeconomic resources
  list(name = "income_difficulty", var = "H1_InkHhMoeite",        domain = "ses",  dir = +1, recode = identity_recode, overlap = FALSE),
  list(name = "work_situation",    var = "H1_WerkSit",            domain = "ses",  dir = +1, recode = adverse_worksit, overlap = TRUE),
  list(name = "education",         var = "H1_Opleid",             domain = "ses",  dir = -1, recode = identity_recode, overlap = FALSE),
  list(name = "labour_part",       var = "H1_Arbeidsparticipatie", domain = "ses", dir = +1, recode = adverse_labour,  overlap = TRUE),
  list(name = "occupation_level",  var = "H1_BeroepsNiveau",      domain = "ses",  dir = -1, recode = identity_recode, overlap = FALSE),
  list(name = "health_literacy",   var = "H1_SBSQ_meanscore",     domain = "ses",  dir = -1, recode = identity_recode, overlap = FALSE),
  list(name = "social_support",    var = "H1_SSQT",               domain = "ses",  dir = -1, recode = identity_recode, overlap = TRUE),
  # domain "str": psychosocial and structural stress
  list(name = "discrimination",    var = "H1_Discr_sumscore",     domain = "str",  dir = +1, recode = identity_recode, overlap = FALSE),
  list(name = "language_diff",     var = "H1_MoeiteNLtaal",       domain = "str",  dir = +1, recode = identity_recode, overlap = FALSE, dutch_zero = TRUE),
  list(name = "work_stress",       var = "H1_StressWerk",         domain = "str",  dir = +1, recode = stress_recode,   overlap = TRUE),
  list(name = "home_stress",       var = "H1_StressHuish",        domain = "str",  dir = +1, recode = stress_recode,   overlap = TRUE),
  # domain "hlt": health and lifestyle
  list(name = "physical_func",     var = "H1_PCS12",              domain = "hlt",  dir = -1, recode = identity_recode, overlap = TRUE),
  list(name = "physical_activity", var = "H1_Squash_totmwk",      domain = "hlt",  dir = -1, recode = identity_recode, overlap = TRUE, winsor = TRUE),
  list(name = "smoking",           var = "H1_PackYears",          domain = "hlt",  dir = +1, recode = identity_recode, overlap = FALSE, winsor = TRUE),
  list(name = "alcohol",           var = "H1_AlcoholConsumption", domain = "hlt",  dir = +1, recode = identity_recode, overlap = FALSE)
)
# Healthcare cultural distance (H1_CultDistMeanScore6_corrected) is deliberately
# excluded: it applies only to experience of healthcare abroad, and Dutch-origin
# participants had been assigned 0, outside its 1-12 range.

MIN_DOMAIN_FRAC <- 0.5   # a domain score needs >= half of its indicators
MIN_DOMAINS     <- 2     # SRL needs >= 2 of 3 domain scores
WINSOR_P        <- c(0.01, 0.99)


# ============================================================================
# 1. Exact (moderated) Ising estimator
# ============================================================================

ising_space <- function(p) {
  S  <- as.matrix(expand.grid(rep(list(0:1), p)))
  pr <- t(utils::combn(p, 2))
  S2 <- S[, pr[, 1], drop = FALSE] * S[, pr[, 2], drop = FALSE]
  sets <- c(as.list(seq_len(p)), lapply(seq_len(nrow(pr)), function(r) pr[r, ]))
  K <- length(sets)
  ukey <- character(0); uidx <- matrix(0L, K, K); Tsub <- list()
  for (k in seq_len(K)) for (l in k:K) {
    u <- sort(unique(c(sets[[k]], sets[[l]]))); key <- paste(u, collapse = ",")
    j <- match(key, ukey)
    if (is.na(j)) {
      ukey <- c(ukey, key); j <- length(ukey)
      Tsub[[j]] <- apply(S[, u, drop = FALSE], 1, prod)
    }
    uidx[k, l] <- uidx[l, k] <- j
  }
  list(p = p, m = nrow(pr), pairs = pr, Tst = cbind(S, S2),   # 2^p x K
       Tsub = do.call(cbind, Tsub), uidx = uidx)
}

pair_products <- function(X, pr) X[, pr[, 1], drop = FALSE] * X[, pr[, 2], drop = FALSE]

row_logsumexp <- function(E) {
  mx <- apply(E, 1, max)
  mx + log(rowSums(exp(E - mx)))
}

# Linear map from free parameters theta to the full parameter vector
# phi = c(vec(Bh), vec(BJ)). Columns of Wh named in h_common share one
# coefficient across all symptoms (used for the common-slope test).
build_map <- function(h_names, J_names, p, m, h_common = character(0)) {
  qh <- length(h_names); qJ <- length(J_names)
  nphi <- qh * p + qJ * m
  lab_phi <- c(outer(h_names, seq_len(p), function(a, i) paste0("h", i, ":", a)),
               outer(J_names, seq_len(m), function(a, k) paste0("J", k, ":", a)))
  theta_id <- integer(nphi); theta_lab <- character(0)
  common_id <- setNames(integer(length(h_common)), h_common)
  for (i in seq_len(p)) for (a in seq_len(qh)) {
    idx <- a + (i - 1) * qh
    nm <- h_names[a]
    if (nm %in% h_common) {
      if (common_id[nm] == 0) {
        theta_lab <- c(theta_lab, paste0("h*:", nm)); common_id[nm] <- length(theta_lab)
      }
      theta_id[idx] <- common_id[nm]
    } else {
      theta_lab <- c(theta_lab, lab_phi[idx]); theta_id[idx] <- length(theta_lab)
    }
  }
  for (k in seq_len(m)) for (a in seq_len(qJ)) {
    idx <- qh * p + a + (k - 1) * qJ
    theta_lab <- c(theta_lab, lab_phi[idx]); theta_id[idx] <- length(theta_lab)
  }
  M <- matrix(0, nphi, length(theta_lab))
  M[cbind(seq_len(nphi), theta_id)] <- 1
  list(M = M, theta_lab = theta_lab, phi_lab = lab_phi, qh = qh, qJ = qJ)
}

# Core: log-likelihood, gradient, and (optionally) exact Hessian.
# For an exponential family with a linear predictor the observed and expected
# information coincide: H = - sum_n D_n' Cov_n(T) D_n.
ising_eval <- function(theta, X2, Wh, WJ, map, sp, w, hessian = FALSE) {
  p <- sp$p; m <- sp$m; K <- p + m
  phi <- drop(map$M %*% theta)
  Bh <- matrix(phi[seq_len(map$qh * p)], map$qh, p)
  BJ <- matrix(phi[map$qh * p + seq_len(map$qJ * m)], map$qJ, m)
  eta <- cbind(Wh %*% Bh, WJ %*% BJ)                  # n x K
  E <- eta %*% t(sp$Tst)                               # n x 2^p
  lz <- row_logsumexp(E)
  ll <- sum(w * (rowSums(eta * X2) - lz))
  P <- exp(E - lz)
  mu <- P %*% sp$Tst                                   # n x K
  R <- (X2 - mu) * w
  g_phi <- c(crossprod(Wh, R[, 1:p, drop = FALSE]),
             crossprod(WJ, R[, p + seq_len(m), drop = FALSE]))
  out <- list(ll = ll, grad = drop(crossprod(map$M, g_phi)))
  if (hessian) {
    D <- function(k) if (k <= p) Wh else WJ
    off <- function(k) if (k <= p) (k - 1) * map$qh else map$qh * p + (k - p - 1) * map$qJ
    q   <- function(k) if (k <= p) map$qh else map$qJ
    nphi <- nrow(map$M)
    H <- matrix(0, nphi, nphi)
    # E[T_k T_l] is the probability that the union of their symptom sets is
    # all present: only 255 distinct unions (size <= 4) for p = 9.
    EU <- P %*% sp$Tsub                                  # n x n_unions
    for (k in seq_len(K)) {
      Dk <- D(k) * w
      for (l in k:K) {
        ckl <- EU[, sp$uidx[k, l]] - mu[, k] * mu[, l]   # Cov_n(T_k, T_l)
        blk <- -crossprod(Dk * ckl, D(l))
        ik <- off(k) + seq_len(q(k)); il <- off(l) + seq_len(q(l))
        H[ik, il] <- blk
        if (l != k) H[il, ik] <- t(blk)
      }
    }
    out$H <- crossprod(map$M, H %*% map$M)
  }
  out
}

# Fit by Newton-Raphson with step halving (the log-likelihood is concave in
# theta). Returns estimates, SEs and fit statistics.
fit_ising <- function(X, Wh, WJ, h_common = character(0), w = NULL,
                      tol = 1e-8, maxit = 100, verbose = FALSE, label = "") {
  X <- as.matrix(X); p <- ncol(X)
  sp <- ising_space(p)
  X2 <- cbind(X, pair_products(X, sp$pairs))
  if (is.null(w)) w <- rep(1, nrow(X))
  map <- build_map(colnames(Wh), colnames(WJ), p, sp$m, h_common)
  theta <- numeric(ncol(map$M))
  # start: intercept of h at marginal logits
  pm <- pmin(pmax(colMeans(X), 1e-3), 1 - 1e-3)
  for (i in seq_len(p)) {
    j <- match(paste0("h", i, ":(Intercept)"), map$theta_lab)
    if (!is.na(j)) theta[j] <- qlogis(pm[i])
  }
  ev <- ising_eval(theta, X2, Wh, WJ, map, sp, w, hessian = TRUE)
  for (it in seq_len(maxit)) {
    step <- tryCatch(solve(-ev$H, ev$grad), error = function(e) MASS_ginv(-ev$H) %*% ev$grad)
    step <- drop(step); s <- 1
    repeat {
      new <- ising_eval(theta + s * step, X2, Wh, WJ, map, sp, w, hessian = FALSE)
      if (is.finite(new$ll) && new$ll >= ev$ll - 1e-10) break
      s <- s / 2; if (s < 1e-8) break
    }
    theta <- theta + s * step
    dec <- sum(step * ev$grad) * s
    ev <- ising_eval(theta, X2, Wh, WJ, map, sp, w, hessian = TRUE)
    if (verbose) cat(sprintf("  [%s] it %d  ll %.4f  |grad| %.2e\n", label, it, ev$ll, max(abs(ev$grad))))
    if (abs(dec) < tol && max(abs(ev$grad)) < 1e-4) break
  }
  V <- tryCatch(solve(-ev$H), error = function(e) MASS_ginv(-ev$H))
  se <- sqrt(pmax(diag(V), 0))
  n <- sum(w); k <- length(theta)
  structure(list(theta = setNames(theta, map$theta_lab), se = setNames(se, map$theta_lab),
                 vcov = V, ll = ev$ll, k = k, n = n,
                 AIC = -2 * ev$ll + 2 * k, BIC = -2 * ev$ll + log(n) * k,
                 max_grad = max(abs(ev$grad)), iter = it, label = label,
                 p = p, m = sp$m, pairs = sp$pairs, symptoms = colnames(X),
                 h_names = colnames(Wh), J_names = colnames(WJ), map = map),
            class = "ising_fit")
}

MASS_ginv <- function(A, tol = sqrt(.Machine$double.eps)) {
  s <- svd(A); d <- s$d; pos <- d > max(tol * d[1], 0)
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / d[pos])
}

lrt <- function(restricted, full, test = "") {
  stat <- 2 * (full$ll - restricted$ll); df <- full$k - restricted$k
  data.frame(test = test, restricted = restricted$label, full = full$label,
             chisq = stat, df = df, p = pchisq(stat, df, lower.tail = FALSE))
}

fit_table <- function(...) {
  fits <- list(...)
  do.call(rbind, lapply(fits, function(f)
    data.frame(model = f$label, k = f$k, LL = f$ll, AIC = f$AIC, BIC = f$BIC,
               max_grad = f$max_grad)))
}

# Coefficients of one covariate on activation, per symptom, with Wald CI,
# odds ratio and Holm-adjusted p.
h_coef <- function(fit, term) {
  lab <- paste0("h", seq_len(fit$p), ":", term)
  if (all(!lab %in% names(fit$theta))) lab <- rep(paste0("h*:", term), fit$p)
  est <- fit$theta[lab]; se <- fit$se[lab]
  z <- est / se; pv <- 2 * pnorm(-abs(z))
  data.frame(symptom = fit$symptoms, term = term, est = est, se = se,
             lo = est - 1.96 * se, hi = est + 1.96 * se,
             OR = exp(est), OR_lo = exp(est - 1.96 * se), OR_hi = exp(est + 1.96 * se),
             p = pv, p_holm = p.adjust(pv, "holm"), row.names = NULL)
}

J_coef <- function(fit, term) {
  lab <- paste0("J", seq_len(fit$m), ":", term)
  est <- fit$theta[lab]; se <- fit$se[lab]; pv <- 2 * pnorm(-abs(est / se))
  data.frame(edge = paste(fit$symptoms[fit$pairs[, 1]], fit$symptoms[fit$pairs[, 2]], sep = "--"),
             term = term, est = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se,
             p = pv, p_holm = p.adjust(pv, "holm"), row.names = NULL)
}

# Model-implied marginal symptom probabilities for a covariate profile.
implied_marginals <- function(fit, wh_row, wJ_row) {
  sp <- ising_space(fit$p)
  phi <- drop(fit$map$M %*% fit$theta)
  Bh <- matrix(phi[seq_len(fit$map$qh * fit$p)], fit$map$qh, fit$p)
  BJ <- matrix(phi[fit$map$qh * fit$p + seq_len(fit$map$qJ * fit$m)], fit$map$qJ, fit$m)
  eta <- c(wh_row %*% Bh, wJ_row %*% BJ)
  e <- drop(sp$Tst %*% eta); pr <- exp(e - max(e)); pr <- pr / sum(pr)
  setNames(drop(pr %*% sp$Tst[, 1:fit$p]), fit$symptoms)
}

# Convenience: design matrices from a data frame
design <- function(df, terms) {
  Wm <- matrix(1, nrow(df), 1, dimnames = list(NULL, "(Intercept)"))
  if (length(terms)) Wm <- cbind(Wm, as.matrix(df[, terms, drop = FALSE]))
  Wm
}

# The four-model family used throughout: moderator(s) on h and/or on J.
# cov_h: covariates always on h; cov_J: covariates always on J.
four_models <- function(X, df, mod, cov_h = character(0), cov_J = character(0), tag = "") {
  f <- function(h, J, lab) fit_ising(X, design(df, h), design(df, J), label = paste0(tag, lab))
  full  <- f(c(cov_h, mod), c(cov_J, mod), "both moderated")
  hOnly <- f(c(cov_h, mod), cov_J,          "activation moderated only")
  JOnly <- f(cov_h,          c(cov_J, mod), "interactions moderated only")
  none  <- f(cov_h,          cov_J,         "neither moderated")
  list(fits = list(full = full, hOnly = hOnly, JOnly = JOnly, none = none),
       table = fit_table(full, hOnly, JOnly, none),
       lrt = rbind(lrt(hOnly, full, "interaction moderation"),
                   lrt(JOnly, full, "activation moderation")))
}


# ============================================================================
# 2. Self-test on simulated data
# ============================================================================
selftest <- function() {
  set.seed(SEED)
  p <- 9; sp <- ising_space(p); m <- sp$m; n <- 6000
  z <- rnorm(n); g <- rbinom(n, 1, 0.5)
  Bh <- rbind(rnorm(p, -1.5, 0.5), runif(p, 0.1, 0.5), rnorm(p, 0, 0.2))   # int, z, g
  BJ <- rbind(rnorm(m, 0.4, 0.4), c(rep(0, m - 3), 0.5, -0.4, 0.3))         # int, z
  Wh <- cbind(1, z, g); WJ <- cbind(1, z)
  eta <- cbind(Wh %*% Bh, WJ %*% BJ)
  E <- eta %*% t(sp$Tst); P <- exp(E - row_logsumexp(E))
  st <- apply(P, 1, function(pr) sample.int(nrow(sp$Tst), 1, prob = pr))
  X <- sp$Tst[st, 1:p]; colnames(X) <- paste0("s", 1:p)
  df <- data.frame(z = z, g = g)

  cat("Self-test 1: gradient vs finite differences\n")
  Whd <- design(df, c("z", "g")); WJd <- design(df, "z")
  map <- build_map(colnames(Whd), colnames(WJd), p, m)
  X2 <- cbind(X, pair_products(X, sp$pairs)); th <- rnorm(ncol(map$M), 0, 0.1)
  ev <- ising_eval(th, X2, Whd, WJd, map, sp, rep(1, n), hessian = TRUE)
  idx <- sample(length(th), 5); fd <- sapply(idx, function(j) {
    e <- numeric(length(th)); e[j] <- 1e-5
    (ising_eval(th + e, X2, Whd, WJd, map, sp, rep(1, n))$ll -
       ising_eval(th - e, X2, Whd, WJd, map, sp, rep(1, n))$ll) / 2e-5 })
  print(round(cbind(analytic = ev$grad[idx], numeric = fd), 4))
  fdH <- (ising_eval(th + replace(numeric(length(th)), idx[1], 1e-5), X2, Whd, WJd, map, sp, rep(1, n))$grad -
          ising_eval(th - replace(numeric(length(th)), idx[1], 1e-5), X2, Whd, WJd, map, sp, rep(1, n))$grad) / 2e-5
  cat(sprintf("  max |Hessian column - finite diff| = %.2e\n", max(abs(ev$H[, idx[1]] - fdH))))

  cat("\nSelf-test 2: parameter recovery\n")
  fit <- fit_ising(X, Whd, WJd, label = "sim")
  hz <- h_coef(fit, "z")
  cat(sprintf("  h slope on z: corr(true, est) = %.3f, mean |error| = %.3f, coverage = %.2f\n",
              cor(Bh[2, ], hz$est), mean(abs(Bh[2, ] - hz$est)),
              mean(Bh[2, ] >= hz$lo & Bh[2, ] <= hz$hi)))
  jz <- J_coef(fit, "z")
  cat(sprintf("  J slope on z: coverage = %.2f; truly non-zero edges recovered (Holm p<.05): %d of 3\n",
              mean(BJ[2, ] >= jz$lo & BJ[2, ] <= jz$hi), sum(jz$p_holm[(m - 2):m] < .05)))

  cat("\nSelf-test 3: binary-dummy model equals two separate group fits\n")
  fg <- fit_ising(X, design(df, "g"), design(df, "g"), label = "dummy")
  f0 <- fit_ising(X[g == 0, ], design(df[g == 0, , drop = FALSE], character(0)),
                  design(df[g == 0, , drop = FALSE], character(0)), label = "g0")
  f1 <- fit_ising(X[g == 1, ], design(df[g == 1, , drop = FALSE], character(0)),
                  design(df[g == 1, , drop = FALSE], character(0)), label = "g1")
  cat(sprintf("  LL dummy = %.4f ; LL g0 + g1 = %.4f\n", fg$ll, f0$ll + f1$ll))

  cat("\nSelf-test 4: LRTs on the four-model family\n")
  fm <- four_models(X, df, "z", cov_h = "g", tag = "sim: ")
  print(fm$table, digits = 8); print(fm$lrt, digits = 4)
  cs <- fit_ising(X, design(df, c("g", "z")), design(df, "z"), h_common = "z", label = "common slope")
  print(lrt(cs, fm$fits$full), digits = 4)
  invisible(TRUE)
}

if (RUN_SELFTEST) { selftest(); quit(save = "no") }


# ============================================================================
# 3. Data and Slow Risk Load (Section A, B)
# ============================================================================
dir.create(OUT_DIR, showWarnings = FALSE)
set.seed(SEED)
options(width = 160)
log_file <- file.path(OUT_DIR, "results_summary.txt")
sink(log_file, split = TRUE)
cat("Revision analyses run", format(Sys.time()), "\n\n")

raw <- read.csv(DATA_FILE)     # TODO: replace with haven::read_sav() for SPSS
cat("Records in source dataset:", nrow(raw), "\n")

clean_codes <- function(x) { x[x %in% MISSING_CODES] <- NA; x }
dutch <- raw[[ETH_VAR]] == DUTCH_CODE

winsorize <- function(x, pr = WINSOR_P) {
  q <- quantile(x, pr, na.rm = TRUE); pmin(pmax(x, q[1]), q[2])
}

build_srl <- function(raw, indicators) {
  Z <- sapply(indicators, function(ind) {
    x <- ind$recode(clean_codes(raw[[ind$var]]))
    if (isTRUE(ind$dutch_zero)) x[dutch & is.na(x)] <- 0
    if (isTRUE(ind$winsor)) x <- winsorize(x)
    ind$dir * (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)   # standardized over source data
  })
  colnames(Z) <- sapply(indicators, `[[`, "name")
  dom <- sapply(indicators, `[[`, "domain")
  doms <- unique(dom)
  D <- sapply(doms, function(d) {
    Zd <- Z[, dom == d, drop = FALSE]
    nobs <- rowSums(!is.na(Zd))
    ifelse(nobs >= MIN_DOMAIN_FRAC * ncol(Zd), rowMeans(Zd, na.rm = TRUE), NA)
  })
  colnames(D) <- doms
  nd <- rowSums(!is.na(D))
  srl <- ifelse(nd >= MIN_DOMAINS, rowMeans(D, na.rm = TRUE), NA)
  miss <- data.frame(indicator = colnames(Z), domain = dom,
                     n_missing = colSums(is.na(Z)), pct_missing = 100 * colMeans(is.na(Z)))
  list(srl = srl, domains = D, Z = Z, n_ind = rowSums(!is.na(Z)), missing = miss)
}

srl_rev <- build_srl(raw, INDICATORS)                                       # Section A
srl_nov <- build_srl(raw, Filter(function(i) !i$overlap, INDICATORS))      # Section B

# Analytic sample
items <- sapply(PHQ_ITEMS, function(v) clean_codes(raw[[v]]))
colnames(items) <- names(PHQ_ITEMS)
d <- data.frame(
  age = raw[[AGE_VAR]], female = as.numeric(raw[[SEX_VAR]] == FEMALE_CODE),
  nonDutch = as.numeric(!dutch),
  srl = srl_rev$srl, srl_nov = srl_nov$srl,
  ses = srl_rev$domains[, "ses"], str = srl_rev$domains[, "str"], hlt = srl_rev$domains[, "hlt"]
)
if (!is.na(ORIG_SRL_VAR)) d$srl_orig <- raw[[ORIG_SRL_VAR]]

flow <- data.frame(step = character(), n = integer())
add_flow <- function(lab, keep) { flow[nrow(flow) + 1, ] <<- list(lab, sum(keep)); keep }
k0 <- add_flow("Source records", rep(TRUE, nrow(d)))
k1 <- add_flow("Revised SRL computable (completion rule)", k0 & !is.na(d$srl))
k2 <- add_flow("+ age, sex, ethnicity observed", k1 & complete.cases(d[, c("age", "female", "nonDutch")]))
k3 <- add_flow("+ all nine PHQ-9 items observed (network sample)", k2 & complete.cases(items))
cat("\n== Section A: exclusion flow ==\n"); print(flow, row.names = FALSE)
write.csv(flow, file.path(OUT_DIR, "A_exclusion_flow.csv"), row.names = FALSE)
write.csv(srl_rev$missing, file.path(OUT_DIR, "A_indicator_missingness_revised.csv"), row.names = FALSE)
cat("\nIndicators observed per participant (revised SRL):\n"); print(summary(srl_rev$n_ind[k3]))

net <- d[k3, ]; Xo <- items[k3, ]
X1 <- (Xo >= 1) * 1; X2c <- (Xo >= 2) * 1           # symptom codings
std <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
net$age_z <- std(net$age)
for (v in c("srl", "srl_nov", "ses", "str", "hlt")) net[[paste0(v, "_z")]] <- std(net[[v]])
COV_H <- c("age_z", "female", "nonDutch")
cat(sprintf("\nNetwork sample n = %d\n", nrow(net)))
cat(sprintf("Correlation revised SRL with non-overlapping SRL: %.3f\n",
            cor(net$srl, net$srl_nov, use = "complete.obs")))
if (!is.null(net$srl_orig))
  cat(sprintf("Correlation revised SRL with original SRL: %.3f\n", cor(net$srl, net$srl_orig, use = "complete.obs")))


# ============================================================================
# 4. Section C0: median-split multigroup models (replication of Table 1)
# ============================================================================
multigroup <- function(X, grp, tag) {
  df <- data.frame(g = grp)
  four_models(X, df, "g", tag = tag)
}
cat("\n== Section C0: median-split multigroup models ==\n")
if (!is.null(d$srl_orig)) {
  # Original sample: SRL available (median taken here, n = 23,689 in the
  # manuscript), networks on complete PHQ items (n = 23,283). Independent of
  # the revised-SRL completion rule. TODO: if the original analysis required a
  # supplied PHQ total, add that condition to s_desc.
  s_desc <- !is.na(d$srl_orig)
  s_net  <- s_desc & complete.cases(items)
  med_o  <- median(d$srl_orig[s_desc])
  cat(sprintf("Original-SRL samples: descriptive n = %d, network n = %d\n", sum(s_desc), sum(s_net)))
  mg_orig <- multigroup((items[s_net, ] >= 1) * 1, as.numeric(d$srl_orig[s_net] > med_o), "orig SRL split: ")
  cat("Original SRL (should reproduce Table 1: LL -93555.20 / -93614.78 / -93611.51 / -95229.03)\n")
  print(mg_orig$table, digits = 8); print(mg_orig$lrt, digits = 4)
}
mg_rev <- multigroup(X1, as.numeric(net$srl > median(d$srl[k1])), "revised SRL split: ")
print(mg_rev$table, digits = 8); print(mg_rev$lrt, digits = 4)
write.csv(mg_rev$table, file.path(OUT_DIR, "C0_multigroup_revisedSRL.csv"), row.names = FALSE)


# ============================================================================
# 5. Section C: exact moderated Ising model, continuous SRL
#    h ~ age + female + nonDutch + SRL ;  J ~ SRL (+ nonDutch, Section F)
# ============================================================================
run_moderated <- function(X, mod, tag) {
  fm <- four_models(X, net, mod, cov_h = COV_H, cov_J = "nonDutch", tag = tag)
  cs <- fit_ising(X, design(net, c(COV_H, mod)), design(net, c("nonDutch", mod)),
                  h_common = mod, label = paste0(tag, "common SRL slope on h"))
  fm$common_slope <- cs
  fm$lrt <- rbind(fm$lrt, lrt(cs, fm$fits$full, "equal SRL slopes on h (symptom specificity)"))
  fm$h <- h_coef(fm$fits$full, mod)
  fm$J <- J_coef(fm$fits$full, mod)
  fm
}
report <- function(fm, file_tag) {
  print(fm$table, digits = 8)
  print(fm$lrt[, c("test", "chisq", "df", "p")], digits = 4)
  cat("SRL effects on activation (per 1 SD; full model):\n"); print(fm$h, digits = 3)
  cat(sprintf("Edges with Holm p < .05 for SRL moderation: %d of %d\n", sum(fm$J$p_holm < .05), nrow(fm$J)))
  write.csv(fm$table, file.path(OUT_DIR, paste0(file_tag, "_fit.csv")), row.names = FALSE)
  write.csv(fm$lrt,   file.path(OUT_DIR, paste0(file_tag, "_lrt.csv")), row.names = FALSE)
  write.csv(fm$h,     file.path(OUT_DIR, paste0(file_tag, "_h_slopes.csv")), row.names = FALSE)
  write.csv(fm$J,     file.path(OUT_DIR, paste0(file_tag, "_J_slopes.csv")), row.names = FALSE)
}

cat("\n== Section C: moderated Ising model, revised SRL, item >= 1 ==\n")
modC <- run_moderated(X1, "srl_z", "C: "); report(modC, "C_moderated_srl")

# Implied marginal probabilities at -1 / +1 SD SRL (reference: mean age, male,
# non-Dutch; change the profile if preferred)                          (R2.4)
prof <- function(z) list(h = c(1, 0, 0, 1, z), J = c(1, 1, z))
imp <- rbind(low = implied_marginals(modC$fits$full, prof(-1)$h, prof(-1)$J),
             high = implied_marginals(modC$fits$full, prof(+1)$h, prof(+1)$J))
cat("\nModel-implied symptom probabilities at SRL -1 SD / +1 SD:\n"); print(round(imp, 3))
write.csv(t(imp), file.path(OUT_DIR, "C_implied_marginals.csv"))

# Section B: same model with non-overlapping SRL                          (R1.3)
cat("\n== Section B: moderated Ising model, non-overlapping SRL ==\n")
modB <- run_moderated(X1, "srl_nov_z", "B: "); report(modB, "B_moderated_srl_nonoverlap")


# ============================================================================
# 6. Section D: domain-specific moderators                                (R1.8)
# ============================================================================
cat("\n== Section D: domains ==\n")
DOMS <- c("ses_z", "str_z", "hlt_z")
dom_sep <- lapply(DOMS, function(v) run_moderated(X1, v, paste0("D[", v, "]: ")))
names(dom_sep) <- DOMS
for (v in DOMS) { cat("\n-- domain", v, "alone --\n"); report(dom_sep[[v]], paste0("D_", v)) }

cat("\n-- domains jointly --\n")
dj_full <- fit_ising(X1, design(net, c(COV_H, DOMS)), design(net, c("nonDutch", DOMS)), label = "D joint: full")
dj_drop <- lapply(DOMS, function(v)
  fit_ising(X1, design(net, c(COV_H, setdiff(DOMS, v))), design(net, c("nonDutch", DOMS)),
            label = paste0("D joint: drop ", v, " from h")))
dj_lrt <- do.call(rbind, Map(function(f, v) lrt(f, dj_full, paste("activation moderation by", v, "| other domains")),
                             dj_drop, DOMS))
print(dj_lrt, digits = 4)
dj_h <- do.call(rbind, lapply(DOMS, function(v) h_coef(dj_full, v)))
print(dj_h, digits = 3)
write.csv(dj_lrt, file.path(OUT_DIR, "D_joint_lrt.csv"), row.names = FALSE)
write.csv(dj_h, file.path(OUT_DIR, "D_joint_h_slopes.csv"), row.names = FALSE)


# ============================================================================
# 7. Section E: alternative cut-off and ordinal check                     (R1.4)
# ============================================================================
cat("\n== Section E: symptom presence = item >= 2 ==\n")
cat("Prevalence at >= 2:\n"); print(round(colMeans(X2c), 3))
modE <- run_moderated(X2c, "srl_z", "E: "); report(modE, "E_moderated_cut2")
mgE <- multigroup(X2c, as.numeric(net$srl > median(d$srl[k1])), "E split: ")
print(mgE$table, digits = 8); print(mgE$lrt, digits = 4)

if (RUN_MGM && requireNamespace("mgm", quietly = TRUE)) {
  cat("\n-- Ordinal cross-check: moderated mixed graphical model (mgm) --\n")
  # PHQ items as 4-category variables, SRL as Gaussian moderator.
  dm <- cbind(Xo, net$srl_z, net$age_z, net$female, net$nonDutch)
  type <- c(rep("c", 9), "g", "g", "c", "c")
  level <- c(rep(4, 9), 1, 1, 2, 2)
  ok <- apply(Xo, 2, function(x) length(unique(x)) == 4)
  if (!all(ok)) warning("Some PHQ items do not use all 4 categories; mgm needs level = observed categories")
  mg <- tryCatch(mgm::mgm(data = dm, type = type, level = level, moderators = 10,
                          lambdaSel = "EBIC", lambdaGam = 0.25, ruleReg = "AND",
                          scale = FALSE, pbar = FALSE, signInfo = FALSE),
                 error = function(e) { message("mgm failed: ", conditionMessage(e)); NULL })
  if (!is.null(mg)) {
    # NOTE: not tested in the preparation environment (mgm unavailable there);
    # check the structure of mg$interactions against ?mgm for your version.
    pw <- mg$interactions$indicator[[1]]; th <- mg$interactions$indicator[[2]]
    lab <- c(names(PHQ_ITEMS), "SRL", "age", "female", "nonDutch")
    srl_rows <- pw[pw[, 1] == 10 | pw[, 2] == 10, , drop = FALSE]
    cat("Symptoms with a selected SRL main effect:\n")
    print(lab[ifelse(srl_rows[, 1] == 10, srl_rows[, 2], srl_rows[, 1])])
    cat("3-way terms (SRL moderating an edge):\n")
    if (!is.null(th) && nrow(th)) print(t(apply(th, 1, function(r) lab[r]))) else cat("  none selected\n")
    saveRDS(mg, file.path(OUT_DIR, "E_mgm_ordinal.rds"))
  }
} else cat("\n(mgm not installed; skipping ordinal cross-check)\n")


# ============================================================================
# 8. Section F: ethnicity                                                (R1.10)
# ============================================================================
cat("\n== Section F: ethnicity ==\n")
# F1: does ethnicity moderate interactions? (full model vs. no nonDutch on J)
f_noEthJ <- fit_ising(X1, design(net, c(COV_H, "srl_z")), design(net, "srl_z"),
                      label = "F: no ethnicity on J")
print(lrt(f_noEthJ, modC$fits$full, "ethnicity moderation of interactions"), digits = 4)
# F2: does the SRL effect on activation differ by ethnicity?
net$srl_x_nonDutch <- net$srl_z * net$nonDutch
f_int <- fit_ising(X1, design(net, c(COV_H, "srl_z", "srl_x_nonDutch")),
                   design(net, c("nonDutch", "srl_z")), label = "F: SRL x ethnicity on h")
print(lrt(modC$fits$full, f_int, "SRL x ethnicity on activation"), digits = 4)
print(h_coef(f_int, "srl_x_nonDutch"), digits = 3)
# F3: multigroup median-split models within each ethnic group
med_all <- median(d$srl[k1])
F3 <- lapply(c(Dutch = 0, nonDutch = 1), function(e) {
  s <- net$nonDutch == e
  multigroup(X1[s, ], as.numeric(net$srl[s] > med_all), paste0("F3 ", ifelse(e, "non-Dutch", "Dutch"), ": "))
})
for (e in names(F3)) { cat("\n--", e, "--\n"); print(F3[[e]]$table, digits = 8); print(F3[[e]]$lrt, digits = 4) }
# F4: moderated model within each ethnic group
F4 <- lapply(c(Dutch = 0, nonDutch = 1), function(e) {
  s <- net$nonDutch == e
  four_models(X1[s, ], net[s, ], "srl_z", cov_h = c("age_z", "female"), tag = paste0("F4 ", e, ": "))
})
for (e in names(F4)) { cat("\n-- moderated within", e, "--\n"); print(F4[[e]]$lrt, digits = 4)
  print(h_coef(F4[[e]]$fits$full, "srl_z"), digits = 3) }


# ============================================================================
# 9. Section G: parameters for the redrawn Fig. 1                         (R2.2)
#    Scenario B: groups differ in one interaction (guilt--concentration).
#    Scenario C: shared interactions; group-2 activation chosen so that its
#    marginal (phi) association network matches scenario B's group 2.
# ============================================================================
cat("\n== Section G: Fig. 1 parameters ==\n")
ising_dist <- function(h, J, sp) {
  e <- drop(sp$Tst %*% c(h, J[sp$pairs])); pr <- exp(e - max(e)); pr / sum(pr)
}
phi_matrix <- function(pr, sp) {
  S <- sp$Tst[, 1:sp$p]; mu <- drop(pr %*% S); E2 <- t(S) %*% (S * pr)
  C <- E2 - outer(mu, mu); C / sqrt(outer(diag(C), diag(C)))
}
spG <- ising_space(9); syms <- names(PHQ_ITEMS)
fitJ <- modC$fits$full      # use the fitted interaction intercepts as a realistic base
phiF <- drop(fitJ$map$M %*% fitJ$theta)
J0 <- matrix(0, 9, 9)
J0[spG$pairs] <- matrix(phiF[fitJ$map$qh * 9 + seq_len(fitJ$map$qJ * 36)], fitJ$map$qJ)[1, ]
J0 <- J0 + t(J0)
h0 <- matrix(phiF[seq_len(fitJ$map$qh * 9)], fitJ$map$qh)[1, ]
gi <- which(syms == "glt"); ci <- which(syms == "con")
JB <- J0; JB[gi, ci] <- JB[ci, gi] <- J0[gi, ci] + 1.0          # extra process in group 2
target <- phi_matrix(ising_dist(h0, JB, spG), spG)
obj <- function(h) { R <- phi_matrix(ising_dist(h, J0, spG), spG) - target; sum(R[upper.tri(R)]^2) }
hC <- optim(h0, obj, method = "BFGS")$par
achieved <- phi_matrix(ising_dist(hC, J0, spG), spG)
cat(sprintf("Fig. 1 match: RMS difference in phi between scenarios B and C (group 2) = %.4f\n",
            sqrt(mean((achieved - target)[upper.tri(target)]^2))))
base <- phi_matrix(ising_dist(h0, J0, spG), spG)
dimnames(base) <- dimnames(target) <- dimnames(achieved) <- list(syms, syms)
write.csv(base, file.path(OUT_DIR, "G_fig1_group1_phi.csv"))
write.csv(target, file.path(OUT_DIR, "G_fig1_scenarioB_group2_phi.csv"))
write.csv(achieved, file.path(OUT_DIR, "G_fig1_scenarioC_group2_phi.csv"))
write.csv(data.frame(symptom = syms, h_group1 = h0, h_scenarioC_group2 = hC),
          file.path(OUT_DIR, "G_fig1_activation.csv"), row.names = FALSE)

saveRDS(list(modC = modC, modB = modB, dom_sep = dom_sep, dj_full = dj_full, modE = modE,
             mg_rev = mg_rev, F3 = F3, F4 = F4), file.path(OUT_DIR, "all_fits.rds"))
cat("\nDone. Outputs in", OUT_DIR, "\n")
sink()

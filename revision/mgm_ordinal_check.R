# =============================================================================
# revision/mgm_ordinal_check.R
#
# Cross-check of the moderated Ising results with mgm. The PHQ-9 items keep
# their 0-3 response categories (as categorical nodes, so the order is not
# used) and continuous Slow Risk Load is the moderator.
#
# Run from the repository root in RStudio:   source("revision/mgm_ordinal_check.R")
#
# Input : data/network_sample_ordinal.csv
#         (same analytic sample as the exact Ising analyses; written by
#         python3 revision/revision_analyses.py --sections export)
# Output: results/revision_analyses/mgm_*.csv and mgm_*.rds
#
# QUICK_TEST <- TRUE: 3,000-person subsample, written to quick_test/
# =============================================================================

QUICK_TEST <- FALSE

if (!requireNamespace("mgm", quietly = TRUE)) stop("Install package mgm before running this script.", call. = FALSE)
library(mgm)

out_dir <- "results/revision_analyses"
if (QUICK_TEST) out_dir <- file.path(out_dir, "quick_test")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
d <- read.csv(file.path("data", "network_sample_ordinal.csv"))
if (QUICK_TEST) { set.seed(1); d <- d[sample.int(nrow(d), min(nrow(d), 3000L)), ] }
symptoms <- c("anh", "dep", "slp", "ene", "app", "glt", "con", "mot", "sui")

run_mgm <- function(moderator, item_type = c("ordinal", "binary")) {
  item_type <- match.arg(item_type)
  X <- d[, symptoms]
  if (item_type == "binary") X[] <- lapply(X, function(x) as.integer(x >= 1))
  dat <- as.matrix(cbind(X, d[, c(moderator, "age_z", "female", "nonDutch")]))
  labels <- c(symptoms, "SRL", "age", "female", "nonDutch")
  type  <- c(rep("c", 9), "g", "g", "c", "c")
  level <- c(sapply(X, function(x) length(unique(x))), 1, 1, 2, 2)
  tag <- paste0("mgm_", item_type, "_", moderator)
  cat("\n==", tag, "==\n")
  set.seed(20261020)
  t0 <- Sys.time()
  if (anyNA(dat) || any(!is.finite(dat))) stop("Non-finite values in ", tag, call. = FALSE)
  fit <- mgm(data = dat, type = type, level = level, moderators = 10,
             lambdaSel = "EBIC", lambdaGam = 0.25, ruleReg = "AND",
             scale = FALSE, pbar = TRUE)
  cat("time:", format(Sys.time() - t0), "\n")
  saveRDS(fit, file.path(out_dir, paste0(tag, ".rds")))

  # pairwise terms: symptom-symptom edges and SRL main effects
  # (weightsAgg are aggregated absolute weights, so no sign)
  pw  <- fit$interactions$indicator[[1]]
  pww <- unlist(fit$interactions$weightsAgg[[1]])
  pair_tab <- data.frame(v1 = labels[pw[, 1]], v2 = labels[pw[, 2]], weight = pww)
  # 3-way terms: SRL moderating a pairwise relation
  tw  <- fit$interactions$indicator[[2]]
  tww <- unlist(fit$interactions$weightsAgg[[2]])
  three_tab <- if (length(tw) && nrow(tw)) data.frame(
    v1 = labels[tw[, 1]], v2 = labels[tw[, 2]], v3 = labels[tw[, 3]], weight = tww) else
    data.frame(v1 = character(), v2 = character(), v3 = character(), weight = numeric())

  srl_main <- subset(pair_tab, (v1 == "SRL" & v2 %in% symptoms) | (v2 == "SRL" & v1 %in% symptoms))
  sym_edges <- subset(pair_tab, v1 %in% symptoms & v2 %in% symptoms)
  srl_mod <- subset(three_tab, (v1 == "SRL" | v2 == "SRL" | v3 == "SRL"))
  srl_mod_sym <- srl_mod[apply(srl_mod[, 1:3], 1, function(r) sum(r %in% symptoms) == 2), ]

  cat(sprintf("Symptom-symptom edges selected: %d of 36\n", nrow(sym_edges)))
  cat(sprintf("Symptoms with an SRL main effect (activation analogue): %d of 9\n", nrow(srl_main)))
  print(srl_main, row.names = FALSE)
  cat(sprintf("Symptom-symptom edges moderated by SRL: %d of 36\n", nrow(srl_mod_sym)))
  print(srl_mod_sym, row.names = FALSE)

  write.csv(pair_tab,  file.path(out_dir, paste0(tag, "_pairwise.csv")), row.names = FALSE)
  write.csv(three_tab, file.path(out_dir, paste0(tag, "_threeway.csv")), row.names = FALSE)
  invisible(fit)
}

sink(file.path(out_dir, "mgm_summary.txt"), split = TRUE)
cat("mgm version:", as.character(packageVersion("mgm")), "\n")
cat("Rows:", nrow(d), " Quick test:", QUICK_TEST, "\n")
fits <- list(
  ordinal_srl     = run_mgm("srl_z",     "ordinal"),
  ordinal_nonover = run_mgm("srl_nov_z", "ordinal"),   # SRL without overlapping indicators
  binary_srl      = run_mgm("srl_z",     "binary")     # regularized version of the binary model
)
sink()
cat("\nDone. Summary in", file.path(out_dir, "mgm_summary.txt"), "\n")

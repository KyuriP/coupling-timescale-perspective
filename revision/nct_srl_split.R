# ==============================================================================
# revision/nct_srl_split.R
#
# Network Comparison Test for the Low and High Slow Risk groups (median split
# of Slow Risk Load), on the same analytic sample as the moderated Ising
# analyses. Settings as in R/05_nct_final_10000.R.
#
# Input : data/network_sample_ordinal.csv
#         (python3 revision/revision_analyses.py --sections export)
# Output: results/revision_analyses/nct/
#
# Run from the repository root:
#
#   source("revision/nct_srl_split.R")
#
# ==============================================================================

SEED <- 20260729L
NCT_ITERATIONS <- 10000L
NCT_GAMMA <- 0.25

if (!requireNamespace("NetworkComparisonTest", quietly = TRUE)) {
  stop("Package 'NetworkComparisonTest' is required.", call. = FALSE)
}

symptoms <- c("anh", "dep", "slp", "ene", "app", "glt", "con", "mot", "sui")

d <- read.csv(file.path("data", "network_sample_ordinal.csv"))
if (!"srl_high" %in% names(d)) {
  stop("Re-run: python3 revision/revision_analyses.py --sections export", call. = FALSE)
}

binary <- as.data.frame(lapply(d[, symptoms], function(x) as.integer(x >= 1)))
group_1 <- binary[d$srl_high == 0, symptoms]
group_2 <- binary[d$srl_high == 1, symptoms]
cat("Low SRL n =", nrow(group_1), " High SRL n =", nrow(group_2), "\n")

output_dir <- file.path("results", "revision_analyses", "nct")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(SEED)
nct_result <- NetworkComparisonTest::NCT(
  data1 = group_1,
  data2 = group_2,
  gamma = NCT_GAMMA,
  it = NCT_ITERATIONS,
  binary.data = TRUE,
  paired = FALSE,
  weighted = TRUE,
  AND = TRUE,
  abs = TRUE,
  test.edges = TRUE,
  edges = "all",
  progressbar = TRUE,
  make.positive.definite = TRUE,
  p.adjust.methods = "holm",
  test.centrality = FALSE,
  verbose = TRUE
)

strengths <- as.numeric(nct_result$glstrinv.sep)

write.csv(
  data.frame(
    n_group_1 = nrow(group_1),
    n_group_2 = nrow(group_2),
    network_structure_statistic_M = as.numeric(nct_result$nwinv.real)[1],
    network_structure_p = as.numeric(nct_result$nwinv.pval)[1],
    global_strength_group_1 = strengths[1],
    global_strength_group_2 = strengths[2],
    global_strength_difference_S = as.numeric(nct_result$glstrinv.real)[1],
    global_strength_p = as.numeric(nct_result$glstrinv.pval)[1],
    permutations = NCT_ITERATIONS
  ),
  file.path(output_dir, "nct_summary.csv"),
  row.names = FALSE
)

if (!is.null(nct_result$einv.pvals)) {
  edge_p <- as.data.frame(nct_result$einv.pvals)
  write.csv(edge_p, file.path(output_dir, "nct_edge_tests.csv"), row.names = FALSE)
}

# Unadjusted edge p-values, using the same (b + 1)/(B + 1) correction as NCT.
saveRDS(nct_result, file.path(output_dir, "nct_result.rds"))
if (is.null(nct_result$einv.perm) || is.null(nct_result$einv.real) ||
    length(dim(nct_result$einv.perm)) != 3L) {
  stop("Expected edge statistics and a three-dimensional permutation array.", call. = FALSE)
}
{
  real <- nct_result$einv.real
  perm <- nct_result$einv.perm
  idx <- which(upper.tri(real), arr.ind = TRUE)
  edge_tab <- data.frame(
    var1 = symptoms[idx[, 1]],
    var2 = symptoms[idx[, 2]],
    difference = real[idx],
    p_unadjusted = vapply(seq_len(nrow(idx)), function(k) {
      (sum(abs(perm[idx[k, 1], idx[k, 2], ]) >= abs(real[idx[k, 1], idx[k, 2]])) + 1) /
        (dim(perm)[3] + 1)
    }, numeric(1))
  )
  edge_tab$p_holm <- p.adjust(edge_tab$p_unadjusted, "holm")
  write.csv(edge_tab, file.path(output_dir, "nct_edge_tests_unadjusted.csv"), row.names = FALSE)
}

write.csv(nct_result$nw1, file.path(output_dir, "nct_network_low.csv"), row.names = TRUE)
write.csv(nct_result$nw2, file.path(output_dir, "nct_network_high.csv"), row.names = TRUE)

message("NCT complete. Results in ", normalizePath(output_dir, mustWork = FALSE))

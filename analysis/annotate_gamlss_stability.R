#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  parsed <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    name <- sub("^--", "", key)
    if (i == length(args)) {
      stop("Missing value for argument: ", key, call. = FALSE)
    }
    i <- i + 1L
    parsed[[name]] <- args[[i]]
    i <- i + 1L
  }
  parsed
}

collapse_flags <- function(flag_names, flag_values) {
  active <- flag_names[vapply(flag_values, isTRUE, FUN.VALUE = logical(1))]
  if (length(active) == 0L) {
    "none"
  } else {
    paste(active, collapse = "; ")
  }
}

safe_quantile <- function(x, prob, fallback) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    fallback
  } else {
    stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE)
  }
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
run_dir <- args[["run-dir"]]
if (is.null(run_dir) || !nzchar(run_dir)) {
  stop("Usage: Rscript annotate_gamlss_stability.R --run-dir path/to/run_dir", call. = FALSE)
}

run_dir <- normalizePath(run_dir, mustWork = TRUE)
tables_dir <- file.path(run_dir, "tables")

results_file <- file.path(tables_dir, "per_gene_results.csv")
counts_file <- file.path(tables_dir, "filtered_counts_matrix.csv")
metadata_file <- file.path(tables_dir, "filtered_sample_metadata.csv")

if (!file.exists(results_file) || !file.exists(counts_file) || !file.exists(metadata_file)) {
  stop("Expected tables are missing under `", tables_dir, "`.", call. = FALSE)
}

results_df <- utils::read.csv(results_file, check.names = FALSE)
counts_df <- utils::read.csv(counts_file, check.names = FALSE)
metadata_df <- utils::read.csv(metadata_file, check.names = FALSE)

if (!all(c("Gene", "fit_engine", "convergence_flag") %in% names(results_df))) {
  stop("Results table is missing required columns.", call. = FALSE)
}
if (!all(c("Gene", "SampleName", "Strain") %in% c(names(counts_df)[1], names(metadata_df)))) {
  stop("Filtered counts/metadata tables do not match expected structure.", call. = FALSE)
}

sample_names <- colnames(counts_df)[-1]
metadata_df <- metadata_df[match(sample_names, metadata_df$SampleName), , drop = FALSE]
if (!identical(sample_names, metadata_df$SampleName)) {
  stop("Filtered count columns do not align to filtered sample metadata.", call. = FALSE)
}

count_matrix <- as.matrix(counts_df[, -1, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
gene_ids <- as.character(counts_df[[1]])

strain_levels <- unique(as.character(metadata_df$Strain))
if (!all(c("N2", "daf2") %in% strain_levels)) {
  stop("Expected N2 and daf2 strains in filtered sample metadata.", call. = FALSE)
}

nonzero_by_strain <- list()
total_by_strain <- list()
max_by_strain <- list()
for (strain in c("N2", "daf2")) {
  idx <- metadata_df$Strain == strain
  strain_counts <- count_matrix[, idx, drop = FALSE]
  nonzero_by_strain[[strain]] <- rowSums(strain_counts > 0)
  total_by_strain[[strain]] <- rowSums(strain_counts)
  max_by_strain[[strain]] <- apply(strain_counts, 1, max)
}

gene_metrics <- data.frame(
  Gene = gene_ids,
  nonzero_samples_N2 = nonzero_by_strain[["N2"]],
  nonzero_samples_daf2 = nonzero_by_strain[["daf2"]],
  total_count_N2 = total_by_strain[["N2"]],
  total_count_daf2 = total_by_strain[["daf2"]],
  max_count_N2 = max_by_strain[["N2"]],
  max_count_daf2 = max_by_strain[["daf2"]],
  stringsAsFactors = FALSE
)

results_df <- merge(results_df, gene_metrics, by = "Gene", all.x = TRUE, sort = FALSE)
results_df <- results_df[match(gene_ids, results_df$Gene), , drop = FALSE]

primary_ok <- results_df$fit_engine == "gamlss" &
  results_df$convergence_flag == "success" &
  (is.na(results_df$failure_reason_simple) | results_df$failure_reason_simple == "none")

mean_effect_pool <- abs(results_df$mean_effect[primary_ok & is.finite(results_df$mean_effect)])
disp_effect_pool <- abs(results_df$disp_effect[primary_ok & is.finite(results_df$disp_effect)])

mean_effect_threshold <- max(log(100), safe_quantile(mean_effect_pool, prob = 0.999, fallback = log(100)))
disp_effect_threshold <- max(log(100), safe_quantile(disp_effect_pool, prob = 0.999, fallback = log(100)))

results_df$flag_fallback_engine <- is.na(results_df$fit_engine) | results_df$fit_engine != "gamlss"
results_df$flag_fit_issue <- results_df$convergence_flag != "success" |
  (!is.na(results_df$warning_text) & nzchar(results_df$warning_text)) |
  (!is.na(results_df$failure_reason_simple) & results_df$failure_reason_simple != "none")
results_df$flag_extreme_mean_effect <- is.finite(results_df$mean_effect) &
  abs(results_df$mean_effect) >= mean_effect_threshold
results_df$flag_extreme_disp_effect <- is.finite(results_df$disp_effect) &
  abs(results_df$disp_effect) >= disp_effect_threshold
results_df$flag_nonfinite_statistics <- !is.finite(results_df$mean_p) |
  !is.finite(results_df$mean_FDR) |
  !is.finite(results_df$disp_p) |
  !is.finite(results_df$disp_FDR)
results_df$flag_zero_adjusted_p <- (is.finite(results_df$mean_FDR) & results_df$mean_FDR == 0) |
  (is.finite(results_df$disp_FDR) & results_df$disp_FDR == 0)
results_df$flag_sparse_gene <- results_df$nonzero_samples_N2 < 10 |
  results_df$nonzero_samples_daf2 < 10 |
  results_df$prop_zero_N2 > 0.95 |
  results_df$prop_zero_daf2 > 0.95 |
  (results_df$mean_count_N2 < 1 & results_df$mean_count_daf2 < 1)
results_df$flag_extreme_sparse_combo <- results_df$flag_sparse_gene &
  (results_df$flag_extreme_mean_effect |
     results_df$flag_extreme_disp_effect |
     results_df$flag_nonfinite_statistics |
     results_df$flag_zero_adjusted_p)

flag_columns <- c(
  "flag_fallback_engine",
  "flag_fit_issue",
  "flag_extreme_mean_effect",
  "flag_extreme_disp_effect",
  "flag_nonfinite_statistics",
  "flag_zero_adjusted_p",
  "flag_sparse_gene",
  "flag_extreme_sparse_combo"
)

results_df$stability_flag_count <- rowSums(results_df[, flag_columns, drop = FALSE], na.rm = TRUE)
results_df$stability_class <- ifelse(
  results_df$flag_fallback_engine |
    results_df$flag_fit_issue |
    results_df$flag_nonfinite_statistics |
    results_df$flag_extreme_sparse_combo,
  "high_risk",
  ifelse(
    results_df$flag_extreme_mean_effect |
      results_df$flag_extreme_disp_effect |
      results_df$flag_zero_adjusted_p |
      results_df$flag_sparse_gene,
    "review",
    "stable"
  )
)
results_df$recommended_use <- ifelse(
  results_df$stability_class == "high_risk",
  "exclude_from_primary_interpretation",
  ifelse(
    results_df$stability_class == "review",
    "review_manually_before_high-confidence_use",
    "retain_for_primary_interpretation"
  )
)

results_df$stability_reasons <- vapply(
  seq_len(nrow(results_df)),
  function(i) {
    collapse_flags(
      flag_names = sub("^flag_", "", flag_columns),
      flag_values = as.list(results_df[i, flag_columns, drop = FALSE])
    )
  },
  FUN.VALUE = character(1)
)

annotated_file <- file.path(tables_dir, "per_gene_results_annotated.csv")
utils::write.csv(results_df, annotated_file, row.names = FALSE)

flag_summary <- data.frame(
  flag = sub("^flag_", "", flag_columns),
  n_genes = vapply(flag_columns, function(col) sum(results_df[[col]], na.rm = TRUE), FUN.VALUE = numeric(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(flag_summary, file.path(tables_dir, "stability_flag_summary.csv"), row.names = FALSE)

class_summary <- data.frame(
  stability_class = names(sort(table(results_df$stability_class), decreasing = TRUE)),
  n_genes = as.integer(sort(table(results_df$stability_class), decreasing = TRUE)),
  stringsAsFactors = FALSE
)
utils::write.csv(class_summary, file.path(tables_dir, "stability_class_summary.csv"), row.names = FALSE)

high_risk_df <- results_df[results_df$stability_class == "high_risk", , drop = FALSE]
review_df <- results_df[results_df$stability_class == "review", , drop = FALSE]
stable_df <- results_df[results_df$stability_class == "stable", , drop = FALSE]

rank_score <- pmin(
  ifelse(is.na(results_df$mean_FDR), Inf, results_df$mean_FDR),
  ifelse(is.na(results_df$disp_FDR), Inf, results_df$disp_FDR)
)
results_df$rank_score <- rank_score
high_risk_df <- results_df[results_df$stability_class == "high_risk", , drop = FALSE]
review_df <- results_df[results_df$stability_class == "review", , drop = FALSE]
stable_df <- results_df[results_df$stability_class == "stable", , drop = FALSE]

high_risk_df <- high_risk_df[order(-high_risk_df$stability_flag_count, high_risk_df$rank_score), , drop = FALSE]
review_df <- review_df[order(-review_df$stability_flag_count, review_df$rank_score), , drop = FALSE]

utils::write.csv(high_risk_df, file.path(tables_dir, "high_risk_genes.csv"), row.names = FALSE)
utils::write.csv(review_df, file.path(tables_dir, "review_genes.csv"), row.names = FALSE)
utils::write.csv(stable_df, file.path(tables_dir, "stable_genes.csv"), row.names = FALSE)

report_lines <- c(
  "# Stability Annotation",
  "",
  paste0("- Run directory: `", run_dir, "`"),
  paste0("- Genes annotated: ", nrow(results_df)),
  paste0("- `gamlss`-success reference pool size: ", sum(primary_ok, na.rm = TRUE)),
  paste0("- Extreme mean-effect threshold: abs(effect) >= ", signif(mean_effect_threshold, 4)),
  paste0("- Extreme dispersion-effect threshold: abs(effect) >= ", signif(disp_effect_threshold, 4)),
  "",
  "## Stability classes",
  paste0("- stable: ", sum(results_df$stability_class == "stable")),
  paste0("- review: ", sum(results_df$stability_class == "review")),
  paste0("- high_risk: ", sum(results_df$stability_class == "high_risk")),
  "",
  "## Flag counts",
  paste0("- fallback_engine: ", sum(results_df$flag_fallback_engine)),
  paste0("- fit_issue: ", sum(results_df$flag_fit_issue)),
  paste0("- extreme_mean_effect: ", sum(results_df$flag_extreme_mean_effect)),
  paste0("- extreme_disp_effect: ", sum(results_df$flag_extreme_disp_effect)),
  paste0("- nonfinite_statistics: ", sum(results_df$flag_nonfinite_statistics)),
  paste0("- zero_adjusted_p: ", sum(results_df$flag_zero_adjusted_p)),
  paste0("- sparse_gene: ", sum(results_df$flag_sparse_gene)),
  paste0("- extreme_sparse_combo: ", sum(results_df$flag_extreme_sparse_combo)),
  "",
  "## Interpretation",
  "- `stable`: no major model-quality or effect-size warning signals from this pass.",
  "- `review`: fitted successfully, but one or more heuristic warning flags suggest manual inspection before high-confidence use.",
  "- `high_risk`: fallback fit, explicit fit issue, non-finite statistics, or extreme effects on sparse genes. These genes are poor candidates for primary biological interpretation without additional validation."
)
writeLines(report_lines, con = file.path(run_dir, "STABILITY_REVIEW.md"))

message("Annotated results written to: ", annotated_file)

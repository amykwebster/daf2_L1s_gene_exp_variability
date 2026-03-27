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

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

capture_conditions <- function(expr) {
  warnings <- character(0)
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(
    value = value,
    warnings = unique(warnings),
    error = inherits(value, "error")
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
run_dir <- normalizePath(args[["run-dir"]] %||% ".", mustWork = TRUE)

tables_dir <- file.path(run_dir, "tables")
robustness_file <- file.path(
  tables_dir,
  "dispersion_replicate_robustness_stable_gamlss_dispersion_sig.csv"
)
counts_file <- file.path(tables_dir, "filtered_counts_matrix.csv")
metadata_file <- file.path(tables_dir, "filtered_sample_metadata.csv")

if (!file.exists(robustness_file)) {
  stop("Robustness table not found: ", robustness_file, call. = FALSE)
}

robustness_df <- utils::read.csv(robustness_file, check.names = FALSE)
counts_df <- utils::read.csv(counts_file, check.names = FALSE)
metadata_df <- utils::read.csv(metadata_file, check.names = FALSE)

concordant_df <- robustness_df[
  robustness_df$robustness_class == "concordant_all_replicates",
  ,
  drop = FALSE
]
if (nrow(concordant_df) == 0L) {
  stop("No concordant genes found in the robustness table.", call. = FALSE)
}

sample_names <- colnames(counts_df)[-1]
metadata_df <- metadata_df[match(sample_names, metadata_df$SampleName), , drop = FALSE]
if (!identical(sample_names, metadata_df$SampleName)) {
  stop("Filtered counts and metadata are not aligned.", call. = FALSE)
}

metadata_df$Strain <- factor(metadata_df$Strain, levels = c("N2", "daf2"))
metadata_df$Replicate <- factor(metadata_df$Replicate)
metadata_df$log_library_size <- log(pmax(metadata_df$library_size, 1))

count_matrix <- as.matrix(counts_df[, -1, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
rownames(count_matrix) <- counts_df[[1]]

out_file <- file.path(tables_dir, "concordant_dispersion_genes_residual_dispersion_long.csv")
summary_file <- file.path(tables_dir, "concordant_dispersion_genes_residual_dispersion_summary.csv")
failed_file <- file.path(tables_dir, "concordant_dispersion_genes_residual_dispersion_failures.csv")

results_long <- vector("list", nrow(concordant_df))
failed_rows <- list()

progress <- utils::txtProgressBar(min = 0, max = nrow(concordant_df), style = 3)
for (i in seq_len(nrow(concordant_df))) {
  gene <- concordant_df$Gene[[i]]
  y <- as.numeric(count_matrix[gene, ])

  working_df <- metadata_df
  working_df$count <- y

  fit_capture <- capture_conditions(
    MASS::glm.nb(
      count ~ Strain + Replicate + offset(log_library_size),
      data = working_df,
      control = stats::glm.control(maxit = 100)
    )
  )

  if (fit_capture$error) {
    failed_rows[[length(failed_rows) + 1L]] <- data.frame(
      Gene = gene,
      error_message = fit_capture$value$message,
      warning_text = if (length(fit_capture$warnings) == 0L) NA_character_ else paste(fit_capture$warnings, collapse = " | "),
      stringsAsFactors = FALSE
    )
    utils::setTxtProgressBar(progress, i)
    next
  }

  mean_fit <- fit_capture$value
  residual_values <- pmax(stats::residuals(mean_fit, type = "pearson")^2, 1e-6)

  results_long[[i]] <- data.frame(
    Gene = gene,
    SampleName = working_df$SampleName,
    Strain = working_df$Strain,
    Replicate = working_df$Replicate,
    Barcode = working_df$Barcode,
    library_size = working_df$library_size,
    raw_count = y,
    fitted_mean = stats::fitted(mean_fit),
    pearson_residual = stats::residuals(mean_fit, type = "pearson"),
    squared_pearson_residual = residual_values,
    log10_squared_pearson_residual = log10(residual_values + 1e-6),
    disp_effect = concordant_df$disp_effect[[i]],
    disp_FDR = concordant_df$disp_FDR[[i]],
    global_disp_direction = concordant_df$global_disp_direction[[i]],
    robustness_class = concordant_df$robustness_class[[i]],
    stringsAsFactors = FALSE
  )

  utils::setTxtProgressBar(progress, i)
}
close(progress)

results_long <- do.call(rbind, results_long[!vapply(results_long, is.null, logical(1))])
utils::write.csv(results_long, out_file, row.names = FALSE)

summary_df <- aggregate(
  cbind(
    mean_raw_count = raw_count,
    mean_fitted_mean = fitted_mean,
    mean_squared_pearson_residual = squared_pearson_residual
  ) ~ Gene + Strain + Replicate + global_disp_direction + disp_effect + disp_FDR,
  data = results_long,
  FUN = mean
)
utils::write.csv(summary_df, summary_file, row.names = FALSE)

failed_df <- if (length(failed_rows) == 0L) {
  data.frame(Gene = character(0), error_message = character(0), warning_text = character(0), stringsAsFactors = FALSE)
} else {
  do.call(rbind, failed_rows)
}
utils::write.csv(failed_df, failed_file, row.names = FALSE)

message("Wrote residual-dispersion values to: ", out_file)

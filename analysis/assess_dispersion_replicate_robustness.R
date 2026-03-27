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

extract_glm_term <- function(model, term_name) {
  coef_table <- coef(summary(model))
  if (!term_name %in% rownames(coef_table)) {
    return(list(effect = NA_real_, p_value = NA_real_, found = FALSE))
  }
  p_col <- intersect(c("Pr(>|z|)", "Pr(>|t|)"), colnames(coef_table))
  list(
    effect = unname(coef_table[term_name, "Estimate"]),
    p_value = if (length(p_col) == 1L) unname(coef_table[term_name, p_col]) else NA_real_,
    found = TRUE
  )
}

classify_direction <- function(effect, tol = 1e-12) {
  if (!is.finite(effect)) {
    return("not_estimable")
  }
  if (effect > tol) {
    return("higher_disp_in_daf2")
  }
  if (effect < -tol) {
    return("lower_disp_in_daf2")
  }
  "no_difference"
}

sanitize_text <- function(x) {
  if (length(x) == 0L || all(is.na(x)) || !any(nzchar(x))) {
    NA_character_
  } else {
    paste(unique(x[!is.na(x) & nzchar(x)]), collapse = " | ")
  }
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
run_dir <- normalizePath(args[["run-dir"]] %||% ".", mustWork = TRUE)
selection <- args[["selection"]] %||% "stable_gamlss_dispersion_sig"

tables_dir <- file.path(run_dir, "tables")
annotated_file <- file.path(tables_dir, "per_gene_results_annotated.csv")
counts_file <- file.path(tables_dir, "filtered_counts_matrix.csv")
metadata_file <- file.path(tables_dir, "filtered_sample_metadata.csv")

if (!file.exists(annotated_file)) {
  stop("Annotated results not found at `", annotated_file, "`.", call. = FALSE)
}

results_df <- utils::read.csv(annotated_file, check.names = FALSE)
counts_df <- utils::read.csv(counts_file, check.names = FALSE)
metadata_df <- utils::read.csv(metadata_file, check.names = FALSE)

gene_col <- names(counts_df)[1]
counts_df$gene_occurrence_id <- ave(seq_len(nrow(counts_df)), counts_df[[gene_col]], FUN = seq_along)
results_df$gene_occurrence_id <- ave(seq_len(nrow(results_df)), results_df$Gene, FUN = seq_along)
counts_df$gene_match_id <- paste(counts_df[[gene_col]], counts_df$gene_occurrence_id, sep = "__occ__")
results_df$gene_match_id <- paste(results_df$Gene, results_df$gene_occurrence_id, sep = "__occ__")

sample_names <- setdiff(colnames(counts_df), c(gene_col, "gene_occurrence_id", "gene_match_id"))
metadata_df <- metadata_df[match(sample_names, metadata_df$SampleName), , drop = FALSE]
if (!identical(sample_names, metadata_df$SampleName)) {
  stop("Filtered counts and metadata are not aligned.", call. = FALSE)
}

metadata_df$Strain <- stats::relevel(factor(metadata_df$Strain), ref = "N2")
metadata_df$Replicate <- factor(metadata_df$Replicate)
metadata_df$log_library_size <- log(pmax(metadata_df$library_size, 1))

count_matrix <- as.matrix(counts_df[, sample_names, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
rownames(count_matrix) <- counts_df$gene_match_id

if (selection == "stable_gamlss_dispersion_sig") {
  selected_results <- results_df[
    results_df$stability_class == "stable" &
      results_df$fit_engine == "gamlss" &
      !is.na(results_df$disp_FDR) &
      results_df$disp_FDR < 0.05,
    ,
    drop = FALSE
  ]
} else if (selection == "all_dispersion_sig") {
  selected_results <- results_df[
    !is.na(results_df$disp_FDR) &
      results_df$disp_FDR < 0.05,
    ,
    drop = FALSE
  ]
} else {
  stop("Unsupported selection: ", selection, call. = FALSE)
}

if (nrow(selected_results) == 0L) {
  stop("No genes matched the requested selection.", call. = FALSE)
}

selected_results$count_row_index <- match(selected_results$gene_match_id, rownames(count_matrix))
if (any(is.na(selected_results$count_row_index))) {
  stop("Some selected genes could not be aligned to the filtered count matrix.", call. = FALSE)
}
selected_results <- selected_results[order(selected_results$count_row_index), , drop = FALSE]

replicate_levels <- levels(metadata_df$Replicate)
results_list <- vector("list", nrow(selected_results))

progress <- utils::txtProgressBar(min = 0, max = nrow(selected_results), style = 3)
for (i in seq_len(nrow(selected_results))) {
  gene_row <- selected_results[i, , drop = FALSE]
  gene_id <- gene_row$Gene[[1]]
  gene_match_id <- gene_row$gene_match_id[[1]]
  count_row_index <- gene_row$count_row_index[[1]]
  y <- as.numeric(count_matrix[count_row_index, ])

  working_df <- metadata_df
  working_df$count_response <- y

  mean_capture <- capture_conditions(
    MASS::glm.nb(
      count_response ~ Strain + Replicate + offset(log_library_size),
      data = working_df,
      control = stats::glm.control(maxit = 100)
    )
  )

  if (mean_capture$error) {
    results_list[[i]] <- data.frame(
      Gene = gene_id,
      gene_match_id = gene_match_id,
      global_disp_direction = classify_direction(gene_row$disp_effect[[1]]),
      interaction_p = NA_real_,
      interaction_support = "model_failure",
      n_replicates_evaluable = 0L,
      n_replicates_same_direction = 0L,
      n_replicates_opposite_direction = 0L,
      replicate_I_direction = NA_character_,
      replicate_II_direction = NA_character_,
      replicate_III_direction = NA_character_,
      replicate_I_effect = NA_real_,
      replicate_II_effect = NA_real_,
      replicate_III_effect = NA_real_,
      replicate_I_p = NA_real_,
      replicate_II_p = NA_real_,
      replicate_III_p = NA_real_,
      robustness_class = "model_failure",
      robustness_reason = paste0("mean_model_failed:", mean_capture$value$message),
      model_warning_text = sanitize_text(mean_capture$warnings),
      stringsAsFactors = FALSE
    )
    utils::setTxtProgressBar(progress, i)
    next
  }

  mean_fit <- mean_capture$value
  working_df$squared_pearson_resid <- pmax(stats::residuals(mean_fit, type = "pearson")^2, 1e-6)

  disp_add_capture <- capture_conditions(
    stats::glm(
      squared_pearson_resid ~ Strain + Replicate,
      family = stats::Gamma(link = "log"),
      data = working_df
    )
  )
  disp_int_capture <- capture_conditions(
    stats::glm(
      squared_pearson_resid ~ Strain * Replicate,
      family = stats::Gamma(link = "log"),
      data = working_df
    )
  )

  interaction_p <- NA_real_
  interaction_support <- "not_tested"

  if (!disp_add_capture$error && !disp_int_capture$error) {
    lrt <- tryCatch(
      stats::anova(disp_add_capture$value, disp_int_capture$value, test = "LRT"),
      error = function(e) NULL
    )
    if (!is.null(lrt) && nrow(lrt) >= 2L && "Pr(>Chi)" %in% colnames(lrt)) {
      interaction_p <- unname(lrt[2, "Pr(>Chi)"])
      interaction_support <- if (is.finite(interaction_p) && interaction_p < 0.05) {
        "interaction_supported"
      } else {
        "no_interaction_support"
      }
    }
  }

  global_direction <- classify_direction(gene_row$disp_effect[[1]])
  replicate_directions <- setNames(rep(NA_character_, length(replicate_levels)), replicate_levels)
  replicate_pvalues <- setNames(rep(NA_real_, length(replicate_levels)), replicate_levels)
  replicate_effects <- setNames(rep(NA_real_, length(replicate_levels)), replicate_levels)
  evaluable <- logical(length(replicate_levels))

  for (rep_name in replicate_levels) {
    rep_idx <- working_df$Replicate == rep_name
    rep_df <- working_df[rep_idx, , drop = FALSE]

    strain_counts <- table(rep_df$Strain)
    if (length(strain_counts) < 2L || any(strain_counts < 5L)) {
      replicate_directions[[rep_name]] <- "insufficient_samples"
      next
    }

    rep_capture <- capture_conditions(
      stats::glm(
        squared_pearson_resid ~ Strain,
        family = stats::Gamma(link = "log"),
        data = rep_df
      )
    )

    if (rep_capture$error) {
      replicate_directions[[rep_name]] <- "fit_failed"
      next
    }

    rep_term <- extract_glm_term(rep_capture$value, "Straindaf2")
    replicate_effects[[rep_name]] <- rep_term$effect
    replicate_pvalues[[rep_name]] <- rep_term$p_value
    replicate_directions[[rep_name]] <- classify_direction(rep_term$effect)
    evaluable[[which(replicate_levels == rep_name)]] <- rep_term$found
  }

  evaluable_dirs <- replicate_directions[evaluable]
  n_same <- sum(evaluable_dirs == global_direction, na.rm = TRUE)
  n_opposite <- sum(
    evaluable_dirs %in% c("higher_disp_in_daf2", "lower_disp_in_daf2") &
      evaluable_dirs != global_direction,
    na.rm = TRUE
  )
  n_evaluable <- sum(evaluable)

  robustness_class <- if (n_evaluable == 0L) {
    "insufficient_replicate_information"
  } else if (!is.na(interaction_p) && interaction_p < 0.05) {
    "interaction_supported"
  } else if (n_same == n_evaluable && n_evaluable == length(replicate_levels)) {
    "concordant_all_replicates"
  } else if (n_same >= 2L && n_opposite == 0L) {
    "concordant_majority_replicates"
  } else if (n_opposite >= 1L) {
    "mixed_direction_across_replicates"
  } else {
    "unclear_replicate_pattern"
  }

  robustness_reason <- switch(
    robustness_class,
    concordant_all_replicates = "All evaluable replicates show the same dispersion direction as the global daf2 effect.",
    concordant_majority_replicates = "At least two replicates agree with the global daf2 dispersion direction and none oppose it.",
    mixed_direction_across_replicates = "At least one replicate shows an opposite dispersion direction from the global daf2 effect.",
    interaction_supported = "A Strain-by-Replicate interaction is supported in the residual-dispersion model, suggesting heterogeneity across replicates.",
    insufficient_replicate_information = "Too few replicate-specific fits could be estimated to evaluate concordance.",
    unclear_replicate_pattern = "Replicate-level pattern was not fully concordant and did not meet the interaction threshold.",
    model_failure = "Mean or residual-dispersion model failed.",
    "Unclassified result."
  )

  results_list[[i]] <- data.frame(
    Gene = gene_id,
    gene_match_id = gene_match_id,
    global_disp_direction = global_direction,
    interaction_p = interaction_p,
    interaction_support = interaction_support,
    n_replicates_evaluable = n_evaluable,
    n_replicates_same_direction = n_same,
    n_replicates_opposite_direction = n_opposite,
    replicate_I_direction = replicate_directions[["I"]] %||% NA_character_,
    replicate_II_direction = replicate_directions[["II"]] %||% NA_character_,
    replicate_III_direction = replicate_directions[["III"]] %||% NA_character_,
    replicate_I_effect = replicate_effects[["I"]] %||% NA_real_,
    replicate_II_effect = replicate_effects[["II"]] %||% NA_real_,
    replicate_III_effect = replicate_effects[["III"]] %||% NA_real_,
    replicate_I_p = replicate_pvalues[["I"]] %||% NA_real_,
    replicate_II_p = replicate_pvalues[["II"]] %||% NA_real_,
    replicate_III_p = replicate_pvalues[["III"]] %||% NA_real_,
    robustness_class = robustness_class,
    robustness_reason = robustness_reason,
    model_warning_text = sanitize_text(c(mean_capture$warnings, disp_add_capture$warnings, disp_int_capture$warnings)),
    stringsAsFactors = FALSE
  )

  utils::setTxtProgressBar(progress, i)
}
close(progress)

robustness_df <- do.call(rbind, results_list)
robustness_df <- merge(selected_results, robustness_df, by = c("Gene", "gene_match_id"), all.x = TRUE, sort = FALSE)
robustness_df <- robustness_df[match(selected_results$gene_match_id, robustness_df$gene_match_id), , drop = FALSE]

selection_tag <- if (selection == "stable_gamlss_dispersion_sig") {
  "stable_gamlss_dispersion_sig"
} else {
  "all_dispersion_sig"
}

out_main <- file.path(tables_dir, paste0("dispersion_replicate_robustness_", selection_tag, ".csv"))
utils::write.csv(robustness_df, out_main, row.names = FALSE)

summary_table <- data.frame(
  robustness_class = names(sort(table(robustness_df$robustness_class), decreasing = TRUE)),
  n_genes = as.integer(sort(table(robustness_df$robustness_class), decreasing = TRUE)),
  stringsAsFactors = FALSE
)
out_summary <- file.path(tables_dir, paste0("dispersion_replicate_robustness_summary_", selection_tag, ".csv"))
utils::write.csv(summary_table, out_summary, row.names = FALSE)

readme_lines <- c(
  "# Replicate robustness of dispersion hits",
  "",
  paste0("- Selection: `", selection, "`"),
  paste0("- Genes assessed: ", nrow(robustness_df)),
  "- Mean model used for residualization: `count ~ Strain + Replicate + offset(log_library_size)` fit by `glm.nb`.",
  "- Residual dispersion proxy: squared Pearson residuals from the mean model.",
  "- Replicate-level direction: estimated from `Gamma(log)` models of squared Pearson residuals within each replicate (`~ Strain`).",
  "- Interaction screen: additive residual-dispersion model (`~ Strain + Replicate`) compared against an interaction model (`~ Strain * Replicate`) using an LRT.",
  "",
  "## Robustness classes",
  "- `concordant_all_replicates`: all evaluable replicates showed the same dispersion direction as the global daf2 effect and no interaction support was detected.",
  "- `concordant_majority_replicates`: at least two replicates agreed with the global direction and none opposed it.",
  "- `mixed_direction_across_replicates`: at least one replicate showed the opposite direction.",
  "- `interaction_supported`: replicate heterogeneity was supported by the interaction model.",
  "- `insufficient_replicate_information`: too few replicate-specific fits could be estimated.",
  "- `unclear_replicate_pattern`: neither clearly concordant nor interaction-supported."
)
writeLines(
  readme_lines,
  con = file.path(run_dir, paste0("DISPERSION_REPLICATE_ROBUSTNESS_", toupper(selection_tag), ".md"))
)

message("Wrote robustness table to: ", out_main)

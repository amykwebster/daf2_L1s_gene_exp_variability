#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

local_r_lib <- "r_libs"
if (dir.exists(local_r_lib)) {
  .libPaths(c(normalizePath(local_r_lib), .libPaths()))
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) {
    y
  } else {
    x
  }
}

# ----------------------------- User-editable settings ------------------------
# Edit `user_config` directly for a local run, or pass command-line overrides:
#   Rscript run_per_gene_gamlss.R --counts path/to/counts.csv \
#     --sample-key path/to/SampleKey.csv --output-dir results_dir
#
# Optional:
#   --config path/to/config_overrides.R    # must define `config_overrides`
#   --max-genes 100                        # useful for smoke tests

default_config <- list(
  counts_file = "umi_count_matrix_all_combined_3sets.csv",
  sample_key_file = "SampleNames.csv",
  output_dir = "per_gene_gamlss_output",
  gene_column = "Gene",
  sample_column = "SampleName",
  strain_column = "Strain",
  replicate_column = "Replicate",
  strain_rep_column = "Strain_Rep",
  strain_reference = "N2",
  significance_fdr = 0.05,
  qc = list(
    libsize_mad_cutoff = 3,
    detected_mad_cutoff = 3,
    min_library_size = 1,
    min_detected_genes = 100
  ),
  gene_filter = list(
    min_nonzero_fraction = 0.05,
    min_nonzero_samples = 10,
    min_total_count = 20
  ),
  model = list(
    mean_terms = c("Strain", "Replicate"),
    dispersion_terms = c("Strain", "Replicate"),
    include_library_offset = TRUE,
    preferred_families = c("NBI", "NBII"),
    gamlss_max_cycles = 200,
    gamlss_trace = FALSE,
    fallback_engine = "glm_nb_gamma_dispersion"
  ),
  plotting = list(
    example_genes_per_category = 2
  ),
  debug = list(
    max_genes = Inf
  )
)

user_config <- list()

# ------------------------------- Helper functions ----------------------------
parse_args <- function(args) {
  parsed <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    if (grepl("=", key, fixed = TRUE)) {
      split_key <- strsplit(sub("^--", "", key), "=", fixed = TRUE)[[1]]
      name <- split_key[[1]]
      value <- paste(split_key[-1], collapse = "=")
    } else {
      name <- sub("^--", "", key)
      if (i == length(args)) {
        stop("Missing value for argument: ", key, call. = FALSE)
      }
      i <- i + 1L
      value <- args[[i]]
    }
    parsed[[name]] <- value
    i <- i + 1L
  }
  parsed
}

load_config_overrides <- function(path) {
  env <- new.env(parent = emptyenv())
  sys.source(path, envir = env)
  if (exists("config_overrides", envir = env, inherits = FALSE)) {
    get("config_overrides", envir = env, inherits = FALSE)
  } else if (exists("config", envir = env, inherits = FALSE)) {
    get("config", envir = env, inherits = FALSE)
  } else {
    stop("Config file must define `config_overrides` or `config`.", call. = FALSE)
  }
}

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

stop_if_missing_package <- function(pkg, reason) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Missing required package `", pkg, "`: ", reason,
      "\nInstall it with install.packages(\"", pkg, "\") and rerun.",
      call. = FALSE
    )
  }
}

collapse_reasons <- function(...) {
  reasons <- unlist(list(...), use.names = FALSE)
  reasons <- reasons[!is.na(reasons) & nzchar(reasons)]
  if (length(reasons) == 0L) {
    "retained"
  } else {
    paste(unique(reasons), collapse = "; ")
  }
}

mad_lower_threshold <- function(x, nmads, offset = 1) {
  x_log <- log10(x + offset)
  med <- stats::median(x_log, na.rm = TRUE)
  mad_value <- stats::mad(x_log, center = med, constant = 1, na.rm = TRUE)
  threshold_log <- med - (nmads * mad_value)
  threshold_raw <- max(0, (10^threshold_log) - offset)
  list(
    median_log10 = med,
    mad_log10 = mad_value,
    threshold_log10 = threshold_log,
    threshold_raw = threshold_raw
  )
}

build_model_formula <- function(response, terms, include_offset, available_data) {
  retained_terms <- character(0)
  for (term in terms) {
    if (!term %in% names(available_data)) {
      next
    }
    values <- available_data[[term]]
    n_unique <- length(unique(values[!is.na(values)]))
    if (n_unique > 1L) {
      retained_terms <- c(retained_terms, term)
    }
  }
  if (include_offset) {
    retained_terms <- c(retained_terms, "offset(log_library_size)")
  }
  rhs <- if (length(retained_terms) == 0L) "1" else paste(retained_terms, collapse = " + ")
  stats::as.formula(paste(response, "~", rhs))
}

rhs_only_formula <- function(formula_object) {
  rhs_terms <- attr(stats::terms(formula_object), "term.labels")
  if (length(rhs_terms) == 0L) {
    stats::as.formula("~ 1")
  } else {
    stats::as.formula(paste("~", paste(rhs_terms, collapse = " + ")))
  }
}

find_strain_term <- function(coef_names, strain_column, strain_levels, reference_level) {
  non_reference <- setdiff(strain_levels, reference_level)
  if (length(non_reference) != 1L) {
    return(NA_character_)
  }
  exact_term <- paste0(strain_column, non_reference)
  if (exact_term %in% coef_names) {
    return(exact_term)
  }
  candidates <- grep(paste0("^", strain_column), coef_names, value = TRUE)
  if (length(candidates) == 1L) {
    return(candidates)
  }
  matched <- candidates[grepl(non_reference, candidates, fixed = TRUE)]
  if (length(matched) >= 1L) {
    return(matched[[1]])
  }
  NA_character_
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
  if (is.na(term_name)) {
    return(list(effect = NA_real_, p_value = NA_real_, se = NA_real_, found = FALSE))
  }
  coef_table <- coef(summary(model))
  if (!term_name %in% rownames(coef_table)) {
    return(list(effect = NA_real_, p_value = NA_real_, se = NA_real_, found = FALSE))
  }
  p_col <- intersect(c("Pr(>|z|)", "Pr(>|t|)"), colnames(coef_table))
  list(
    effect = unname(coef_table[term_name, "Estimate"]),
    p_value = if (length(p_col) == 1L) unname(coef_table[term_name, p_col]) else NA_real_,
    se = unname(coef_table[term_name, "Std. Error"]),
    found = TRUE
  )
}

extract_gamlss_coefficients <- function(fit, what) {
  extracted <- tryCatch(stats::coef(fit, what = what), error = function(e) NULL)
  if (!is.null(extracted)) {
    return(extracted)
  }
  fit[[paste0(what, ".coefficients")]] %||% NULL
}

extract_gamlss_vcov <- function(fit, what) {
  extracted <- tryCatch(stats::vcov(fit, what = what), error = function(e) NULL)
  if (!is.null(extracted)) {
    coef_mu <- extract_gamlss_coefficients(fit, "mu")
    coef_sigma <- extract_gamlss_coefficients(fit, "sigma")
    mu_len <- length(coef_mu %||% numeric(0))
    sigma_len <- length(coef_sigma %||% numeric(0))
    total_len <- mu_len + sigma_len

    if (is.matrix(extracted) && nrow(extracted) >= total_len && total_len > 0L) {
      if (what == "mu" && mu_len > 0L) {
        idx <- seq_len(mu_len)
        block <- extracted[idx, idx, drop = FALSE]
        rownames(block) <- names(coef_mu)
        colnames(block) <- names(coef_mu)
        return(block)
      }
      if (what == "sigma" && sigma_len > 0L) {
        idx <- mu_len + seq_len(sigma_len)
        block <- extracted[idx, idx, drop = FALSE]
        rownames(block) <- names(coef_sigma)
        colnames(block) <- names(coef_sigma)
        return(block)
      }
    }
    return(extracted)
  }
  fit[[paste0(what, ".var")]] %||% NULL
}

extract_gamlss_term <- function(fit, what, term_name) {
  if (is.na(term_name)) {
    return(list(effect = NA_real_, p_value = NA_real_, se = NA_real_, found = FALSE))
  }
  coefficients <- extract_gamlss_coefficients(fit, what)
  vcov_mat <- extract_gamlss_vcov(fit, what)
  if (is.null(coefficients) || is.null(vcov_mat) || !term_name %in% names(coefficients)) {
    return(list(effect = NA_real_, p_value = NA_real_, se = NA_real_, found = FALSE))
  }
  vcov_names <- rownames(vcov_mat) %||% colnames(vcov_mat)
  if (is.null(vcov_names) || !term_name %in% vcov_names) {
    return(list(effect = NA_real_, p_value = NA_real_, se = NA_real_, found = FALSE))
  }
  se <- sqrt(diag(vcov_mat))[term_name]
  estimate <- unname(coefficients[term_name])
  if (!is.finite(se) || se <= 0) {
    p_value <- NA_real_
  } else {
    z_value <- estimate / se
    p_value <- 2 * stats::pnorm(abs(z_value), lower.tail = FALSE)
  }
  list(
    effect = estimate,
    p_value = p_value,
    se = unname(se),
    found = TRUE
  )
}

safe_variance <- function(x) {
  if (length(x) <= 1L) {
    return(NA_real_)
  }
  stats::var(x)
}

fit_gene_gamlss <- function(y, sample_df, mean_formula, sigma_formula, strain_term, config) {
  if (!requireNamespace("gamlss", quietly = TRUE) ||
      !requireNamespace("gamlss.dist", quietly = TRUE)) {
    return(NULL)
  }

  sample_df$count_response <- as.numeric(y)
  failure_notes <- character(0)
  data_symbol <- ".codex_gamlss_fit_data"
  assign(data_symbol, sample_df, envir = .GlobalEnv)
  on.exit(
    if (exists(data_symbol, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = data_symbol, envir = .GlobalEnv)
    },
    add = TRUE
  )

  for (family_name in config$model$preferred_families) {
    family_object <- tryCatch(
      get(family_name, envir = asNamespace("gamlss.dist")),
      error = function(e) NULL
    )
    if (is.null(family_object)) {
      failure_notes <- c(failure_notes, paste0("family_not_found:", family_name))
      next
    }

    condition_capture <- capture_conditions(
      eval(
        substitute(
          gamlss::gamlss(
            formula = MEAN_FORMULA,
            sigma.formula = SIGMA_FORMULA,
            family = FAMILY_OBJECT,
            data = DATA_OBJECT,
            control = gamlss::gamlss.control(
              n.cyc = MAX_CYCLES,
              trace = TRACE_FLAG
            )
          ),
          list(
            MEAN_FORMULA = mean_formula,
            SIGMA_FORMULA = sigma_formula,
            FAMILY_OBJECT = family_object,
            DATA_OBJECT = as.name(data_symbol),
            MAX_CYCLES = config$model$gamlss_max_cycles,
            TRACE_FLAG = config$model$gamlss_trace
          )
        ),
        envir = .GlobalEnv
      )
    )

    if (condition_capture$error) {
      failure_notes <- c(
        failure_notes,
        paste0("gamlss_", family_name, "_error:", condition_capture$value$message)
      )
      next
    }

    fit <- condition_capture$value
    mean_test <- extract_gamlss_term(fit, "mu", strain_term)
    disp_test <- extract_gamlss_term(fit, "sigma", strain_term)
    warnings <- condition_capture$warnings
    terms_found <- isTRUE(mean_test$found) && isTRUE(disp_test$found)

    fit_result <- list(
      engine = "gamlss",
      family = family_name,
      mean_effect = mean_test$effect,
      mean_p = mean_test$p_value,
      disp_effect = disp_test$effect,
      disp_p = disp_test$p_value,
      mean_term_found = mean_test$found,
      disp_term_found = disp_test$found,
      convergence_flag = if (length(warnings) == 0L) "success" else "warning",
      warning_text = if (length(warnings) == 0L) NA_character_ else paste(warnings, collapse = " | "),
      failure_reason = if (terms_found) {
        NA_character_
      } else {
        "gamlss_term_extraction_failed"
      }
    )

    if (terms_found) {
      return(fit_result)
    }
    failure_notes <- c(failure_notes, fit_result$failure_reason)
  }

  list(
    engine = "gamlss",
    family = NA_character_,
    mean_effect = NA_real_,
    mean_p = NA_real_,
    disp_effect = NA_real_,
    disp_p = NA_real_,
    mean_term_found = FALSE,
    disp_term_found = FALSE,
    convergence_flag = "failed",
    warning_text = if (length(failure_notes) == 0L) NA_character_ else paste(unique(failure_notes), collapse = " | "),
    failure_reason = if (length(failure_notes) == 0L) "gamlss_failed" else paste(unique(failure_notes), collapse = " | ")
  )
}

fit_gene_fallback <- function(y, sample_df, mean_formula, disp_formula, strain_term) {
  sample_df$count_response <- y

  mean_capture <- capture_conditions(
    MASS::glm.nb(
      formula = mean_formula,
      data = sample_df,
      control = stats::glm.control(maxit = 100)
    )
  )

  if (mean_capture$error) {
    return(list(
      engine = "fallback",
      family = "glm.nb + gamma_dispersion",
      mean_effect = NA_real_,
      mean_p = NA_real_,
      disp_effect = NA_real_,
      disp_p = NA_real_,
      convergence_flag = "failed",
      warning_text = if (length(mean_capture$warnings) == 0L) NA_character_ else paste(mean_capture$warnings, collapse = " | "),
      failure_reason = paste0("fallback_mean_error:", mean_capture$value$message)
    ))
  }

  mean_fit <- mean_capture$value
  mean_test <- extract_glm_term(mean_fit, strain_term)

  squared_residuals <- pmax(stats::residuals(mean_fit, type = "pearson")^2, 1e-6)
  sample_df$squared_pearson_resid <- squared_residuals
  disp_formula_fallback <- stats::update(disp_formula, squared_pearson_resid ~ .)

  disp_capture <- capture_conditions(
    stats::glm(
      formula = disp_formula_fallback,
      family = stats::Gamma(link = "log"),
      data = sample_df
    )
  )

  if (disp_capture$error) {
    return(list(
      engine = "fallback",
      family = "glm.nb + gamma_dispersion",
      mean_effect = mean_test$effect,
      mean_p = mean_test$p_value,
      disp_effect = NA_real_,
      disp_p = NA_real_,
      convergence_flag = "partial_failure",
      warning_text = paste(c(mean_capture$warnings, disp_capture$warnings), collapse = " | "),
      failure_reason = paste0("fallback_dispersion_error:", disp_capture$value$message)
    ))
  }

  disp_fit <- disp_capture$value
  disp_test <- extract_glm_term(disp_fit, strain_term)
  warnings <- unique(c(mean_capture$warnings, disp_capture$warnings))

  list(
    engine = "fallback",
    family = "glm.nb + gamma_dispersion",
    mean_effect = mean_test$effect,
    mean_p = mean_test$p_value,
    disp_effect = disp_test$effect,
    disp_p = disp_test$p_value,
    convergence_flag = if (length(warnings) == 0L) "success" else "warning",
    warning_text = if (length(warnings) == 0L) NA_character_ else paste(warnings, collapse = " | "),
    failure_reason = if (mean_test$found && disp_test$found) NA_character_ else "fallback_term_extraction_failed"
  )
}

fit_single_gene <- function(y, sample_df, mean_formula, sigma_formula, strain_term, config) {
  sample_df$count_response <- as.numeric(y)

  if (all(sample_df$count_response == 0)) {
    return(list(
      mean_effect = NA_real_,
      mean_p = NA_real_,
      disp_effect = NA_real_,
      disp_p = NA_real_,
      convergence_flag = "failed",
      model_family = NA_character_,
      fit_engine = NA_character_,
      failure_reason = "all_zero_after_filtering",
      warning_text = NA_character_
    ))
  }

  if (length(unique(sample_df$count_response)) == 1L) {
    return(list(
      mean_effect = NA_real_,
      mean_p = NA_real_,
      disp_effect = NA_real_,
      disp_p = NA_real_,
      convergence_flag = "failed",
      model_family = NA_character_,
      fit_engine = NA_character_,
      failure_reason = "no_expression_variation_after_filtering",
      warning_text = NA_character_
    ))
  }

  gamlss_result <- fit_gene_gamlss(
    y = sample_df$count_response,
    sample_df = sample_df,
    mean_formula = mean_formula,
    sigma_formula = sigma_formula,
    strain_term = strain_term,
    config = config
  )

  if (!is.null(gamlss_result) &&
      !identical(gamlss_result$convergence_flag, "failed") &&
      is.na(gamlss_result$failure_reason)) {
    return(list(
      mean_effect = gamlss_result$mean_effect,
      mean_p = gamlss_result$mean_p,
      disp_effect = gamlss_result$disp_effect,
      disp_p = gamlss_result$disp_p,
      convergence_flag = gamlss_result$convergence_flag,
      model_family = gamlss_result$family,
      fit_engine = gamlss_result$engine,
      failure_reason = gamlss_result$failure_reason,
      warning_text = gamlss_result$warning_text
    ))
  }

  fallback_result <- fit_gene_fallback(
    y = sample_df$count_response,
    sample_df = sample_df,
    mean_formula = mean_formula,
    disp_formula = sigma_formula,
    strain_term = strain_term
  )

  combined_failure <- unique(c(
    gamlss_result$failure_reason %||% character(0),
    fallback_result$failure_reason %||% character(0)
  ))
  combined_failure <- combined_failure[!is.na(combined_failure) & nzchar(combined_failure)]

  list(
    mean_effect = fallback_result$mean_effect,
    mean_p = fallback_result$mean_p,
    disp_effect = fallback_result$disp_effect,
    disp_p = fallback_result$disp_p,
    convergence_flag = fallback_result$convergence_flag,
    model_family = fallback_result$family,
    fit_engine = fallback_result$engine,
    failure_reason = if (length(combined_failure) == 0L) NA_character_ else paste(combined_failure, collapse = " | "),
    warning_text = fallback_result$warning_text
  )
}

pick_example_genes <- function(results_df, per_category) {
  categories <- c("mean_only", "dispersion_only", "both", "neither")
  selected <- character(0)
  for (category in categories) {
    subset_df <- results_df[results_df$classification == category, , drop = FALSE]
    if (nrow(subset_df) == 0L) {
      next
    }
    ranking_score <- pmin(
      ifelse(is.na(subset_df$mean_FDR), Inf, subset_df$mean_FDR),
      ifelse(is.na(subset_df$disp_FDR), Inf, subset_df$disp_FDR)
    )
    subset_df$ranking_score <- ranking_score
    subset_df <- subset_df[order(subset_df$ranking_score, subset_df$Gene), , drop = FALSE]
    selected <- c(selected, utils::head(subset_df$Gene_internal, per_category))
  }
  unique(selected)
}

write_filtered_counts <- function(count_matrix, gene_names, output_file) {
  filtered_counts_df <- data.frame(
    Gene = gene_names,
    count_matrix,
    check.names = FALSE
  )
  utils::write.csv(filtered_counts_df, file = output_file, row.names = FALSE, quote = TRUE)
}

write_qc_readme <- function(output_file, context) {
  lines <- c(
    "# Per-gene mean and dispersion analysis",
    "",
    "## Inputs",
    paste0("- Counts file: `", context$counts_file, "`"),
    paste0("- Sample key file: `", context$sample_key_file, "`"),
    "",
    "## Sample alignment",
    paste0("- Samples in count matrix: ", context$n_count_samples),
    paste0("- Samples in sample key: ", context$n_metadata_samples),
    paste0("- Samples retained after intersection: ", context$n_intersection_samples),
    paste0("- Samples excluded because they were absent from the sample key: ", context$n_missing_in_key),
    paste0("- Samples excluded because they were absent from the count matrix: ", context$n_missing_in_counts),
    "",
    "## Sample QC filtering",
    paste0(
      "- Low-library threshold: samples with library size < ",
      format(round(context$libsize_threshold_raw, 2), trim = TRUE),
      " (computed as median(log10(libsize + 1)) - ",
      context$libsize_mad_cutoff,
      " * MAD)"
    ),
    paste0(
      "- Low-detection threshold: samples with detected genes < ",
      format(round(context$detected_threshold_raw, 2), trim = TRUE),
      " (computed as median(log10(detected genes + 1)) - ",
      context$detected_mad_cutoff,
      " * MAD)"
    ),
    paste0("- Minimum library size hard filter: ", context$min_library_size),
    paste0("- Minimum detected genes hard filter: ", context$min_detected_genes),
    paste0("- Samples retained after QC: ", context$n_samples_after_qc),
    paste0("- Samples removed by QC: ", context$n_samples_removed),
    "",
    "## Gene filtering",
    paste0("- Minimum nonzero samples: ", context$min_nonzero_samples),
    paste0("- Minimum total count across retained worms: ", context$min_total_count),
    paste0("- Genes before filtering: ", context$n_genes_before),
    paste0("- Genes retained for modeling: ", context$n_genes_after),
    "",
    "## Models",
    paste0("- Mean model: `", deparse(context$mean_formula), "`"),
    paste0("- Dispersion model: `", deparse(context$disp_formula), "`"),
    paste0("- Preferred GAMLSS families: ", paste(context$preferred_families, collapse = ", ")),
    paste0("- Fallback model: `", context$fallback_engine, "`"),
    "",
    "## Interpretation",
    "- `mean_effect` is the daf-2 vs N2 genotype coefficient in the mean model on the model link scale.",
    "- `disp_effect` is the daf-2 vs N2 genotype coefficient in the dispersion model. Positive values indicate higher residual dispersion in daf-2 after accounting for mean expression, library size, and replicate.",
    "- `mean_FDR` and `disp_FDR` are Benjamini-Hochberg adjusted p-values across genes for the respective genotype tests.",
    "",
    "## Diagnostics",
    paste0("- Genes fit successfully without warnings: ", context$n_fit_success),
    paste0("- Genes fit with warnings or partial fallbacks: ", context$n_fit_warning),
    paste0("- Genes that failed: ", context$n_fit_failed),
    paste0("- Preferred GAMLSS package available: ", context$gamlss_available),
    "",
    "## Caveats",
    "- Per-gene count models can be slow for large gene sets. Use `--max-genes` for a smoke test before a full run.",
    "- If GAMLSS is unavailable or unstable, the script falls back to a two-step approximation: a negative-binomial mean model (`glm.nb`) followed by a gamma model on squared Pearson residuals for dispersion. That fallback is useful diagnostically but is not identical to a full joint GAMLSS likelihood fit.",
    "- Extremely sparse genes can still fail even after filtering; those failures are preserved in the output tables instead of being dropped silently."
  )
  writeLines(lines, con = output_file)
}

save_qc_distribution_plot <- function(sample_qc, output_file, lib_thr, det_thr) {
  grDevices::png(output_file, width = 1800, height = 1400, res = 180)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    grDevices::dev.off()
  }, add = TRUE)
  par(mfrow = c(2, 2), mar = c(5, 5, 4, 2) + 0.1)

  hist(log10(sample_qc$library_size + 1), breaks = 30, col = "grey80",
       main = "Library size", xlab = "log10(library size + 1)")
  abline(v = log10(lib_thr + 1), col = "firebrick", lwd = 2, lty = 2)

  hist(log10(sample_qc$detected_genes + 1), breaks = 30, col = "grey80",
       main = "Detected genes", xlab = "log10(detected genes + 1)")
  abline(v = log10(det_thr + 1), col = "firebrick", lwd = 2, lty = 2)

  plot(log10(sample_qc$library_size + 1), log10(sample_qc$detected_genes + 1),
       pch = 16, cex = 0.8,
       col = ifelse(sample_qc$keep_sample, "#1B9E7790", "#D95F0290"),
       xlab = "log10(library size + 1)",
       ylab = "log10(detected genes + 1)",
       main = "QC scatter")
  abline(v = log10(lib_thr + 1), col = "firebrick", lwd = 2, lty = 2)
  abline(h = log10(det_thr + 1), col = "firebrick", lwd = 2, lty = 2)

  barplot(table(sample_qc$keep_sample), col = c("firebrick", "steelblue"),
          names.arg = c("Removed", "Retained"),
          main = "QC decision counts", ylab = "Samples")
}

save_qc_group_plot <- function(sample_qc, output_file) {
  sample_qc$group <- paste(sample_qc$Strain, sample_qc$Replicate, sep = "_")
  group_levels <- unique(sample_qc$group)
  group_colors <- setNames(
    rep(c("#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E7298A", "#E6AB02"), length.out = length(group_levels)),
    group_levels
  )

  grDevices::png(output_file, width = 2200, height = 1400, res = 180)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    grDevices::dev.off()
  }, add = TRUE)
  par(mfrow = c(1, 2), mar = c(9, 5, 4, 2) + 0.1)

  boxplot(log10(library_size + 1) ~ group, data = sample_qc,
          las = 2, col = group_colors[group_levels],
          main = "Library size by strain/replicate",
          ylab = "log10(library size + 1)")

  boxplot(log10(detected_genes + 1) ~ group, data = sample_qc,
          las = 2, col = group_colors[group_levels],
          main = "Detected genes by strain/replicate",
          ylab = "log10(detected genes + 1)")
}

save_volcano_plot <- function(results_df, effect_col, p_col, title, output_file) {
  p_values <- results_df[[p_col]]
  y_values <- -log10(p_values)
  y_values[!is.finite(y_values)] <- NA_real_
  sig_flag <- !is.na(results_df[[sub("_p$", "_FDR", p_col)]]) &
    results_df[[sub("_p$", "_FDR", p_col)]] < 0.05

  grDevices::png(output_file, width = 1600, height = 1400, res = 180)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    grDevices::dev.off()
  }, add = TRUE)
  par(mar = c(5, 5, 4, 2) + 0.1)

  plot(results_df[[effect_col]], y_values,
       pch = 16, cex = 0.7,
       col = ifelse(sig_flag, "#D95F0290", "#4D4D4D60"),
       xlab = effect_col,
       ylab = "-log10(p-value)",
       main = title)
  abline(h = -log10(0.05), col = "steelblue", lty = 2, lwd = 2)
  abline(v = 0, col = "grey60", lty = 3)
}

save_effect_scatter <- function(results_df, output_file) {
  grDevices::png(output_file, width = 1600, height = 1400, res = 180)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    grDevices::dev.off()
  }, add = TRUE)
  par(mar = c(5, 5, 4, 2) + 0.1)

  class_cols <- c(
    mean_only = "#1B9E77A0",
    dispersion_only = "#D95F02A0",
    both = "#7570B3A0",
    neither = "#4D4D4D50"
  )
  plot(results_df$mean_effect, results_df$disp_effect,
       pch = 16, cex = 0.7,
       col = class_cols[results_df$classification],
       xlab = "Mean effect",
       ylab = "Dispersion effect",
       main = "Mean vs dispersion effect")
  abline(h = 0, v = 0, col = "grey60", lty = 3)
  legend("topright", legend = names(class_cols), col = class_cols, pch = 16, bty = "n")
}

save_mean_variance_plot <- function(summary_stats, strain_levels, output_file) {
  grDevices::png(output_file, width = 1800, height = 900, res = 180)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    grDevices::dev.off()
  }, add = TRUE)
  par(mfrow = c(1, length(strain_levels)), mar = c(5, 5, 4, 2) + 0.1)

  for (strain in strain_levels) {
    mean_col <- paste0("mean_count_", sanitize_name(strain))
    var_col <- paste0("variance_count_", sanitize_name(strain))
    x <- log10(summary_stats[[mean_col]] + 1e-6)
    y <- log10(summary_stats[[var_col]] + 1e-6)
    plot(x, y, pch = 16, cex = 0.45, col = "#1F78B480",
         xlab = paste0("log10 mean count (", strain, ")"),
         ylab = paste0("log10 variance (", strain, ")"),
         main = paste("Mean-variance:", strain))
    smooth_ok <- is.finite(x) & is.finite(y)
    if (sum(smooth_ok) > 10L) {
      smooth_line <- stats::lowess(x[smooth_ok], y[smooth_ok], f = 0.3)
      lines(smooth_line, col = "firebrick", lwd = 2)
    }
  }
}

save_example_gene_plots <- function(example_gene_ids, count_matrix, gene_ids_internal, gene_labels, sample_df, results_df, output_file) {
  if (length(example_gene_ids) == 0L) {
    return(invisible(NULL))
  }

  gene_index <- match(example_gene_ids, gene_ids_internal)
  keep <- !is.na(gene_index)
  example_gene_ids <- example_gene_ids[keep]
  gene_index <- gene_index[keep]
  if (length(example_gene_ids) == 0L) {
    return(invisible(NULL))
  }

  replicate_levels <- levels(sample_df$Replicate)
  replicate_cols <- setNames(
    rep(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), length.out = length(replicate_levels)),
    replicate_levels
  )
  strain_positions <- seq_along(levels(sample_df$Strain))
  strain_map <- setNames(strain_positions, levels(sample_df$Strain))

  grDevices::pdf(output_file, width = 10, height = 8)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    grDevices::dev.off()
  }, add = TRUE)

  per_page <- 4L
  for (page_start in seq(1L, length(example_gene_ids), by = per_page)) {
    page_end <- min(page_start + per_page - 1L, length(example_gene_ids))
    page_idx <- page_start:page_end
    par(mfrow = c(2, 2), mar = c(5, 5, 4, 2) + 0.1)

    for (i in page_idx) {
      gene_internal <- example_gene_ids[[i]]
      gene <- gene_labels[[gene_index[[i]]]]
      counts <- as.numeric(count_matrix[gene_index[[i]], ])
      x <- unname(strain_map[as.character(sample_df$Strain)]) +
        stats::runif(length(counts), min = -0.16, max = 0.16)
      y <- log10(counts + 1)

      plot(x, y, xaxt = "n", pch = 16,
           col = replicate_cols[as.character(sample_df$Replicate)],
           xlab = "Strain", ylab = "log10(count + 1)",
           main = gene)
      axis(1, at = strain_positions, labels = levels(sample_df$Strain))
      boxplot(log10(counts + 1) ~ sample_df$Strain, add = TRUE, axes = FALSE,
              border = "grey30", col = NA)

      result_row <- results_df[results_df$Gene_internal == gene_internal, , drop = FALSE]
      subtitle_text <- paste0(
        "class=", result_row$classification[[1]],
        " | mean_FDR=", signif(result_row$mean_FDR[[1]], 3),
        " | disp_FDR=", signif(result_row$disp_FDR[[1]], 3)
      )
      mtext(subtitle_text, side = 3, line = 0.2, cex = 0.8)
    }

    if (length(page_idx) < per_page) {
      for (unused in seq_len(per_page - length(page_idx))) {
        plot.new()
      }
    }
  }
}

# -------------------------------- Main analysis ------------------------------
cli_args <- parse_args(commandArgs(trailingOnly = TRUE))
config <- modifyList(default_config, user_config)

if (!is.null(cli_args$config)) {
  config <- modifyList(config, load_config_overrides(cli_args$config))
}
if (!is.null(cli_args$counts)) {
  config$counts_file <- cli_args$counts
}
if (!is.null(cli_args$`sample-key`)) {
  config$sample_key_file <- cli_args$`sample-key`
}
if (!is.null(cli_args$`output-dir`)) {
  config$output_dir <- cli_args$`output-dir`
}
if (!is.null(cli_args$`max-genes`)) {
  config$debug$max_genes <- as.numeric(cli_args$`max-genes`)
}

stop_if_missing_package("MASS", "used for the negative-binomial fallback model.")

output_dir <- normalizePath(config$output_dir, winslash = "/", mustWork = FALSE)
tables_dir <- file.path(output_dir, "tables")
plots_dir <- file.path(output_dir, "plots")
logs_dir <- file.path(output_dir, "logs")
safe_dir_create(output_dir)
safe_dir_create(tables_dir)
safe_dir_create(plots_dir)
safe_dir_create(logs_dir)

message("Reading inputs...")
counts_df <- utils::read.csv(
  file = config$counts_file,
  check.names = FALSE
)
sample_key <- utils::read.csv(
  file = config$sample_key_file,
  check.names = FALSE
)

required_sample_cols <- c(
  config$sample_column,
  config$strain_column,
  config$replicate_column
)
missing_sample_cols <- setdiff(required_sample_cols, names(sample_key))
if (length(missing_sample_cols) > 0L) {
  stop("Sample key is missing required columns: ", paste(missing_sample_cols, collapse = ", "), call. = FALSE)
}
if (!config$gene_column %in% names(counts_df)) {
  stop("Counts file is missing the gene column `", config$gene_column, "`.", call. = FALSE)
}

gene_ids_original <- as.character(counts_df[[config$gene_column]])
gene_ids_internal <- make.unique(gene_ids_original)
count_sample_names <- setdiff(names(counts_df), config$gene_column)

count_matrix <- as.matrix(counts_df[, count_sample_names, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
rownames(count_matrix) <- gene_ids_internal

metadata_sample_names <- as.character(sample_key[[config$sample_column]])
samples_in_both <- intersect(count_sample_names, metadata_sample_names)
missing_in_key <- setdiff(count_sample_names, metadata_sample_names)
missing_in_counts <- setdiff(metadata_sample_names, count_sample_names)

if (length(samples_in_both) == 0L) {
  stop("No overlapping samples were found between the count matrix and sample key.", call. = FALSE)
}

sample_key_aligned <- sample_key[match(samples_in_both, metadata_sample_names), , drop = FALSE]
count_matrix_aligned <- count_matrix[, samples_in_both, drop = FALSE]

if (!identical(colnames(count_matrix_aligned), as.character(sample_key_aligned[[config$sample_column]]))) {
  stop("Sample order mismatch remains after alignment.", call. = FALSE)
}

sample_key_aligned[[config$strain_column]] <- factor(sample_key_aligned[[config$strain_column]])
if (!config$strain_reference %in% levels(sample_key_aligned[[config$strain_column]])) {
  stop("Strain reference `", config$strain_reference, "` is not present after sample alignment.", call. = FALSE)
}
sample_key_aligned[[config$strain_column]] <- stats::relevel(
  sample_key_aligned[[config$strain_column]],
  ref = config$strain_reference
)
sample_key_aligned[[config$replicate_column]] <- factor(sample_key_aligned[[config$replicate_column]])

sample_key_aligned$library_size <- colSums(count_matrix_aligned)
sample_key_aligned$detected_genes <- colSums(count_matrix_aligned > 0)
sample_key_aligned$log_library_size <- log(pmax(sample_key_aligned$library_size, 1))

message("Computing sample QC thresholds...")
libsize_threshold <- mad_lower_threshold(
  sample_key_aligned$library_size,
  nmads = config$qc$libsize_mad_cutoff
)
detected_threshold <- mad_lower_threshold(
  sample_key_aligned$detected_genes,
  nmads = config$qc$detected_mad_cutoff
)

sample_key_aligned$exclude_low_library <- sample_key_aligned$library_size < max(
  config$qc$min_library_size,
  libsize_threshold$threshold_raw
)
sample_key_aligned$exclude_low_detected <- sample_key_aligned$detected_genes < max(
  config$qc$min_detected_genes,
  detected_threshold$threshold_raw
)
sample_key_aligned$exclude_missing_design <- is.na(sample_key_aligned[[config$strain_column]]) |
  is.na(sample_key_aligned[[config$replicate_column]])
sample_key_aligned$keep_sample <- !(
  sample_key_aligned$exclude_low_library |
    sample_key_aligned$exclude_low_detected |
    sample_key_aligned$exclude_missing_design
)
sample_key_aligned$exclusion_reason <- vapply(
  seq_len(nrow(sample_key_aligned)),
  function(i) {
    collapse_reasons(
      if (sample_key_aligned$exclude_low_library[[i]]) "low_library_size" else NULL,
      if (sample_key_aligned$exclude_low_detected[[i]]) "low_detected_genes" else NULL,
      if (sample_key_aligned$exclude_missing_design[[i]]) "missing_design_metadata" else NULL
    )
  },
  FUN.VALUE = character(1)
)

sample_qc_table <- sample_key_aligned
sample_qc_table$alignment_status <- ifelse(
  sample_qc_table[[config$sample_column]] %in% samples_in_both,
  "matched",
  "unmatched"
)

excluded_samples <- sample_qc_table[!sample_qc_table$keep_sample, , drop = FALSE]
sample_filtering_summary <- as.data.frame(
  table(
    Strain = sample_qc_table[[config$strain_column]],
    Replicate = sample_qc_table[[config$replicate_column]],
    Keep = sample_qc_table$keep_sample
  )
)

retained_sample_mask <- sample_qc_table$keep_sample
sample_metadata_filtered <- droplevels(sample_qc_table[retained_sample_mask, , drop = FALSE])
count_matrix_filtered_samples <- count_matrix_aligned[, retained_sample_mask, drop = FALSE]

if (nrow(sample_metadata_filtered) < 4L) {
  stop("Too few samples remain after QC filtering to fit the requested models.", call. = FALSE)
}
if (length(unique(sample_metadata_filtered[[config$strain_column]])) < 2L) {
  stop("Only one strain remains after QC filtering; the genotype test is no longer identifiable.", call. = FALSE)
}

message("Filtering genes...")
min_nonzero_samples <- max(
  config$gene_filter$min_nonzero_samples,
  ceiling(config$gene_filter$min_nonzero_fraction * ncol(count_matrix_filtered_samples))
)
nonzero_counts <- rowSums(count_matrix_filtered_samples > 0)
total_counts <- rowSums(count_matrix_filtered_samples)
gene_keep <- nonzero_counts >= min_nonzero_samples &
  total_counts >= config$gene_filter$min_total_count

count_matrix_gene_filtered <- count_matrix_filtered_samples[gene_keep, , drop = FALSE]
gene_ids_model <- gene_ids_original[gene_keep]
gene_ids_model_internal <- rownames(count_matrix_filtered_samples)[gene_keep]
if (is.finite(config$debug$max_genes)) {
  keep_n <- min(nrow(count_matrix_gene_filtered), as.integer(config$debug$max_genes))
  count_matrix_gene_filtered <- count_matrix_gene_filtered[seq_len(keep_n), , drop = FALSE]
  gene_ids_model <- gene_ids_model[seq_len(keep_n)]
  gene_ids_model_internal <- gene_ids_model_internal[seq_len(keep_n)]
}

if (nrow(count_matrix_gene_filtered) == 0L) {
  stop("No genes passed the post-QC filtering criteria.", call. = FALSE)
}

mean_formula <- build_model_formula(
  response = "count_response",
  terms = config$model$mean_terms,
  include_offset = isTRUE(config$model$include_library_offset),
  available_data = sample_metadata_filtered
)
disp_formula <- build_model_formula(
  response = "count_response",
  terms = config$model$dispersion_terms,
  include_offset = FALSE,
  available_data = sample_metadata_filtered
)

strain_levels <- levels(sample_metadata_filtered[[config$strain_column]])
strain_term <- find_strain_term(
  coef_names = colnames(stats::model.matrix(rhs_only_formula(mean_formula), data = sample_metadata_filtered)),
  strain_column = config$strain_column,
  strain_levels = strain_levels,
  reference_level = config$strain_reference
)
if (is.na(strain_term)) {
  stop("Unable to identify the strain coefficient for daf-2 vs N2.", call. = FALSE)
}

message("Computing per-gene summary statistics...")
summary_stats <- data.frame(
  Gene_internal = gene_ids_model_internal,
  Gene = gene_ids_model,
  stringsAsFactors = FALSE
)
for (strain in strain_levels) {
  idx <- sample_metadata_filtered[[config$strain_column]] == strain
  strain_counts <- count_matrix_gene_filtered[, idx, drop = FALSE]
  suffix <- sanitize_name(strain)
  summary_stats[[paste0("mean_count_", suffix)]] <- rowMeans(strain_counts)
  summary_stats[[paste0("variance_count_", suffix)]] <- apply(strain_counts, 1, safe_variance)
  summary_stats[[paste0("prop_zero_", suffix)]] <- rowMeans(strain_counts == 0)
}

message("Fitting per-gene models...")
results_list <- vector("list", nrow(count_matrix_gene_filtered))
progress <- utils::txtProgressBar(min = 0, max = nrow(count_matrix_gene_filtered), style = 3)
for (i in seq_len(nrow(count_matrix_gene_filtered))) {
  gene_fit <- fit_single_gene(
    y = count_matrix_gene_filtered[i, ],
    sample_df = sample_metadata_filtered,
    mean_formula = mean_formula,
    sigma_formula = disp_formula,
    strain_term = strain_term,
    config = config
  )
  results_list[[i]] <- c(
    list(
      Gene_internal = gene_ids_model_internal[[i]],
      Gene = gene_ids_model[[i]],
      n_samples_used = nrow(sample_metadata_filtered)
    ),
    gene_fit
  )
  utils::setTxtProgressBar(progress, i)
}
close(progress)

results_df <- do.call(rbind, lapply(results_list, as.data.frame, stringsAsFactors = FALSE))
rownames(results_df) <- NULL

results_df$mean_FDR <- stats::p.adjust(results_df$mean_p, method = "BH")
results_df$disp_FDR <- stats::p.adjust(results_df$disp_p, method = "BH")
results_df$classification <- ifelse(
  !is.na(results_df$mean_FDR) & results_df$mean_FDR < config$significance_fdr &
    !is.na(results_df$disp_FDR) & results_df$disp_FDR < config$significance_fdr,
  "both",
  ifelse(
    !is.na(results_df$mean_FDR) & results_df$mean_FDR < config$significance_fdr,
    "mean_only",
    ifelse(
      !is.na(results_df$disp_FDR) & results_df$disp_FDR < config$significance_fdr,
      "dispersion_only",
      "neither"
    )
  )
)

results_df <- merge(results_df, summary_stats, by = c("Gene_internal", "Gene"), all.x = TRUE, sort = FALSE)
results_df <- results_df[match(gene_ids_model_internal, results_df$Gene_internal), , drop = FALSE]

results_df$failure_reason_simple <- vapply(
  seq_len(nrow(results_df)),
  function(i) {
    reason <- results_df$failure_reason[[i]]
    if (is.na(reason) || !nzchar(reason)) {
      "none"
    } else {
      sub(":.*$", "", strsplit(reason, " \\| ", perl = TRUE)[[1]][1])
    }
  },
  FUN.VALUE = character(1)
)

diagnostic_summary <- data.frame(
  metric = c(
    "genes_attempted",
    "genes_fit_success",
    "genes_fit_warning",
    "genes_fit_failed",
    "gamlss_available"
  ),
  value = c(
    nrow(results_df),
    sum(results_df$convergence_flag == "success"),
    sum(results_df$convergence_flag %in% c("warning", "partial_failure")),
    sum(results_df$convergence_flag == "failed"),
    requireNamespace("gamlss", quietly = TRUE) && requireNamespace("gamlss.dist", quietly = TRUE)
  ),
  stringsAsFactors = FALSE
)

failure_reason_counts <- sort(table(results_df$failure_reason_simple), decreasing = TRUE)
failure_reason_table <- data.frame(
  failure_reason = names(failure_reason_counts),
  n_genes = as.integer(failure_reason_counts),
  stringsAsFactors = FALSE
)

message("Writing output tables...")
utils::write.csv(sample_metadata_filtered, file.path(tables_dir, "filtered_sample_metadata.csv"), row.names = FALSE)
write_filtered_counts(
  count_matrix = count_matrix_gene_filtered,
  gene_names = gene_ids_model,
  output_file = file.path(tables_dir, "filtered_counts_matrix.csv")
)
utils::write.csv(sample_qc_table, file.path(tables_dir, "qc_summary_table.csv"), row.names = FALSE)
utils::write.csv(excluded_samples, file.path(tables_dir, "excluded_samples.csv"), row.names = FALSE)
utils::write.csv(sample_filtering_summary, file.path(tables_dir, "sample_filtering_summary.csv"), row.names = FALSE)
utils::write.csv(results_df, file.path(tables_dir, "per_gene_results.csv"), row.names = FALSE)
utils::write.csv(diagnostic_summary, file.path(tables_dir, "diagnostic_summary.csv"), row.names = FALSE)
utils::write.csv(failure_reason_table, file.path(tables_dir, "failure_reason_summary.csv"), row.names = FALSE)
alignment_exclusions <- rbind(
  data.frame(
    source = rep("counts_only", length(missing_in_key)),
    sample_name = missing_in_key,
    stringsAsFactors = FALSE
  ),
  data.frame(
    source = rep("sample_key_only", length(missing_in_counts)),
    sample_name = missing_in_counts,
    stringsAsFactors = FALSE
  )
)
utils::write.csv(alignment_exclusions, file.path(tables_dir, "alignment_exclusions.csv"), row.names = FALSE)

message("Saving plots...")
save_qc_distribution_plot(
  sample_qc = sample_qc_table,
  output_file = file.path(plots_dir, "sample_qc_distributions.png"),
  lib_thr = max(config$qc$min_library_size, libsize_threshold$threshold_raw),
  det_thr = max(config$qc$min_detected_genes, detected_threshold$threshold_raw)
)
save_qc_group_plot(
  sample_qc = sample_qc_table,
  output_file = file.path(plots_dir, "sample_qc_by_strain_replicate.png")
)
save_volcano_plot(
  results_df = results_df,
  effect_col = "mean_effect",
  p_col = "mean_p",
  title = "Volcano plot: mean effect",
  output_file = file.path(plots_dir, "volcano_mean.png")
)
save_volcano_plot(
  results_df = results_df,
  effect_col = "disp_effect",
  p_col = "disp_p",
  title = "Volcano plot: dispersion effect",
  output_file = file.path(plots_dir, "volcano_dispersion.png")
)
save_effect_scatter(
  results_df = results_df,
  output_file = file.path(plots_dir, "mean_vs_dispersion_effect.png")
)
save_mean_variance_plot(
  summary_stats = results_df,
  strain_levels = strain_levels,
  output_file = file.path(plots_dir, "mean_variance_by_genotype.png")
)
save_example_gene_plots(
  example_gene_ids = pick_example_genes(results_df, config$plotting$example_genes_per_category),
  count_matrix = count_matrix_gene_filtered,
  gene_ids_internal = gene_ids_model_internal,
  gene_labels = gene_ids_model,
  sample_df = sample_metadata_filtered,
  results_df = results_df,
  output_file = file.path(plots_dir, "example_genes.pdf")
)

write_qc_readme(
  output_file = file.path(output_dir, "README.md"),
  context = list(
    counts_file = config$counts_file,
    sample_key_file = config$sample_key_file,
    n_count_samples = length(count_sample_names),
    n_metadata_samples = nrow(sample_key),
    n_intersection_samples = length(samples_in_both),
    n_missing_in_key = length(missing_in_key),
    n_missing_in_counts = length(missing_in_counts),
    libsize_threshold_raw = max(config$qc$min_library_size, libsize_threshold$threshold_raw),
    detected_threshold_raw = max(config$qc$min_detected_genes, detected_threshold$threshold_raw),
    libsize_mad_cutoff = config$qc$libsize_mad_cutoff,
    detected_mad_cutoff = config$qc$detected_mad_cutoff,
    min_library_size = config$qc$min_library_size,
    min_detected_genes = config$qc$min_detected_genes,
    n_samples_after_qc = nrow(sample_metadata_filtered),
    n_samples_removed = sum(!sample_qc_table$keep_sample),
    min_nonzero_samples = min_nonzero_samples,
    min_total_count = config$gene_filter$min_total_count,
    n_genes_before = nrow(count_matrix_filtered_samples),
    n_genes_after = nrow(count_matrix_gene_filtered),
    mean_formula = mean_formula,
    disp_formula = disp_formula,
    preferred_families = config$model$preferred_families,
    fallback_engine = config$model$fallback_engine,
    n_fit_success = sum(results_df$convergence_flag == "success"),
    n_fit_warning = sum(results_df$convergence_flag %in% c("warning", "partial_failure")),
    n_fit_failed = sum(results_df$convergence_flag == "failed"),
    gamlss_available = requireNamespace("gamlss", quietly = TRUE) && requireNamespace("gamlss.dist", quietly = TRUE)
  )
)

utils::capture.output(sessionInfo(), file = file.path(logs_dir, "session_info.txt"))

message("Analysis complete.")

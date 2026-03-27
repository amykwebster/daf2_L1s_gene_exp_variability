# ===============================
# PCA on filtered samples
# ===============================

# 0) Setup ===============================
library(tidyverse)
library(edgeR)
library(ggrepel)

setwd("/Users/sophi/Documents/CELseq_tutorial/pipeline/")
counts_path <- "umi_count_matrix_all_combined_3sets.csv"
annot_path  <- "WS273_geneNames.csv"

n_top_genes     <- 1000
expr_thresh_cpm <- 1
expressed_frac  <- 0.60
min_reads       <- 200000

# 1) Load counts ===============================
CountFile <- read.csv(counts_path, header = TRUE, check.names = TRUE)
stopifnot("Gene" %in% names(CountFile))

# 2) Parse set & sample ===============================
all_samp <- setdiff(names(CountFile), "Gene")
clean    <- sub("^X", "", all_samp)

spl      <- strsplit(clean, "_")
set_vec  <- vapply(spl, function(x) paste(x[-length(x)], collapse = "_"), character(1))
samp_vec <- vapply(spl, function(x) x[length(x)], character(1))
samp_num <- suppressWarnings(as.numeric(sub("S$", "", samp_vec)))

sample_info <- tibble(
  orig = all_samp,
  clean = clean,
  set = set_vec,
  sample = samp_vec,
  sample_num = samp_num
)

ord <- order(sample_info$set, sample_info$sample_num)
sample_info <- sample_info[ord, ]
col_order   <- sample_info$orig

# 3) Build counts matrix===============================
CountFile <- CountFile %>% dplyr::select(Gene, all_of(col_order))

counts <- as.matrix(CountFile[, -1, drop = FALSE])
rownames(counts) <- CountFile$Gene
colnames(counts) <- sub("^X", "", colnames(counts))
row.names(sample_info) <- sample_info$clean

sets <- sample_info$set
names(sets) <- sample_info$clean

# 4) Create functions ===============================
logcpm_from_counts <- function(counts_mat) {
  y <- DGEList(counts = counts_mat)
  y <- calcNormFactors(y, method = "TMM")
  cpm(y, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
}

plot_pca <- function(logCPM, sets, title_prefix, n_top = 1000, out_prefix = "PCA_plot") {
  if (is.finite(n_top)) {
    gvar <- apply(logCPM, 1, var, na.rm = TRUE)
    top_idx <- order(gvar, decreasing = TRUE)[seq_len(min(n_top, nrow(logCPM)))]
    logCPM <- logCPM[top_idx, , drop = FALSE]
  }
  
  pca <- prcomp(t(logCPM), center = TRUE, scale. = TRUE)
  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
  
  df <- data.frame(
    Sample = colnames(logCPM),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    Set = sets[colnames(logCPM)],
    stringsAsFactors = FALSE
  )
  
  plot_one <- function(pc_x, pc_y) {
    x_col <- paste0("PC", pc_x)
    y_col <- paste0("PC", pc_y)
    x_lab <- paste0("PC", pc_x, " (", round(100 * var_exp[pc_x], 1), "%)")
    y_lab <- paste0("PC", pc_y, " (", round(100 * var_exp[pc_y], 1), "%)")
    
    p <- ggplot(df, aes_string(x = x_col, y = y_col, color = "Set", label = "Sample")) +
      geom_point(size = 2.5) +
      labs(
        title = paste0(title_prefix, " (", x_col, " vs ", y_col, ")"),
        x = x_lab,
        y = y_lab
      ) +
      theme_classic(base_size = 11) +
      guides(color = guide_legend(override.aes = list(size = 4)))
    
    print(p)
    
    fn <- sprintf("%s_%s_vs_%s.png", out_prefix, x_col, y_col)
    ggsave(fn, p, width = 8, height = 6, dpi = 300)
  }
  
  plot_one(1, 2)
  plot_one(2, 3)
  plot_one(1, 3)
}

# 5) Apply filters to counts ===============================

#a) filter > 200,000
totals <- colSums(counts, na.rm = TRUE)
keep_samples <- totals > min_reads
counts_filt <- counts[, keep_samples, drop = FALSE]

#b) filter for protein-coding genes
annot <- read.csv(annot_path, header = TRUE)
stopifnot(all(c("WB_id", "type", "live_dead") %in% names(annot)))

pc_ids <- annot$WB_id[annot$type == "protein_coding_gene" & annot$live_dead == "Live"]
counts_filt <- counts_filt[rownames(counts_filt) %in% pc_ids, , drop = FALSE]

#c) filter for expressed genes
y_tmp <- DGEList(counts = counts_filt)
cpm_tmp <- cpm(y_tmp)

min_expr <- ceiling(ncol(counts_filt) * expressed_frac)
keep_genes_expr <- rowSums(cpm_tmp > expr_thresh_cpm) >= min_expr

counts_filt <- counts_filt[keep_genes_expr, , drop = FALSE]

# 6) Generate PCA ===============================
sets_filt <- sets[colnames(counts_filt)]
logCPM_filt <- logcpm_from_counts(counts_filt)

plot_pca(
  logCPM_filt,
  sets = sets_filt[colnames(logCPM_filt)],
  title_prefix = sprintf(
    "PCA: >200k reads + protein-coding + expressed (CPM>%g in >=%d%% samples)",
    expr_thresh_cpm, round(100 * expressed_frac)
  ),
  n_top = n_top_genes,
  out_prefix = "PCA_filtered_gt200k_proteinCoding_expressed_by_set"
)

# 7) Summary ===============================
cat("Final filtered PCA genes:", nrow(logCPM_filt), " samples:", ncol(logCPM_filt), "\n")


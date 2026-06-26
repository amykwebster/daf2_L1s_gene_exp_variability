#start analysis with variance-stabilized residuals (R_units)
get_pairwise_gene_cors <- function(df, min_worm_frac = 0.8, method = "pearson") {
  
  df2 <- df %>%
    dplyr::select(Strain, SampleIDs, gene, R_units) %>%
    filter(!is.na(R_units))
  
  df2 %>%
    group_split(Strain) %>%
    map_dfr(function(strain_df) {
      
      strain_name <- unique(strain_df$Strain)
      n_worms <- n_distinct(strain_df$SampleIDs)
      
      keep_genes <- strain_df %>%
        group_by(gene) %>%
        summarize(
          n_worms_measured = n_distinct(SampleIDs),
          .groups = "drop"
        ) %>%
        filter(n_worms_measured >= min_worm_frac * n_worms) %>%
        pull(gene)
      
      wide_df <- strain_df %>%
        filter(gene %in% keep_genes) %>%
        pivot_wider(
          names_from = gene,
          values_from = R_units
        )
      
      mat <- wide_df %>%
        dplyr::select(-Strain, -SampleIDs) %>%
        as.matrix()
      
      cor_mat <- cor(
        mat,
        use = "pairwise.complete.obs",
        method = method
      )
      
      pairwise_cors <- cor_mat[upper.tri(cor_mat)]
      
      tibble(
        Strain = strain_name,
        n_worms = n_worms,
        n_genes = ncol(mat),
        r = pairwise_cors
      )
    })
}

pairwise_gene_cor_by_strain <- get_pairwise_gene_cors(
  GroupBySample2,
  min_worm_frac = 0.8,
  method = "pearson"
)

cor_summary_by_strain <- pairwise_gene_cor_by_strain %>%
  group_by(Strain) %>%
  summarize(
    n_worms = dplyr::first(n_worms),
    n_genes = dplyr::first(n_genes),
    n_gene_pairs = n(),
    mean_r = mean(r, na.rm = TRUE),
    median_r = median(r, na.rm = TRUE),
    prop_positive = mean(r > 0, na.rm = TRUE),
    prop_gt_0.1 = mean(r > 0.1, na.rm = TRUE),
    prop_gt_0.2 = mean(r > 0.2, na.rm = TRUE),
    prop_lt_neg_0.1 = mean(r < -0.1, na.rm = TRUE),
    prop_lt_neg_0.2 = mean(r < -0.2, na.rm = TRUE),
    .groups = "drop"
  )

cor_summary_by_strain


ggplot(pairwise_gene_cor_by_strain, aes(x = r, fill = Strain)) +
  geom_density(alpha = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(
    x = "Pairwise gene-gene correlation across single worms",
    y = "Density",
    title = "Gene-set co-expression across individuals",
    subtitle = "Each value is the correlation between two genes across worms"
  )



library(tidyverse)

set.seed(1)

n_boot <- 2000
gene_set <- GeneLists$Tepper_ClassI

# Make wide matrices once from all genes in R_long
wide_all <- R_long %>%
  dplyr::select(Strain, SampleIDs, gene, R_units) %>%
  filter(!is.na(R_units)) %>%
  pivot_wider(names_from = gene, values_from = R_units)

mats <- wide_all %>%
  group_split(Strain) %>%
  set_names(map_chr(., ~ unique(.x$Strain))) %>%
  map(~ {
    mat <- .x %>%
      dplyr::select(-Strain, -SampleIDs) %>%
      as.matrix()
    rownames(mat) <- .x$SampleIDs
    mat
  })

# Only use genes present as matrix columns for all strains
common_genes <- Reduce(intersect, map(mats, colnames))
gene_set <- intersect(gene_set, common_genes)
background_genes <- setdiff(common_genes, gene_set)
n_genes <- length(gene_set)

get_cor_summary <- function(genes, label = "observed", iter = NA) {
  map_dfr(names(mats), function(strain_name) {
    mat <- mats[[strain_name]][, genes, drop = FALSE]
    cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
    r <- cor_mat[upper.tri(cor_mat)]
    
    tibble(
      label = label,
      iter = iter,
      Strain = strain_name,
      n_worms = nrow(mat),
      n_genes = ncol(mat),
      n_gene_pairs = length(r),
      mean_r = mean(r, na.rm = TRUE),
      median_r = median(r, na.rm = TRUE),
      prop_positive = mean(r > 0, na.rm = TRUE),
      prop_gt_0.1 = mean(r > 0.1, na.rm = TRUE),
      prop_gt_0.2 = mean(r > 0.2, na.rm = TRUE),
      prop_lt_neg_0.1 = mean(r < -0.1, na.rm = TRUE),
      prop_lt_neg_0.2 = mean(r < -0.2, na.rm = TRUE)
    )
  })
}

observed_summary <- get_cor_summary(gene_set, label = "observed", iter = 0)

random_summary <- map_dfr(seq_len(n_boot), function(i) {
  random_genes <- sample(background_genes, n_genes, replace = FALSE)
  get_cor_summary(random_genes, label = "random", iter = i)
})

all_summary <- bind_rows(observed_summary, random_summary)

summary_diff <- all_summary %>%
  dplyr::select(label, iter, Strain, mean_r, median_r, prop_positive, prop_gt_0.1, prop_gt_0.2) %>%
  pivot_wider(
    names_from = Strain,
    values_from = c(mean_r, median_r, prop_positive, prop_gt_0.1, prop_gt_0.2)
  ) %>%
  mutate(
    diff_mean_r = mean_r_N2 - mean_r_daf2,
    diff_median_r = median_r_N2 - median_r_daf2,
    diff_prop_positive = prop_positive_N2 - prop_positive_daf2,
    diff_prop_gt_0.1 = prop_gt_0.1_N2 - prop_gt_0.1_daf2,
    diff_prop_gt_0.2 = prop_gt_0.2_N2 - prop_gt_0.2_daf2
  )

observed_summary

summary_diff %>%
  filter(label == "observed")

summary_diff %>%
  filter(label == "random") %>%
  summarize(
    random_mean_diff_mean_r = mean(diff_mean_r, na.rm = TRUE),
    random_95_low = quantile(diff_mean_r, 0.025, na.rm = TRUE),
    random_95_high = quantile(diff_mean_r, 0.975, na.rm = TRUE),
    observed_diff_mean_r = summary_diff$diff_mean_r[summary_diff$label == "observed"],
    empirical_p_N2_minus_daf2 = mean(diff_mean_r >= summary_diff$diff_mean_r[summary_diff$label == "observed"], na.rm = TRUE)
  )




library(tidyverse)

set.seed(1)

n_boot <- 2000
gene_set <- GeneLists$Tepper_ClassI

wide_all <- R_long %>%
  dplyr::select(Strain, SampleIDs, gene, R_units) %>%
  filter(!is.na(R_units), gene %in% gene_set) %>%
  pivot_wider(names_from = gene, values_from = R_units)

mats <- wide_all %>%
  group_split(Strain) %>%
  set_names(map_chr(., ~ unique(.x$Strain))) %>%
  map(~ {
    mat <- .x %>%
      dplyr::select(-Strain, -SampleIDs) %>%
      as.matrix()
    rownames(mat) <- .x$SampleIDs
    mat
  })

common_genes <- Reduce(intersect, map(mats, colnames))
mats <- map(mats, ~ .x[, common_genes, drop = FALSE])

cor_summary <- function(mat) {
  cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
  r <- cor_mat[upper.tri(cor_mat)]
  tibble(
    mean_r = mean(r, na.rm = TRUE),
    median_r = median(r, na.rm = TRUE),
    prop_positive = mean(r > 0, na.rm = TRUE),
    prop_gt_0.1 = mean(r > 0.1, na.rm = TRUE),
    prop_gt_0.2 = mean(r > 0.2, na.rm = TRUE)
  )
}

observed <- map_dfr(names(mats), function(strain_name) {
  cor_summary(mats[[strain_name]]) %>%
    mutate(iter = 0, Strain = strain_name, type = "observed")
})

boot <- map_dfr(seq_len(n_boot), function(i) {
  map_dfr(names(mats), function(strain_name) {
    
    mat0 <- mats[[strain_name]]
    
    sampled_worms <- sample(rownames(mat0), nrow(mat0), replace = TRUE)
    sampled_genes <- sample(colnames(mat0), ncol(mat0), replace = TRUE)
    
    mat_boot <- mat0[sampled_worms, sampled_genes, drop = FALSE]
    
    cor_summary(mat_boot) %>%
      mutate(iter = i, Strain = strain_name, type = "bootstrap")
  })
})

boot_diff <- bind_rows(observed, boot) %>%
  pivot_wider(
    names_from = Strain,
    values_from = c(mean_r, median_r, prop_positive, prop_gt_0.1, prop_gt_0.2)
  ) %>%
  mutate(
    diff_mean_r = mean_r_N2 - mean_r_daf2,
    diff_median_r = median_r_N2 - median_r_daf2,
    diff_prop_positive = prop_positive_N2 - prop_positive_daf2,
    diff_prop_gt_0.1 = prop_gt_0.1_N2 - prop_gt_0.1_daf2,
    diff_prop_gt_0.2 = prop_gt_0.2_N2 - prop_gt_0.2_daf2
  )

boot_diff %>%
  filter(type == "bootstrap") %>%
  summarize(
    mean_diff_mean_r = mean(diff_mean_r, na.rm = TRUE),
    lower_95_mean_r = quantile(diff_mean_r, 0.025, na.rm = TRUE),
    upper_95_mean_r = quantile(diff_mean_r, 0.975, na.rm = TRUE),
    prop_diff_mean_r_gt_0 = mean(diff_mean_r > 0, na.rm = TRUE),
    
    mean_diff_prop_gt_0.2 = mean(diff_prop_gt_0.2, na.rm = TRUE),
    lower_95_prop_gt_0.2 = quantile(diff_prop_gt_0.2, 0.025, na.rm = TRUE),
    upper_95_prop_gt_0.2 = quantile(diff_prop_gt_0.2, 0.975, na.rm = TRUE),
    prop_diff_prop_gt_0.2_gt_0 = mean(diff_prop_gt_0.2 > 0, na.rm = TRUE)
  )

ggplot(boot_diff, aes(x = diff_mean_r, fill = type)) +
  geom_histogram(data = filter(boot_diff, type == "bootstrap"), bins = 30, alpha = 0.6) +
  geom_vline(
    xintercept = boot_diff$diff_mean_r[boot_diff$type == "observed"],
    linewidth = 1
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(
    x = "Bootstrap N2 - daf-2 difference in mean pairwise correlation",
    y = "Bootstrap count",
    title = "Stability of Tepper Class I coordination difference"
  )




ggplot(summary_diff, aes(x = diff_mean_r, fill = label)) +
  geom_histogram(data = filter(summary_diff, label == "random"), bins = 30, alpha = 0.6) +
  geom_vline(
    xintercept = summary_diff$diff_mean_r[summary_diff$label == "observed"],
    linewidth = 1
  ) +
  theme_classic() +
  labs(
    x = "N2 - daf-2 difference in mean pairwise gene correlation",
    y = "Random gene-set count",
    title = "Observed Tepper Class I coordination compared with random gene sets"
  )







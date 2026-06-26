suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

qcut <- 0.05
dir.create("output", showWarnings = FALSE)

paper <- read.csv("data/paper_age_day8_variance_inputs.csv", check.names = FALSE)
stv <- read.csv("data/starvation_gamlss_inputs.csv", check.names = FALSE) |>
  mutate(
    stv_sig = !is.na(starvation_dispersion_fdr) &
      starvation_dispersion_fdr < qcut &
      recommended_use == "retain_for_primary_interpretation",
    stv_dir = if_else(starvation_dispersion_effect < 0,
                      "lower dispersion in daf-2", "higher dispersion in daf-2")
  )

age <- paper |>
  filter(comparison == "Day_8_vs_1", group %in% c("n2_sw", "daf2_sw")) |>
  transmute(
    Gene, Symbol,
    query_group = recode(group, n2_sw = "N2", daf2_sw = "daf-2"),
    variance_effect, variance_fdr,
    sig = !is.na(variance_fdr) & variance_fdr < qcut,
    dir = if_else(variance_effect > 0, "higher variance with age", "lower variance with age")
  )

day8 <- paper |>
  filter(comparison == "Strain_QZ120_vs_QZ0", group == "d8_sw") |>
  transmute(
    Gene, Symbol, variance_effect, variance_fdr,
    sig = !is.na(variance_fdr) & variance_fdr < qcut,
    dir = if_else(variance_effect > 0, "higher variance in daf-2", "lower variance in daf-2")
  )

query_sets <- bind_rows(
  crossing(age, query_dir = c("higher variance with age", "lower variance with age")) |>
    transmute(Gene, query_set = paste(query_group, query_dir, sep = ": "),
              query_sig = sig & dir == query_dir),
  crossing(day8, query_dir = c("higher variance in daf-2", "lower variance in daf-2")) |>
    transmute(Gene, query_set = paste("Day 8", query_dir, sep = ": "),
              query_sig = sig & dir == query_dir)
)

enrich_one <- function(query_set_i, stv_dir_i) {
  d <- inner_join(filter(query_sets, query_set == query_set_i), stv, by = "Gene")
  a <- sum(d$query_sig & d$stv_sig & d$stv_dir == stv_dir_i, na.rm = TRUE)
  m <- sum(d$query_sig, na.rm = TRUE)
  n <- sum(d$stv_sig & d$stv_dir == stv_dir_i, na.rm = TRUE)
  N <- nrow(d)
  data.frame(
    query_set = query_set_i, starvation_set = stv_dir_i, universe = N,
    overlap = a, query_n = m, starvation_n = n, expected = m * n / N,
    odds_ratio = ((a + .5) * (N - m - n + a + .5)) / ((m - a + .5) * (n - a + .5)),
    p = phyper(a - 1, m, N - m, n, lower.tail = FALSE)
  )
}

query_order <- c("daf-2: lower variance with age", "daf-2: higher variance with age",
                 "N2: lower variance with age", "N2: higher variance with age",
                 "Day 8: lower variance in daf-2", "Day 8: higher variance in daf-2")

ov <- expand_grid(
  query_set = rev(query_order),
  stv_dir = c("lower dispersion in daf-2", "higher dispersion in daf-2")
) |>
  rowwise() |>
  do(enrich_one(.$query_set, .$stv_dir)) |>
  ungroup() |>
  mutate(
    fdr = p.adjust(p, "BH"),
    neglog10_fdr = -log10(fdr),
    query_set = factor(query_set, levels = query_order),
    starvation_set = factor(starvation_set, c("lower dispersion in daf-2", "higher dispersion in daf-2"))
  )

write.csv(ov, "output/directional_overlap_enrichment.csv", row.names = FALSE)

counts <- bind_rows(
  age |> filter(sig) |> count(set = paste(query_group, dir, sep = ": ")),
  day8 |> filter(sig) |> count(set = paste("Day 8", dir, sep = ": ")),
  stv |> filter(stv_sig) |> count(set = paste("Starvation", stv_dir, sep = ": "))
) |>
  mutate(set = factor(set, levels = set[order(n)]))

p1 <- ggplot(counts, aes(n, set)) +
  geom_col(fill = "#4C78A8", width = .72) +
  geom_text(aes(label = scales::comma(n)), hjust = -0.05, size = 2.7) +
  scale_x_continuous(expand = expansion(mult = c(0, .16))) +
  labs(x = "Significant genes", y = NULL, title = "A. Differential-variance call sets") +
  theme_classic(base_size = 9)

p2 <- ggplot(ov, aes(starvation_set, query_set, fill = neglog10_fdr)) +
  geom_tile(color = "white", linewidth = .6) +
  geom_text(aes(label = paste0("n=", overlap, "\nFDR=", signif(fdr, 2))), size = 2.35) +
  scale_fill_gradient(low = "#E8E8E8", high = "#2F6F9F") +
  labs(x = NULL, y = NULL, fill = "-log10 FDR", title = "B. Directional overlap enrichment") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

fig <- (p1 | p2) +
  plot_annotation(
    title = "Variance gene sets compared with starvation N2 vs daf-2 dispersion",
    subtitle = "Panel A shows full significant call-set sizes; panel B tests directional overlap within the common tested gene universe."
  )

ggsave("output/variance_age_starvation_AB_suppfig.pdf", fig, width = 11, height = 4.6, device = "pdf")

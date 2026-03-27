library(gplots)
library(dplyr)
library(readr)
library(eulerr)
library(grid)
library(tidyverse)
library(gridExtra)

GAMLSS_concordantGenes<-read.csv("concordant_dispersion_genes_residual_dispersion_summary.csv",header = T) #Condordant genes are found in Supplementary File 1 tab
GeneLists<-read.csv("GeneListsForVarAnalysis.csv",header = T) #For gene lists, import DAF-16 target lists from Supplementary File 1 tab
background<-read.csv("per_gene_results_annotated.csv",header = T)
backgroundGenes<-background$Gene

GAMLSS_concordantGenes%>%
  filter(Strain=="N2" & Replicate=="I")%>%
  group_by(global_disp_direction)%>%
  summarize(count=n())

LowerDispersionDaf2<-GAMLSS_concordantGenes%>%
  filter(Strain=="N2" & Replicate=="I" & global_disp_direction=="lower_disp_in_daf2")

HigherDispersionDaf2<-GAMLSS_concordantGenes%>%
  filter(Strain=="N2" & Replicate=="I" & global_disp_direction=="higher_disp_in_daf2")

HigherDispersionDaf2Genes<-HigherDispersionDaf2$Gene
LowerDispersionDaf2Genes<-LowerDispersionDaf2$Gene

run_hypergeom_overlap <- function(target_genes, gene_list_df, background_genes, target_name) {
  
  N <- length(background_genes)
  n <- length(target_genes)
  
  results_list <- lapply(names(gene_list_df), function(colname) {
    
    test_genes <- gene_list_df[[colname]] %>%
      na.omit() %>%
      unique() %>%
      intersect(background_genes)
    
    K <- length(test_genes)
    overlap_genes <- intersect(target_genes, test_genes)
    k <- length(overlap_genes)
    
    pval <- if (K == 0 || n == 0) {
      NA
    } else {
      phyper(q = k - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
    }
    
    data.frame(
      target_list = target_name,
      gene_list_name = colname,
      background_size = N,
      target_size = n,
      comparison_size = K,
      overlap_size = k,
      expected_overlap = (n * K) / N,
      fold_enrichment = ifelse((n * K) / N > 0, k / ((n * K) / N), NA),
      p_value = pval,
      overlap_genes = paste(sort(overlap_genes), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  
  results <- bind_rows(results_list)
  results$FDR <- p.adjust(results$p_value, method = "BH")
  results %>% arrange(p_value)
}

LowerDispersion_overlap_results <- run_hypergeom_overlap(
  target_genes = LowerDispersionDaf2Genes,
  gene_list_df = GeneLists,
  background_genes = backgroundGenes,
  target_name = "LowerDispersionDaf2Genes"
)

HigherDispersion_overlap_results <- run_hypergeom_overlap(
  target_genes = HigherDispersionDaf2Genes,
  gene_list_df = GeneLists,
  background_genes = backgroundGenes,
  target_name = "HigherDispersionDaf2Genes"
)

All_overlap_results <- bind_rows(LowerDispersion_overlap_results,
                                 HigherDispersion_overlap_results)

All_overlap_results$FDR_all_tests <- p.adjust(All_overlap_results$p_value, method = "BH")

write.csv(LowerDispersion_overlap_results,
          "LowerDispersionDaf2_hypergeometric_overlaps.csv",
          row.names = FALSE)

write.csv(HigherDispersion_overlap_results,
          "HigherDispersionDaf2_hypergeometric_overlaps.csv",
          row.names = FALSE)

write.csv(All_overlap_results,
          "All_dispersion_hypergeometric_overlaps.csv",
          row.names = FALSE)

ConcordantDispersionDaf2<-GAMLSS_concordantGenes%>%
  filter(Strain=="N2" & Replicate=="I")


ConcordantDispersionDaf2_WithBackground<-inner_join(ConcordantDispersionDaf2,background,join_by("Gene"))

ConcordantDispersionDaf2_WithBackground2<-right_join(ConcordantDispersionDaf2,background,join_by("Gene"))
tail(ConcordantDispersionDaf2_WithBackground2)

PlotA<-ConcordantDispersionDaf2_WithBackground2%>%
  ggplot(aes(x=mean_effect,y=disp_effect.y,color=classification))+
  geom_point(alpha=0.3)+
  scale_color_manual(values = c("purple","red","blue","gray"))+
  xlim(-6,6)+ylim(-6,6)+
  theme_classic(base_size = 20)+
  theme(aspect.ratio = 1)+geom_hline(yintercept = 0,linetype="dotted")+
  labs(x="Mean effect",y="Dispersion effect")


ConcordantDispersionDaf2_WithBackground2%>%
  group_by(classification)%>%
  summarize(count=n())
  

PlotB<-ConcordantDispersionDaf2_WithBackground2%>%
  ggplot(aes(x=mean_effect))+
  geom_density()+ 
  xlim(-5,5)+geom_vline(xintercept = 0,linetype="dashed")+
  theme_classic(base_size = 20)+
  theme(aspect.ratio = 1)+
  labs(x="Mean effect", y="Density")

PlotC<-ConcordantDispersionDaf2_WithBackground2%>%
  ggplot(aes(x=disp_effect.y))+
  geom_density()+ 
  xlim(-5,5)+geom_vline(xintercept = 0,linetype="dashed")+
  theme_classic(base_size = 20)+
  theme(aspect.ratio = 1)+
  labs(x="Dispersion effect",y="Density")

head(ConcordantDispersionDaf2_WithBackground2)

PlotD<-ConcordantDispersionDaf2_WithBackground2%>%
  filter(mean_FDR<0.05)%>%
  ggplot(aes(x=mean_effect))+
  geom_density(color="blue")+ 
  xlim(-5,5)+geom_vline(xintercept = 0,linetype="dashed")+
  theme_classic(base_size = 20)+
  theme(aspect.ratio = 1)+
  labs(x="Mean effect", y="Density")

ConcordantDispersionDaf2_WithBackground2%>%
  filter(disp_FDR.y<0.05)%>%
  mutate(direction=if_else(disp_effect.y<0,"negative","positive"))%>%
  group_by(direction)%>%
  summarize(count=n())%>%
  mutate(prop=count/sum(count))


PlotE<-ConcordantDispersionDaf2_WithBackground2%>%
  filter(disp_FDR.y<0.05)%>%
  ggplot(aes(x=disp_effect.y))+
  geom_density(color="red")+ 
  xlim(-5,5)+geom_vline(xintercept = 0,linetype="dashed")+
  theme_classic(base_size = 20)+
  theme(aspect.ratio = 1)+
  labs(x="Dispersion effect",y="Density")


PlotA2 <- PlotA + theme(legend.position = "none")
grid.arrange(PlotA2, PlotB, PlotC, ncol = 3, widths = c(1.2, 1, 1))
grid.arrange(PlotD,PlotE, ncol = 2, widths = c(1, 1))


wilcox.test(ConcordantDispersionDaf2_WithBackground2$disp_effect.y, mu = 0, alternative = "less")
wilcox.test(ConcordantDispersionDaf2_WithBackground2$mean_effect, mu = 0, alternative = "two.sided")

median(ConcordantDispersionDaf2_WithBackground2$mean_effect,na.rm = TRUE)
mean(ConcordantDispersionDaf2_WithBackground2$mean_effect,na.rm = TRUE)

median(ConcordantDispersionDaf2_WithBackground2$disp_effect.y,na.rm = TRUE)
mean(ConcordantDispersionDaf2_WithBackground2$disp_effect.y,na.rm = TRUE)


x <- sum(ConcordantDispersionDaf2_WithBackground2$disp_effect.y < 0, na.rm = TRUE)
n <- sum(!is.na(ConcordantDispersionDaf2_WithBackground2$disp_effect.y))
binom.test(x, n, p = 0.5, alternative = "greater")


binom_mean <- binom.test(
  sum(ConcordantDispersionDaf2_WithBackground2$mean_effect < 0, na.rm = TRUE),
  sum(!is.na(ConcordantDispersionDaf2_WithBackground2$mean_effect)),
  p = 0.5,
  alternative = "two.sided"
)

ConcordantDispersionDaf2_WithBackground%>%
  mutate(mean_FDR_sig=if_else(mean_FDR<0.05,"sig mean","NS mean"))%>%
  filter(mean_FDR_sig=="NS mean" & disp_effect.x< -2 & prop_zero_N2<0.05 & prop_zero_daf2<0.05)
  
sig_lists <- LowerDispersion_overlap_results %>%
  filter(FDR < 0.05) %>%
  pull(gene_list_name)

for (this_list in names(GeneLists)) {
  
  overlap_string <- LowerDispersion_overlap_results %>%
    filter(gene_list_name == this_list) %>%
    pull(overlap_genes) %>%
    .[1]
  
  overlap_vec <- if (is.na(overlap_string) || overlap_string == "") {
    character(0)
  } else {
    str_split(overlap_string, pattern = ",", simplify = FALSE)[[1]] %>%
      str_trim()
  }
  
  p <- ConcordantDispersionDaf2_WithBackground %>%
    mutate(
      overlap_status = if_else(Gene %in% overlap_vec, "in overlap", "not in overlap"),
      overlap_status = factor(overlap_status, levels = c("not in overlap", "in overlap"))
      ) %>%
    ggplot(aes(x = disp_effect.x, color = overlap_status)) +
    stat_ecdf(linewidth=1)+
    theme_classic(base_size = 12) +
    theme(aspect.ratio = 1) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    ggtitle(this_list)
  
  ggsave(
    filename = paste0("CumulDist_Var_", make.names(this_list), ".pdf"),
    plot = p,
    width = 5,
    height = 5
  )
}

ConcordantDispersionDaf2_WithBackground %>%
  mutate(
    overlap_status = if_else(Gene %in% overlap_vec, "in overlap", "not in overlap"),
    overlap_status = factor(overlap_status, levels = c("not in overlap", "in overlap"))
  ) %>%
  group_by(overlap_status)%>%
  summarize(count=n())

sig_lists <- LowerDispersion_overlap_results %>%
  filter(FDR < 0.05) %>%
  pull(gene_list_name)

for (this_list in names(GeneLists)) {
  
  gene_list_vec <- GeneLists[[this_list]] %>%
    na.omit() %>%
    unique()
  
  gene_list_vec_bg <- intersect(gene_list_vec, backgroundGenes)
  
  plot_df <- ConcordantDispersionDaf2_WithBackground2 %>%
    filter(Gene %in% backgroundGenes) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    mutate(
      list_status = if_else(Gene %in% gene_list_vec_bg, "in gene list", "background only"),
      list_status = factor(list_status, levels = c("background only", "in gene list"))
    )
  
  x_in_mean <- plot_df %>%
    filter(list_status == "in gene list") %>%
    pull(mean_effect)
  
  x_out_mean <- plot_df %>%
    filter(list_status == "background only") %>%
    pull(mean_effect)
  
  ks_mean <- ks.test(x_in_mean, x_out_mean)
  
  x_in_disp <- plot_df %>%
    filter(list_status == "in gene list") %>%
    pull(disp_effect.y)
  
  x_out_disp <- plot_df %>%
    filter(list_status == "background only") %>%
    pull(disp_effect.y)
  
  ks_disp <- ks.test(x_in_disp, x_out_disp)
  
  p_mean <- ggplot(plot_df, aes(x = mean_effect, color = list_status)) +
    stat_ecdf(linewidth = 1) +
    coord_cartesian(xlim = c(-4, 4)) +
    scale_color_manual(values = c("gray", "blue")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic(base_size = 12) +
    theme(aspect.ratio = 1) +
    ggtitle(
      paste("Cumulative distribution of mean effect:", this_list),
      subtitle = paste0(
        "KS p = ", signif(ks_mean$p.value, 3),
        "; n = ", length(x_in_mean), " vs ", length(x_out_mean)
      )
    ) +
    xlab("mean_effect") +
    ylab("Cumulative fraction")
  
  p_disp <- ggplot(plot_df, aes(x = disp_effect.y, color = list_status)) +
    stat_ecdf(linewidth = 1) +
    coord_cartesian(xlim = c(-4, 4)) +
    scale_color_manual(values = c("gray", "red")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic(base_size = 12) +
    theme(aspect.ratio = 1) +
    ggtitle(
      paste("Cumulative distribution of dispersion effect:", this_list),
      subtitle = paste0(
        "KS p = ", signif(ks_disp$p.value, 3),
        "; n = ", length(x_in_disp), " vs ", length(x_out_disp)
      )
    ) +
    xlab("disp_effect") +
    ylab("Cumulative fraction")
  
  ggsave(paste0("Updated_ECDF_mean_", make.names(this_list), ".pdf"), p_mean, width = 5, height = 5)
  ggsave(paste0("Updated_ECDF_disp_", make.names(this_list), ".pdf"), p_disp, width = 5, height = 5)
}




overlap_results <- read_csv("LowerDispersionDaf2_hypergeometric_overlaps.csv")

overlap_results <- read_csv("LowerDispersionDaf2_hypergeometric_overlaps.csv")

make_euler_from_results <- function(results_df, list1, list2,
                                    target_name = "Reduced dispersion in daf-2",
                                    save_file = NULL) {
  
  row1 <- results_df %>% filter(gene_list_name == list1)
  row2 <- results_df %>% filter(gene_list_name == list2)
  
  if (nrow(row1) != 1 || nrow(row2) != 1) {
    stop("Could not uniquely find one or both gene lists in results_df.")
  }
  
  A  <- row1$target_size[1]
  B  <- row1$comparison_size[1]
  C  <- row2$comparison_size[1]
  AB <- row1$overlap_size[1]
  AC <- row2$overlap_size[1]
  BC <- 0
  ABC <- 0
  
  A_only <- A - AB - AC + ABC
  B_only <- B - AB - BC + ABC
  C_only <- C - AC - BC + ABC
  
  if (A_only < 0 || B_only < 0 || C_only < 0) {
    stop("One or more unique region sizes are negative. Check the input overlaps.")
  }
  
  fit <- euler(c(
    "ReducedDisp" = A_only,
    "ClassI" = B_only,
    "ClassII" = C_only,
    "ReducedDisp&ClassI" = AB - ABC,
    "ReducedDisp&ClassII" = AC - ABC,
    "ClassI&ClassII" = BC - ABC,
    "ReducedDisp&ClassI&ClassII" = ABC
  ))
  
  if (!is.null(save_file)) {
    png(save_file, width = 1800, height = 1800, res = 300)
    on.exit(dev.off(), add = TRUE)
  }
  
  plot(
    fit,
    fills = list(fill = c("gray70", "steelblue", "tomato"), alpha = 0.5),
    edges = list(col = "black"),
    labels = list(labels = c(target_name, list1, list2)),
    quantities = TRUE,
    main = paste(target_name, "\nvs", list1, "and", list2)
  )
}


make_euler_from_results(
  overlap_results,
  list1 = "Tepper_ClassI",
  list2 = "Tepper_ClassII",
  save_file = "Euler_Tepper_ClassI_ClassII.png"
)

make_euler_from_results(
  overlap_results,
  list1 = "DAF16_ClassI_DOWN",
  list2 = "DAF16_ClassII_UP",
  save_file = "Euler_Tepper_ClassI_ClassII.png"
)

make_euler_from_results(
  overlap_results,
  list1 = "DOWN_daf16_starvationL1_Kaplan",
  list2 = "UP_daf16_starvationL1_Kaplan",
  save_file = "Euler_Tepper_ClassI_ClassII.png"
)



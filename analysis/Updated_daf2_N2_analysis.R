library(tidyverse)
library(edgeR)
library(clusterProfiler)
library(org.Ce.eg.db)
library(dplyr)
library(tidyr)
library(readr)
library(forcats)
library(stringr)


CELseqUMICounts<-read.csv("umi_count_matrix_all_combined_3sets.csv",header = T)
CELseqUMICountsLongFormat<-CELseqUMICounts%>%pivot_longer(!Gene,names_to = "sample_id",values_to = "counts")

SampleKey<-read.csv("SampleNames.csv",header = T)

CELseqCountsLongKey<-merge(CELseqUMICountsLongFormat,SampleKey,by.x="sample_id",by.y="SampleName")
length(unique(CELseqCountsLongKey$sample_id)) #288

#count how many reads per sample within this group

stats1<-CELseqCountsLongKey%>%
  group_by(sample_id,Strain,Replicate,Strain_Rep,Barcode)%>%
  mutate(CountsPerSample=sum(counts))%>%
  filter(CountsPerSample>250000)

stats2<-stats1%>%
  group_by(Strain,Replicate,Strain_Rep)%>%
  filter(Gene=="WBGene00000001")%>%
  summarize(numberSamples=n())

stats3<-stats1%>%
  group_by(Strain,Replicate,Strain_Rep,sample_id)%>%
  filter(Gene=="WBGene00000001")%>%
  select(-Gene,-counts)

counts_filtered<-stats1

head(counts_filtered)
counts_filtered2<-counts_filtered%>%
  ungroup()%>%
  select(sample_id,Gene,counts)


counts_filtered3<-counts_filtered2%>%
  pivot_wider(names_from = sample_id,values_from = counts)

counts_filtered3<-as.data.frame(counts_filtered3)

rownames(counts_filtered3)<-counts_filtered3$Gene

WS273_geneNames<-read.csv("WS273_geneNames.csv",header = T)

y <- DGEList(counts=counts_filtered3, genes=rownames(counts_filtered3))


# Filter for protein coding genes
#y <- y[y$genes$genes %in% protein_coding_genes,]
#print(paste("After filtering:", nrow(y), "genes remaining"))

#CPM matrix
cpm_matrix<-cpm(y)
dim(cpm_matrix)[2]


prop_of_samples<-0.6

# Filter for expressed genes (CPM > 1 in at least 4 samples)
keep <- rowSums(cpm_matrix > 1) >= dim(cpm_matrix)[2]*prop_of_samples
y <- y[keep,]
print(paste("After filtering:", nrow(y), "genes remaining"))

# Normalize
print("Normalizing data...")
y <- calcNormFactors(y, method="TMM")


# Calculate CPM
print("Calculating CPM...")
cpm_matrix2 <- cpm(y)
head(cpm_matrix2)

# Create design matrix
print("Creating design matrix...")
stats3$Strain_Rep
group <- factor(stats3$Strain_Rep)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Estimate dispersion
print("Estimating dispersion...")
y <- estimateDisp(y, design)

# Fit model
print("Fitting model...")
fit <- glmQLFit(y, design)


dim(cpm_matrix2)[2]

cpm_d2_df<-data.frame(cpm_matrix2)

cpm_d2_df_orig<-cpm_d2_df[,1:dim(cpm_matrix2)[2]]
cpm_d2_df_orig$gene<-rownames(cpm_d2_df_orig)

cpm_d2_df_long<-cpm_d2_df_orig%>%pivot_longer(1:dim(cpm_matrix2)[2],names_to = "SampleIDs",values_to = "CPMs")

cpm_d2_df_long2<-merge(cpm_d2_df_long,SampleKey,by.x="SampleIDs",by.y="SampleName")

cpm_d2_df_long2<-cpm_d2_df_long2%>%
  select(SampleIDs,gene,CPMs,Strain_Rep,Strain,Replicate)

StrainRep_MeanVar<-cpm_d2_df_long2%>%
  group_by(Strain_Rep,gene,Strain,Replicate)%>%
  summarize(Mean=mean(CPMs),Variance=var(CPMs))%>%
  mutate(Log10mean=log10(Mean),Log10variance=log10(Variance))%>%
  filter(Log10mean>0)%>%
  filter(Log10variance>0)


nested_data <- StrainRep_MeanVar %>%
  ungroup()%>%
  group_by(Replicate)%>% #CHANGED THIS TO REPLICATE INSTEAD OF STRAIN_REP
  nest()

fitted_models <- nested_data %>%
  mutate(model = map(data, ~ loess(Log10variance ~ Log10mean, data = ., span = 0.75)))

results <- fitted_models %>%
  mutate(fitted_values = map(model, "fitted")) 

final_data <- results %>%
  unnest(cols = c(data, fitted_values))%>%
  mutate()

write.table(final_data[,c(1:8,10:14)],"MeanVar_Fit_Updated021126.txt",quote = F,sep = "\t")

final_data<-final_data%>%
  mutate(VarDistance=Log10variance - fitted_values)%>%
  mutate(meanDistance=mean(VarDistance),SDDistance=sd(VarDistance))%>%
  mutate(Zscore=(VarDistance-meanDistance)/SDDistance)

final_data2<-final_data%>%
  select(-model)

# expected log10 variance per gene & group from final_data
exp_sd_tbl <- final_data %>%
  dplyr::select(Strain_Rep, gene, exp_log10var = fitted_values) %>%
  dplyr::distinct() %>%
  dplyr::mutate(sd_exp = 10^(exp_log10var/2))

# Join expected SD, center by group mean CPM, divide by expected SD
R_long <- cpm_d2_df_long2 %>%
  dplyr::inner_join(exp_sd_tbl, by = c("Strain_Rep","gene")) %>%
  dplyr::group_by(Strain_Rep, gene) %>%
  dplyr::mutate(mu = mean(CPMs, na.rm = TRUE),
                R_units = (CPMs - mu) / sd_exp) %>%
  dplyr::ungroup()

head(R_long)

R_long<-R_long%>%
  mutate(abs_R=abs(R_units))

head(R_long_new)
length(unique(R_long_new$SampleIDs)) #248

ggplot(filter(R_long, gene %in% c("WBGene00001852")),
       aes(x = abs_R, color = Strain)) +
  stat_ecdf(geom = "step", linewidth = 0.5, pad = FALSE,alpha=0.7) +
  facet_grid(Replicate.x~ gene, scales = "fixed") + 
  coord_cartesian(xlim = c(0, 5)) +theme(aspect.ratio = 1)+
  labs(x = "|R_units| (distance from mean in expected-SD)",
       y = "F(|R|)", color = "Group",
       title = "Right-tail ECDFs") +
  theme_bw()+theme(aspect.ratio = 1)


write.table(R_long,"ExpectedSD_units_UpdatedByReps_021126.txt",quote = F,sep = "\t")

R_long_new<-read.table("ExpectedSD_units_UpdatedByReps_021126.txt")

#genes from paper
R_long_subset<-R_long_new%>%
  filter(gene %in% c("WBGene00012544","WBGene00001852","WBGene00020297","WBGene00021253","WBGene00000230","WBGene00011273","WBGene00012354","WBGene00017121"))

ggplot(filter(R_long_new, gene %in% c("WBGene00000230","WBGene00011273","WBGene00012354","WBGene00017121")),
              aes(x = abs_R, color = Strain)) +
  stat_ecdf(geom = "step", linewidth = 0.5, pad = FALSE,alpha=0.7) +
  facet_wrap(~ gene, scales = "fixed",ncol=4) +
  scale_color_manual(values = c("red","black"))+
  coord_cartesian(xlim = c(0, 5)) +theme(aspect.ratio = 1)+
  labs(x = "|R_units| (distance from mean in expected-SD)",
       y = "F(|R|)", color = "Group") +
  theme_classic()+theme(aspect.ratio = 1)

ggplot(filter(R_long_new, gene %in% c("WBGene00000230","WBGene00011273","WBGene00012354","WBGene00017121")),
       aes(x = Strain_Rep,y=R_units, color = Strain)) +
  geom_violin()+geom_jitter(width = 0.15,alpha=0.8)+
  facet_wrap(~ gene, scales = "free_y",ncol=4) +geom_hline(yintercept = 0,linetype="dashed")+
  scale_color_manual(values = c("red","black"))+
  theme_classic()+theme(aspect.ratio = 1)

grid.arrange(PlotF,PlotG,ncol=1)

ggplot(filter(R_long_new, gene %in% c("WBGene00000230","WBGene00011273","WBGene00012354")),
       aes(x = abs_R, color = Strain)) +
  stat_ecdf(geom = "step", linewidth = 0.5, pad = FALSE,alpha=0.7) +
  facet_wrap(gene ~ Replicate.x, scales = "fixed") + 
  coord_cartesian(xlim = c(0, 5)) +theme(aspect.ratio = 1)+
  labs(x = "|R_units| (distance from mean in expected-SD)",
       y = "F(|R|)", color = "Group",
       title = "Right-tail ECDFs") +
  theme_bw()+theme(aspect.ratio = 1)

Gene_q90_Mean<-R_long_new%>%
  group_by(gene,Strain_Rep,Strain,Replicate.x)%>%
  summarize(mean_CPM=mean(CPMs),q90=quantile(abs_R, probs = 0.9, na.rm = TRUE))

Gene_q90_Mean2<-R_long_new%>%
  group_by(gene,Strain)%>%
  summarize(mean_CPM=mean(CPMs),q90=quantile(abs_R, probs = 0.9, na.rm = TRUE))

subset(Gene_q90_Mean,gene=="WBGene00001852")
subset(Gene_q90_Mean2,gene=="WBGene00001852")





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



lower_genes <- unique(na.omit(LowerDispersionDaf2Genes))
write.table(lower_genes,"LowDispersionGenes.txt",quote = F,sep = "\t")
higher_genes <- unique(na.omit(HigherDispersionDaf2Genes))
universe_genes <- unique(na.omit(backgroundGenes))

ego_lower_BP <- enrichGO(
  gene          = lower_genes,
  universe      = universe_genes,
  OrgDb         = org.Ce.eg.db,
  keyType       = "WORMBASE",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = FALSE
)

ego_higher_BP <- enrichGO(
  gene          = higher_genes,
  universe      = universe_genes,
  OrgDb         = org.Ce.eg.db,
  keyType       = "WORMBASE",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = FALSE
)


ego_lower_MF <- enrichGO(lower_genes, universe = universe_genes, OrgDb = org.Ce.eg.db,
                         keyType = "WORMBASE", ont = "MF", pAdjustMethod = "BH")
ego_lower_CC <- enrichGO(lower_genes, universe = universe_genes, OrgDb = org.Ce.eg.db,
                         keyType = "WORMBASE", ont = "CC", pAdjustMethod = "BH")

ego_higher_MF <- enrichGO(higher_genes, universe = universe_genes, OrgDb = org.Ce.eg.db,
                          keyType = "WORMBASE", ont = "MF", pAdjustMethod = "BH")
ego_higher_CC <- enrichGO(higher_genes, universe = universe_genes, OrgDb = org.Ce.eg.db,
                          keyType = "WORMBASE", ont = "CC", pAdjustMethod = "BH")

lower_BP_res <- as.data.frame(ego_lower_BP)
higher_BP_res <- as.data.frame(ego_higher_BP)
higher_BP_res

write.csv(lower_BP_res, "GO_BP_lower_dispersion_daf2.csv", row.names = FALSE)

dotplot(ego_lower_BP, showCategory = 15, title = "Lower dispersion in daf-2: GO BP")
dotplot(ego_higher_BP, showCategory = 15, title = "Higher dispersion in daf-2: GO BP")
dotplot(ego_lower_MF, showCategory = 15, title = "Lower dispersion in daf-2: GO MF")
dotplot(ego_lower_CC, showCategory = 15, title = "Lower dispersion in daf-2: GO CC")

lower_MF_res<-as.data.frame(ego_lower_MF)
write.csv(lower_MF_res, "GO_MF_lower_dispersion_daf2.csv", row.names = FALSE)

lower_CC_res<-as.data.frame(ego_lower_CC)
write.csv(lower_CC_res, "GO_CC_lower_dispersion_daf2.csv", row.names = FALSE)


Muscle_Age_Down <- read_csv("Muscle_Age_Down.csv") #Gene lists significant in each scRNA-seq analysis, FDR<0.05, Supp File 1
Muscle_Age_Up   <- read_csv("Muscle_Age_Up.csv")
Neuron_Age_Down <- read_csv("Neuron_Age_Down.csv")
Neuron_Age_Up   <- read_csv("Neuron_Age_Up.csv")

lookup_long <- WS273_geneNames %>%
  dplyr::select(WB_id, symbol, sequence) %>%
  pivot_longer(
    cols = c(symbol, sequence),
    names_to = "name_type",
    values_to = "gene_name"
  ) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  distinct(gene_name, .keep_all = TRUE)

Muscle_Age_Down_WB <- Muscle_Age_Down %>%
  left_join(lookup_long, by = c("gene" = "gene_name"))

Muscle_Age_Up_WB <- Muscle_Age_Up %>%
  left_join(lookup_long, by = c("gene" = "gene_name"))

Neuron_Age_Down_WB <- Neuron_Age_Down %>%
  left_join(lookup_long, by = c("gene" = "gene_name"))

Neuron_Age_Up_WB <- Neuron_Age_Up %>%
  left_join(lookup_long, by = c("gene" = "gene_name"))

Muscle_Age_Down_Genes <- Muscle_Age_Down_WB %>%
  pull(WB_id) %>%
  na.omit() %>%
  unique() %>%
  intersect(backgroundGenes)

Muscle_Age_Up_Genes <- Muscle_Age_Up_WB %>%
  pull(WB_id) %>%
  na.omit() %>%
  unique() %>%
  intersect(backgroundGenes)

Neuron_Age_Down_Genes <- Neuron_Age_Down_WB %>%
  pull(WB_id) %>%
  na.omit() %>%
  unique() %>%
  intersect(backgroundGenes)

Neuron_Age_Up_Genes <- Neuron_Age_Up_WB %>%
  pull(WB_id) %>%
  na.omit() %>%
  unique() %>%
  intersect(backgroundGenes)

mapping_summary <- tibble(
  list = c("Muscle_Age_Down", "Muscle_Age_Up", "Neuron_Age_Down", "Neuron_Age_Up"),
  input_n = c(nrow(Muscle_Age_Down), nrow(Muscle_Age_Up), nrow(Neuron_Age_Down), nrow(Neuron_Age_Up)),
  mapped_n = c(sum(!is.na(Muscle_Age_Down_WB$WB_id)),
               sum(!is.na(Muscle_Age_Up_WB$WB_id)),
               sum(!is.na(Neuron_Age_Down_WB$WB_id)),
               sum(!is.na(Neuron_Age_Up_WB$WB_id))),
  in_background_n = c(length(Muscle_Age_Down_Genes),
                      length(Muscle_Age_Up_Genes),
                      length(Neuron_Age_Down_Genes),
                      length(Neuron_Age_Up_Genes))
)

print(mapping_summary)

AgeGeneSets <- list(
  Muscle_Age_Down = Muscle_Age_Down_Genes,
  Muscle_Age_Up   = Muscle_Age_Up_Genes,
  Neuron_Age_Down = Neuron_Age_Down_Genes,
  Neuron_Age_Up   = Neuron_Age_Up_Genes
)

run_hypergeom_overlap_list <- function(target_genes, gene_sets, background_genes, target_name) {
  
  N <- length(background_genes)
  n <- length(target_genes)
  
  results_list <- lapply(names(gene_sets), function(set_name) {
    
    test_genes <- unique(na.omit(gene_sets[[set_name]]))
    test_genes <- intersect(test_genes, background_genes)
    
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
      gene_list_name = set_name,
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

LowerDispersion_Aging_overlap_results <- run_hypergeom_overlap_list(
  target_genes = LowerDispersionDaf2Genes,
  gene_sets = AgeGeneSets,
  background_genes = backgroundGenes,
  target_name = "LowerDispersionDaf2Genes"
)

HigherDispersion_Aging_overlap_results <- run_hypergeom_overlap_list(
  target_genes = HigherDispersionDaf2Genes,
  gene_sets = AgeGeneSets,
  background_genes = backgroundGenes,
  target_name = "HigherDispersionDaf2Genes"
)

write.csv(LowerDispersion_Aging_overlap_results,
          "LowerDispersion_Aging_overlap_results.csv",
          row.names = FALSE)

write.csv(HigherDispersion_Aging_overlap_results,
          "HigherDispersion_Aging_overlap_results.csv",
          row.names = FALSE)
LowerDispersion_Aging_overlap_results

LowerDispersion_Aging_overlap_results <- read_csv("LowerDispersion_Aging_overlap_results.csv")

plot_df <- LowerDispersion_Aging_overlap_results %>%
  mutate(
    gene_list_name = factor(
      gene_list_name,
      levels = c("Muscle_Age_Down", "Muscle_Age_Up", "Neuron_Age_Down", "Neuron_Age_Up")
    ),
    neglog10FDR = -log10(FDR),
    obs_exp = overlap_size / expected_overlap
  )

plot_df <- plot_df %>%
  mutate(gene_list_name = reorder(gene_list_name, fold_enrichment))

ggplot(plot_df, aes(x = fold_enrichment, y = gene_list_name)) +
  geom_point(aes(size = overlap_size, color = neglog10FDR)) +
  scale_color_gradientn(
    colours = rev(c("#DE6A6A", "#B07AA1", "#4A90D9")),
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(name = "Overlap size") +
  theme_bw(base_size = 12) +
  xlab("Fold enrichment") +
  ylab("") +
  ggtitle("Aging gene sets enriched among lower-dispersion daf-2 genes")

plot_df <- LowerDispersion_Aging_overlap_results %>%
  mutate(
    gene_list_name = factor(
      gene_list_name,
      levels = c("Muscle_Age_Down", "Muscle_Age_Up", "Neuron_Age_Down", "Neuron_Age_Up")
    )
  ) %>%
  tidyr::pivot_longer(
    cols = c(overlap_size, expected_overlap),
    names_to = "type",
    values_to = "count"
  ) %>%
  mutate(
    type = recode(type,
                  overlap_size = "Observed overlap",
                  expected_overlap = "Expected overlap")
  )

ggplot(plot_df, aes(x = gene_list_name, y = count, fill = type)) +
  geom_col(position = "dodge") +
  theme_classic(base_size = 12) +
  xlab("") +
  ylab("Number of genes") +
  ggtitle("Observed vs expected overlap with lower-dispersion daf-2 genes") +
  coord_flip()



WormBase_TissueEnrichment<-read.csv("WormBase_TissueEnrichment.csv",header = T)
head(WormBase_TissueEnrichment)
tissue_plot_df <- WormBase_TissueEnrichment %>%
  filter(Q.value<0.05)%>%
  arrange(Q.value) %>%
  slice_head(n = 12) %>%
  mutate(
    neglog10FDR = -log10(Q.value),
    Term = reorder(Term, Enrichment.Fold.Change)
  )

ggplot(tissue_plot_df, aes(x = Enrichment.Fold.Change, y = Term)) +
  geom_point(aes(size = Observed, color = neglog10FDR)) +
  scale_color_gradientn(
    colours = c("#4A90D9", "#B07AA1", "#DE6A6A"),
    name = expression(-log[10](Q.value))
  ) +
  scale_size_continuous(name = "Overlap size") +
  theme_classic(base_size = 12) +
  xlab("Fold enrichment") +
  ylab("") +
  ggtitle("WormBase tissue enrichments among lower-dispersion daf-2 genes")


TopTermSelectiveCollapse <- WormBase_TissueEnrichment %>%
  mutate(
    Term = str_trim(Term),
    DisplayTerm = case_when(
      str_detect(Term, "^(DA|DB|DD|VA|VB|VC|VD)\\d+$") ~ "ventral cord motor neurons",
      str_detect(Term, "^AS\\d+$") ~ "ventral cord motor neurons",
      Term %in% c("dorso-rectal ganglion", "retrovesicular ganglion", "preanal ganglion", "posterior lateral ganglion") ~ "ganglia",
      TRUE ~ Term
    ),
    neglog10Q = -log10(Q.value)
  ) %>%
  group_by(DisplayTerm) %>%
  slice_min(Q.value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(DisplayTerm = reorder(DisplayTerm, Enrichment.Fold.Change))

TopTermSelectiveCollapse%>%
arrange(Q.value) %>%
  slice_head(n = 20) %>%
ggplot( aes(x = Enrichment.Fold.Change, y = DisplayTerm)) +
  geom_point(aes(size = Observed, color = neglog10Q)) +
  scale_color_gradientn(
    colours = c("#4A90D9", "#B07AA1", "#DE6A6A"),
    name = expression(-log[10](Q.value))
  ) +
  scale_size_continuous(name = "Overlap size") +
  theme_bw(base_size = 12) +
  xlab("Fold enrichment") +
  ylab("") +
  ggtitle("WormBase tissue enrichments among lower-dispersion daf-2 genes")




AgeGeneSets <- list(
  Muscle_Age_Down = Muscle_Age_Down_Genes,
  Muscle_Age_Up   = Muscle_Age_Up_Genes,
  Neuron_Age_Down = Neuron_Age_Down_Genes,
  Neuron_Age_Up   = Neuron_Age_Up_Genes
)

for (this_list in names(AgeGeneSets)) {
  
  gene_list_vec_bg <- intersect(unique(na.omit(AgeGeneSets[[this_list]])), backgroundGenes)
  
  plot_df <- ConcordantDispersionDaf2_WithBackground2 %>%
    filter(Gene %in% backgroundGenes) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    filter(is.finite(disp_effect.y)) %>%
    mutate(
      list_status = if_else(Gene %in% gene_list_vec_bg, "in gene list", "background only"),
      list_status = factor(list_status, levels = c("background only", "in gene list"))
    )
  
  x_in  <- plot_df %>%
    filter(list_status == "in gene list") %>%
    pull(disp_effect.y)
  
  x_out <- plot_df %>%
    filter(list_status == "background only") %>%
    pull(disp_effect.y)
  
  ks_res <- ks.test(x_in, x_out)
  wilcox_res <- wilcox.test(x_in, x_out, alternative = "less")
  
  p <- ggplot(plot_df, aes(x = disp_effect.y, color = list_status)) +
    stat_ecdf(linewidth = 1) +
    coord_cartesian(xlim = c(-5, 5)) +
    scale_color_manual(values = c("background only" = "gray70", "in gene list" = "red")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic(base_size = 12) +
    theme(aspect.ratio = 1) +
    ggtitle(
      paste("Dispersion effects:", this_list),
      subtitle = paste0(
        "KS p = ", signif(ks_res$p.value, 3),
        "; Wilcoxon p = ", signif(wilcox_res$p.value, 3),
        "; n = ", length(x_in), " vs ", length(x_out)
      )
    ) +
    xlab("Dispersion effect") +
    ylab("Cumulative fraction")
  
  ggsave(
    filename = paste0("ECDF_disp_", make.names(this_list), ".pdf"),
    plot = p,
    width = 5,
    height = 5
  )
}


AgeGeneSets <- list(
  Muscle_Age_Down = Muscle_Age_Down_Genes,
  Muscle_Age_Up   = Muscle_Age_Up_Genes,
  Neuron_Age_Down = Neuron_Age_Down_Genes,
  Neuron_Age_Up   = Neuron_Age_Up_Genes
)

combined_df <- map_dfr(names(AgeGeneSets), function(this_list) {
  gene_list_vec_bg <- intersect(unique(na.omit(AgeGeneSets[[this_list]])), backgroundGenes)
  
  ConcordantDispersionDaf2_WithBackground2 %>%
    filter(Gene %in% backgroundGenes) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    filter(is.finite(disp_effect.y)) %>%
    mutate(
      gene_set = this_list,
      list_status = if_else(Gene %in% gene_list_vec_bg, "in gene list", "background only"),
      list_status = factor(list_status, levels = c("background only", "in gene list"))
    )
})

library(dplyr)
library(ggplot2)

pval_df <- combined_df %>%
  group_by(gene_set) %>%
  summarize(
    p_value = wilcox.test(
      disp_effect.y[list_status == "in gene list"],
      disp_effect.y[list_status == "background only"],
      alternative = "less"
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("Wilcoxon p = ", signif(p_value, 3)),
    y = 1.35
  )

ggplot(combined_df, aes(x = gene_set, y = disp_effect.y, color = list_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(
    data = pval_df,
    aes(x = gene_set, y = y, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_color_manual(values = c("background only" = "gray70", "in gene list" = "red")) +
  coord_cartesian(ylim = c(-3, 1.5)) +
  theme_classic(base_size = 12) +
  theme(aspect.ratio = 1, legend.position = "none") +
  xlab("") +
  ylab("Dispersion effect")



get_go_genes <- function(go_df, term_name, background_genes) {
  genes <- go_df %>%
    filter(Description == term_name) %>%
    pull(geneID) %>%
    .[1]
  
  if (is.na(genes) || genes == "") {
    return(character(0))
  }
  
  genes <- strsplit(genes, "/")[[1]]
  intersect(unique(genes), background_genes)
}

make_upset_df <- function(set_list, background_genes) {
  upset_df <- tibble(Gene = background_genes)
  for (nm in names(set_list)) {
    upset_df[[nm]] <- as.integer(background_genes %in% set_list[[nm]])
  }
  as.data.frame(upset_df)
}

LowerDisp_set <- intersect(unique(LowerDispersionDaf2Genes), backgroundGenes)
MuscleAgeDown_set <- intersect(unique(Muscle_Age_Down_Genes), backgroundGenes)
NeuronAgeDown_set <- intersect(unique(Neuron_Age_Down_Genes), backgroundGenes)

Energy_set <- get_go_genes(
  lower_BP_res,
  "generation of precursor metabolites and energy",
  backgroundGenes
)

plot_sets_energy <- list(
  LowerDispersion = LowerDisp_set,
  Muscle_Age_Down = MuscleAgeDown_set,
  Neuron_Age_Down = NeuronAgeDown_set,
  Energy_GO = Energy_set
)

upset_df_energy <- make_upset_df(plot_sets_energy, backgroundGenes)

head(upset_df_energy)

calc_atleast_intersection <- function(df, included_sets) {
  
  N <- nrow(df)
  
  observed <- df %>%
    filter(if_all(all_of(included_sets), ~ . == 1)) %>%
    nrow()
  
  p_included <- sapply(included_sets, function(s) mean(df[[s]] == 1))
  p_expected <- prod(p_included)
  expected <- N * p_expected
  
  pval <- binom.test(observed, N, p = p_expected, alternative = "greater")$p.value
  
  tibble(
    intersection = paste(included_sets, collapse = " & "),
    observed = observed,
    expected = expected,
    fold_enrichment = observed / expected,
    p_value = pval
  )
}



intersection_results <- bind_rows(
  calc_atleast_intersection(upset_df_energy, c("LowerDispersion", "Muscle_Age_Down")),
  calc_atleast_intersection(upset_df_energy, c("LowerDispersion", "Neuron_Age_Down")),
  calc_atleast_intersection(upset_df_energy, c("LowerDispersion", "Energy_GO")),
  calc_atleast_intersection(upset_df_energy, c("LowerDispersion", "Muscle_Age_Down", "Energy_GO")),
  calc_atleast_intersection(upset_df_energy, c("LowerDispersion", "Neuron_Age_Down", "Energy_GO")),
  calc_atleast_intersection(upset_df_energy, c("LowerDispersion", "Muscle_Age_Down", "Neuron_Age_Down")),
  calc_atleast_intersection(upset_df_energy, c("LowerDispersion", "Muscle_Age_Down", "Neuron_Age_Down", "Energy_GO"))
) %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

print(intersection_results)

intersection_results %>%
  mutate(intersection = reorder(intersection, fold_enrichment)) %>%
  ggplot(aes(x = fold_enrichment, y = intersection)) +
  geom_col(fill = "gray40") +
  geom_text(
    aes(label = paste0("obs=", observed,
                       ", exp=", round(expected, 1),
                       ", FDR=", signif(FDR, 2))),
    hjust = -0.1,
    size = 3
  ) +
  theme_classic(base_size = 12) +
  xlab("Fold enrichment of overlap (at least these sets)") +
  ylab("") +
  coord_cartesian(xlim = c(0, max(intersection_results$fold_enrichment) * 1.4))



plot_df <- intersection_results %>%
  mutate(
    intersection_short = c(
      "LowDisp + MuscleDown",
      "LowDisp + NeuronDown",
      "LowDisp + Energy",
      "LowDisp + MuscleDown + Energy",
      "LowDisp + NeuronDown + Energy",
      "LowDisp + MuscleDown + NeuronDown",
      "All four"
    ),
    label = paste0(
      "obs=", observed,
      ", exp=", signif(expected, 2),
      ", FDR=", signif(FDR, 2)
    ),
    intersection_short = factor(
      intersection_short,
      levels = intersection_short[order(fold_enrichment)]
    )
  )

ggplot(plot_df, aes(x = fold_enrichment, y = intersection_short)) +
  geom_col(fill = "gray35", width = 0.7) +
  geom_text(
    aes(label = label),
    hjust = -0.05,
    size = 3.2
  ) +
  scale_x_log10() +
  theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    plot.margin = margin(10, 120, 10, 10)
  ) +
  xlab("Fold enrichment of overlap\n(at least these sets, log10 scale)") +
  ggtitle("Convergence of lower-dispersion, aging-down, and energy-metabolism genes") +
  coord_cartesian(
    xlim = c(min(plot_df$fold_enrichment) * 0.9,
             max(plot_df$fold_enrichment) * 2)
  )


FourWayOverlapGenes<-upset_df_energy%>%
  mutate(Sum=LowerDispersion+Muscle_Age_Down+Neuron_Age_Down+Energy_GO)%>%
  filter(Sum==4)

FourWayOverlap<-FourWayOverlapGenes$Gene



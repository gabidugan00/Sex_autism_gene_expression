
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(stringr)
  library(cowplot)
  library(ggplot2)
  library(readxl)
  library(tibble)
  library(WriteXLS)
  library(ggVennDiagram)
  library(SuperExactTest)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(clusterProfiler)
})

# ============================================================
# Load data
# ============================================================

sex_dge <- readRDS("/vf/users/HGBomics/GabiDugan/OxCam/Human.Transcriptomics/Sex_ASD/Full_CT_BR/rCCA/Redo_114/sex_dge")

asd_dge <- read_xlsx("~/Data_S1_CELLTYPE-SPECIFIC_RAW_DEG.xlsx", col_names = TRUE)
asd_dge  <- asd_dge[!grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", asd_dge$gene), ]  # remove date-like gene names

scz_sheets <- excel_sheets("~/Copy_SCZ_Suptable4_DGE_results.xlsx")
scz_dge    <- lapply(scz_sheets, function(s) {
  df <- read_xlsx("~/Copy_SCZ_Suptable4_DGE_results.xlsx", col_names = TRUE, sheet = s)
  df$Celltype <- s
  df
})
SCZ_dge <- bind_rows(scz_dge)


# ============================================================
# Define DEG sets
# ============================================================

ASD_DEGs <- asd_dge |>
  filter(abs(logFC_ASD_vs_CTL) > 0.2, FDR_pVal_CTL_vs_ASD < 0.05) |>
  distinct(gene, .keep_all = TRUE)

SEX_DEGs <- sex_dge |>
  filter(adj.P.Val < 0.1) |>
  distinct(gene, .keep_all = TRUE)

SCZ_DEGS <- SCZ_dge |>
  filter(abs(Meta_logFC) > 0.1, Meta_adj.P.Val < 0.05) |>
  distinct(gene, .keep_all = TRUE)



# ============================================================
# Panel 1: Venn Diagram
# ============================================================

gene_lists <- list(Sex = SEX_DEGs$gene, Autism = ASD_DEGs$gene, Schizophrenia = SCZ_DEGS$gene)

p1 = ggVennDiagram(gene_lists) + theme(legend.position = "none")


# ============================================================
# Fisher's exact tests: sex vs. condition (all genes + autosomal only)
# ============================================================

run_fisher <- function(dge1, dge2, deg1, deg2) {
  bg  <- intersect(dge1$gene, dge2$gene)
  tbl <- table(
    DE_1 = bg %in% deg1$gene,
    DE_2 = bg %in% deg2$gene
  )
  fisher.test(tbl, alternative = "two.sided")
}

fisher_results <- list(
  "Sex & Autism"               = run_fisher(asd_dge,      sex_dge,      ASD_DEGs,      SEX_DEGs),
  "Sex & Schizophrenia"        = run_fisher(SCZ_dge,      sex_dge,      SCZ_DEGS,      SEX_DEGs)
  # "Autism & Schizophrenia"    = run_fisher(asd_dge,      sex_dge,      ASD_DEGs,      SEX_DEGs),
  #   "Auto Sex & Autism"          = run_fisher(asd_dge_auto, sex_dge_auto, ASD_DEGs_auto, SEX_DEGs_auto),
  #   "Auto Sex & Schizophrenia"   = run_fisher(SCZ_dge_auto, sex_dge_auto, SCZ_DEGS_auto, SEX_DEGs_auto)
)

fisher_df <- bind_rows(
  lapply(names(fisher_results), function(nm) {
    data.frame(OR = fisher_results[[nm]]$estimate, p = fisher_results[[nm]]$p.value, comparison = nm)
  })
) |>
  mutate(
    padj       = p.adjust(p,method="fdr"),
    comparison = factor(comparison, levels = rev(c(
      "Sex & Autism", 
      "Sex & Schizophrenia"
    )))
  )

# ============================================================
# Panel 2: odds ratio dotplot
# ============================================================

venn_colors <- c(
  "Sex & Autism"             = "#F0C674",
  "Sex & Schizophrenia"      = "#D1EEEE"
)

fisher_df$significant = "Yes"
fisher_df$significant[fisher_df$FDR > 0.05] = "No"

fisher_df$comparison = factor(fisher_df$comparison, levels =c("Sex & Autism", "Sex & Schizophrenia"))

p2_or_dotplot = ggplot(fisher_df, aes(x = comparison, y = OR, fill = comparison)) +
  geom_hline(yintercept = 1, linetype = "dashed",color = "grey50", linewidth = 0.4) +
  geom_segment( aes(x=comparison, xend=comparison, y=0, yend=OR)) +
  geom_point(aes(size = -log10(padj), color = significant, alpha = significant),shape =21, stroke = 1.5) +
  scale_color_manual(values = c("Yes" = "black", "No" = "white")) +
  scale_alpha_manual(values = c("Yes" = 1, "No" = 0.5)) +
  scale_size_continuous(name   = "-log10(FDR q)", range  = c(5, 30),limits=c(0,100)) +
  scale_fill_manual(name   = "Shared DEGs",values = venn_colors) +
  labs(x = "Comparison",y = "Odds Ratio",title    = "") +
  theme_cowplot()+
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 18), 
        axis.text.y      = element_text(size = 18), 
        axis.title      = element_text(size = 20), 
        plot.margin = margin(t = 10, r = 10, b = 20, l = 40, unit = "pt"),
        panel.grid.major.x = element_line(color = "grey92"),
        panel.grid.minor   = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.position  = "right")+
  ylim(0,5)

p2_or_dotplot



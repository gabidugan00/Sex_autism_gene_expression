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

sex_dge <- readRDS("sex_dge") #done in sex- limma voom analysis

asd_dge <- read_xlsx("~/Data_S1_CELLTYPE-SPECIFIC_RAW_DEG.xlsx", col_names = TRUE) #read in from Wamsley et al. DGE
asd_dge  <- asd_dge[!grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", asd_dge$gene), ]  # remove date-like gene names

scz_sheets <- excel_sheets("~/Copy_SCZ_Suptable4_DGE_results.xlsx") #read in from Ruzicka et al. DGE
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

fisher_results <- list(
  "Sex & Autism"               = run_fisher(asd_dge,      sex_dge,      ASD_DEGs,      SEX_DEGs),
  "Sex & Schizophrenia"        = run_fisher(SCZ_dge,      sex_dge,      SCZ_DEGS,      SEX_DEGs)
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
  "Auto Sex & Autism"        = "#F0C674",
  "Sex & Schizophrenia"      = "#D1EEEE",
  "Auto Sex & Schizophrenia" = "#D1EEEE"
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

# ============================================================
# Panel 3: Chromosome of origin for overlapping DEGs
# ============================================================

autism_sex_degs        <- data.frame(gene = intersect(ASD_DEGs$gene, SEX_DEGs$gene)) |> annotate_chromosome()
schizophrenia_sex_degs <- data.frame(gene = intersect(SCZ_DEGS$gene, SEX_DEGs$gene)) |> annotate_chromosome()

chromosome_degs <- bind_rows(
  mutate(autism_sex_degs,        Dataset = "Sex & Autism"),
  mutate(schizophrenia_sex_degs, Dataset = "Sex & Schizophrenia")
)

chrom_colors <- c(Auto = "grey40", "X-Linked" = "#ddac8e", "Y-Linked" = "#b35806", PAR = "burlywood3")

p3_barplot_degs_chromosome <- ggplot(chromosome_degs, aes(x = 1, fill = chromosome, y = Dataset)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = chrom_colors) +
  xlim(0, 600) +
  coord_flip() +
  theme_cowplot() +
  theme(
    axis.text.x   = element_text(size = 18, angle =45, hjust=1),
    axis.text.y   = element_text(size = 18),
    axis.title  = element_text(size = 20)
  ) +
  xlab("Overlapping DEG Count") +
  ylab("Comparison")

# Triple overlap
bg_all     <- Reduce(intersect, list(asd_dge$gene, sex_dge$gene, SCZ_dge$gene))

shared_all <- Reduce(intersect, list(ASD_DEGs$gene, SEX_DEGs$gene, SCZ_DEGS$gene))

tbl_triple <- table(
  DE_all = bg_all %in% shared_all,
  SC     = bg_all %in% xy_genes_nopar
)

fisher.test(tbl_triple)

# ============================================================
# Panel 3: Enrichment of Sex-chromosome genes 
# ============================================================

fisher_results_sc <- list(
  "X-Chromosome Sex & Autism"               = run_sc_fisher(ASD_DEGs, SEX_DEGs, asd_dge,  sex_dge, x_genes_nopar),
  "X-Chromosome Sex & Schizophrenia"        = run_sc_fisher(SCZ_DEGS, SEX_DEGs, SCZ_dge,  sex_dge,x_genes_nopar),
  "Y-Chromosome Sex & Autism"               = run_sc_fisher(ASD_DEGs, SEX_DEGs, asd_dge,  sex_dge,y_genes_nopar),
  "Y-Chromosome Sex & Schizophrenia"        = run_sc_fisher(SCZ_DEGS, SEX_DEGs, SCZ_dge,  sex_dge,y_genes_nopar),
  "PAR Sex & Autism"               = run_sc_fisher(ASD_DEGs, SEX_DEGs, asd_dge,  sex_dge,par_genes),
  "PAR Sex & Schizophrenia"        = run_sc_fisher(SCZ_DEGS, SEX_DEGs, SCZ_dge,  sex_dge,par_genes),
  "Sex-Chromosomes Sex & Autism" = run_sc_fisher(ASD_DEGs, SEX_DEGs, asd_dge,  sex_dge,xy_genes_nopar),
  "Sex-Chromosomes Sex & Schizophrenia" = run_sc_fisher(SCZ_DEGS, SEX_DEGs, SCZ_dge,  sex_dge,xy_genes_nopar)
)

fisher_df_sc <- bind_rows(
  lapply(names(fisher_results_sc), function(nm) {
    data.frame(OR = fisher_results_sc[[nm]]$estimate, p = fisher_results_sc[[nm]]$p.value, comparison = nm, datasets = ifelse(grepl("Autism",nm ),"Sex & Autism","Sex & Schizophrenia"), sex_chromosome = sub(" Sex &.*", "",nm))
  })
) |>
  mutate(
    padj       = p.adjust(p, method="fdr"),
    datasets = factor(datasets, levels = rev(c(
      "Sex & Autism", 
      "Sex & Schizophrenia"
    )))
  )

fisher_df_sc$significant = "Yes"
fisher_df_sc$significant[fisher_df_sc$padj > 0.05] = "No"

fisher_df_sc$datasets = factor(fisher_df_sc$datasets, levels =c("Sex & Autism", "Sex & Schizophrenia"))
fisher_df_sc$sex_chromosome = factor(fisher_df_sc$sex_chromosome, levels =c("X-Chromosome", "Y-Chromosome","PAR","Sex-Chromosomes"))

fisher_sc_figure_df = fisher_df_sc[fisher_df_sc$sex_chromosome =="Sex-Chromosomes",]

p4_sex_chromosome_enrichment = ggplot(fisher_sc_figure_df, aes(x = datasets, y = OR, fill = datasets)) +
  geom_hline(yintercept = 1, linetype = "dashed",color = "grey50", linewidth = 0.4) +
  geom_segment( aes(x=datasets, xend=datasets, y=0, yend=OR, color = significant)) +
  geom_point(aes(size = -log10(padj), color = significant, alpha = significant),shape =21, stroke = 0.6) +
  scale_color_manual(values = c("Yes" = "black", "No" = "white")) +
  scale_alpha_manual(values = c("Yes" = 1, "No" = 0.9)) +
  scale_size_continuous(name   = "-log10(FDR q)", range  = c(5, 30),limits=c(0,100)) +
  scale_fill_manual(name   = "Shared DEGs",values = venn_colors) +
  labs(x = "Comparison",y = "Odds Ratio",title    = "Sex-Chromosome Genes") +
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


# ============================================================
# Sup Analysis: SuperExactTest: triple overlap significance
# ============================================================

gene_lists <- list(Sex = SEX_DEGs$gene, Autism = ASD_DEGs$gene, Schizophrenia = SCZ_DEGS$gene)

gene_lists_bg <- lapply(gene_lists, function(x) x[x %in% bg_all])
overlap_test  <- supertest(gene_lists_bg, n = length(unique(bg_all)))
overlap_test$padj = p.adjust(overlap_test$P.value,method = "fdr")

summary(overlap_test)

# ============================================================
# Sup analysis: GO enrichment: pairwise overlaps
# ============================================================

aut_sex_enrich <- run_enrichGO(ASD_DEGs, SEX_DEGs, asd_dge, sex_dge)

# ============================================================
# Panel 5 GO enrichment: plots
# ============================================================
p5_sex_aut_enrich = plot_enrichGO(aut_sex_enrich$Full_Enrich, title="Sex & Autism", fill_color ="#F0C674" )

##Alternative optin if i want fill to follow a gradient
# fill = scale_fill_gradient2(low = 'white', mid = '#f7f7f7', high = '#F0C674', midpoint = 0,limits= c(0,12)

p5_sex_aut_enrich

# =============================================================================
# Panel 6 Enrichment of risk gene sets on global union of DEGs
# =============================================================================

# Autism and sex
aut_sex_shared <- intersect(unique(ASD_DEGs$gene), unique(SEX_DEGs$gene))
bg <- intersect(unique(asd_dge$gene), unique(sex_dge$gene))

riskgene_enrichment_autism_sex <- run_fisher_enrichment(aut_sex_shared, bg, riskgene_pathways)

riskgene_enrichment_autism_sex$FDR <- p.adjust(riskgene_enrichment_autism_sex$P_Value, method = "fdr")

pathway_order = str_wrap(c( "Combined Autism","Autism - SFARI Score 1","Autism - SFARI Score 2",
                            "Autism - Fu et al.",  "Autism - Sattestrom et al."), width =20)

riskgene_enrichment_autism_sex$pathway = str_wrap(riskgene_enrichment_autism_sex$Pathway, width =20)

riskgene_enrichment_autism_sex  = riskgene_enrichment_autism_sex[riskgene_enrichment_autism_sex$pathway  %in% pathway_order,]
riskgene_enrichment_autism_sex$pathway <- factor(riskgene_enrichment_autism_sex$pathway, levels = pathway_order)
riskgene_enrichment_autism_sex$comparison = "Sex & Autism"
riskgene_enrichment_autism_sex$significant = riskgene_enrichment_autism_sex$FDR < 0.05 

p_risk_gene_sex_autism <- ggplot(riskgene_enrichment_autism_sex, aes(x = pathway, y = Odds_Ratio, fill = comparison)) +
  geom_hline(yintercept = 1, linetype = "dashed",color = "grey50", linewidth = 0.4) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=Odds_Ratio, color = significant)) +
  geom_line(aes(group = pathway), color = "grey80", linewidth = 0.4) +
  geom_point(aes(size = -log10(FDR), color = significant, alpha = significant),shape =21, stroke = 1.5) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) +
  scale_size_continuous(name   = "-log10(FDR q)",range  = c(5, 30),limits=c(0,100)) +
  scale_fill_manual(name   = "Shared DEGs",values = venn_colors) +
  labs(x = "Risk Gene List",
       y = "Odds Ratio",
       title    = "") +
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

p_risk_gene_sex_autism

# ============================================================
# Save main text figure
# ============================================================

pdf("MainText_Figure1_global_DEGs_dropSCZ_redo.pdf", width = 20, height = 16)

plot_grid(
  
  plot_grid(p1, p2_or_dotplot + theme(legend.position = "bottom"),
            p3_barplot_degs_chromosome + theme(legend.position = "bottom"),
            p4_sex_chromosome_enrichment+ theme(legend.position = "bottom"),
            align = "h", axis="bt",
            ncol = 4, rel_widths = c(1.2,0.8,0.8,0.8)
  ),
  plot_grid(p5_sex_aut_enrich+ theme(legend.position = "none"), p_risk_gene_sex_autism, ncol = 2, align = "h", axis="bt", rel_widths=c(1.2,1)), nrow=2, rel_heights =c(1.3,1)
)
dev.off()


# ============================================================
# GO enrichment: pairwise + unique overlaps
# ============================================================

# GO enrichment: excluding overlapping DEGs
unique_enrich_aut_sex <- run_enrichGO(ASD_DEGs, SEX_DEGs, asd_dge, sex_dge, remove_genes = SCZ_DEGS$gene)
# "before removal"
# [1] 599
# [1] "after removal"
# [1] 324
unique_enrich_scz_sex <- run_enrichGO(SCZ_DEGS, SEX_DEGs, SCZ_dge, sex_dge, remove_genes = ASD_DEGS$gene)
# [1] "before removal"
# [1] 327
# [1] "after removal"
# [1] 52

plot_enrichGO(unique_enrich_aut_sex$Full_Enrich, title="Unique Sex & Autism" , fill = scale_fill_gradient2(low = 'white', mid = '#f7f7f7', high = 'skyblue', midpoint = 0, limits= c(0,12)))

plot_enrichGO(unique_enrich_scz_sex$Full_Enrich, title="Unique Sex & Schizophrenia" , fill = scale_fill_gradient2(low = 'white', mid = '#f7f7f7', high = 'skyblue', midpoint = 0, limits= c(0,12)))


# Unique pathways per overlap 
unique_pathways <- function(target, other1, other2) {
  all_other <- union(other1$Sig_Enrich$Description, other2$Sig_Enrich$Description)
  target$Sig_Enrich[!target$Sig_Enrich$Description %in% all_other, ]
}

unique_sex_aut <- unique_pathways(aut_sex_enrich, aut_scz_enrich, scz_sex_enrich)
unique_sex_scz <- unique_pathways(scz_sex_enrich, aut_scz_enrich, aut_sex_enrich)
unique_aut_scz <- unique_pathways(aut_scz_enrich, aut_sex_enrich, scz_sex_enrich)



# ============================================================
# GO enrichment: global overlap
# ============================================================

# GO enrichment for the triple overlap (275 shared DEGs)
triple_enrich <- enrichGO_function(
  data.frame(gene = shared_all), data.frame(gene = shared_all),
  data.frame(gene = bg_all),     data.frame(gene = bg_all)
)

WriteXLS(triple_enrich$Sig_Enrich, "July2026_Overlapping275DEGs_sig_enrich.xlsx")


# SH3 domain binding genes in unique autism-sex overlap
sh3_row   <- unique_enrich_aut_sex$Sig_Enrich[unique_enrich_aut_sex$Sig_Enrich$Description == "SH3 domain binding", ]
sh3_entrez <- unlist(str_split(sh3_row$geneID, "/"))
sh3_genes  <- mapIds(org.Hs.eg.db, keys = sh3_entrez, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")


# ============================================================
# Comparing odds ratios
# ============================================================
odds_ratio = lapply(names(fisher_results), function(nm){
  res = fisher_results[[nm]]
  data.frame(
    OR = as.numeric(res$estimate),
    Upper = as.numeric(res$conf.int[2]),
    Lower = as.numeric(res$conf.int[1]),
    Comparison = nm
  )
})

odds_ratio = Reduce(rbind, odds_ratio)

odds_ratio$logOR = log(odds_ratio$OR)

odds_ratio$SE_logOR = (log(odds_ratio$Upper) - log(odds_ratio$Lower)) / (2 * 1.96)

logOR_diff = odds_ratio$logOR[1] - odds_ratio$logOR[2]
SE_diff = sqrt(odds_ratio$SE_logOR[1]^2 + odds_ratio$SE_logOR[2]^2)

z.score = logOR_diff / SE_diff
p.value = 2 * (1 - pnorm(abs(z.score)))




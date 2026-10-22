# ============================================================
# GO enrichment: pairwise overlaps
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


p_risk_gene_sex_autism


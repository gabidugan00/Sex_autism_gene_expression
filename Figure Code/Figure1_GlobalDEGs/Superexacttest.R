
# ============================================================
# SuperExactTest: triple overlap significance
# ============================================================

gene_lists <- list(Sex = SEX_DEGs$gene, Autism = ASD_DEGs$gene, Schizophrenia = SCZ_DEGS$gene)

gene_lists_bg <- lapply(gene_lists, function(x) x[x %in% bg_all])
overlap_test  <- supertest(gene_lists_bg, n = length(unique(bg_all)))
overlap_test$padj = p.adjust(overlap_test$P.value,method = "fdr")

summary(overlap_test)

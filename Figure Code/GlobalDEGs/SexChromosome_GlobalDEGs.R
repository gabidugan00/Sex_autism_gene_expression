# ============================================================
# Chromosome of origin for overlapping DEGs
# ============================================================

annotate_chromosome <- function(df) {
  df |>
    mutate(chromosome = case_when(
      gene %in% y_genes_nopar ~ "Y-Linked",
      gene %in% x_genes_nopar ~ "X-Linked",
      gene %in% par_genes     ~ "PAR",
      TRUE                    ~ "Auto"
    ))
}

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

p3_barplot_degs_chromosome

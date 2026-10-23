
#Functions
filter_DEGs <- function(df, logFC_col, p_col, logFC_cut = 0.2, p_cut = 0.05, remove_xy = FALSE, xy_genes = NULL) {
  out <- df %>%
    dplyr::filter(abs(.data[[logFC_col]]) > logFC_cut, .data[[p_col]] < p_cut) %>%
    dplyr::filter(!grepl("^(RPL|RPS|MRPL|MRPS)", gene), !grepl("^MT-", gene))
  if (remove_xy && !is.null(xy_genes)) out <- out %>% dplyr::filter(!gene %in% xy_genes)
  out
}

run_fishers_conconcate = function(asd_dge, sex_dge, looser = FALSE, remove_xy = FALSE, xy_genes = NULL, asd_thresh = list(logFC = 0.2, p = 0.05)) {
  
  sex_dge$Sex_Broad_Celltype_Br = paste(sex_dge$Sex_Broad_Celltype, sex_dge$brainregion_shorten, sep=".")
  
  ASD_DEGs <- filter_DEGs(asd_dge, "logFC_ASD_vs_CTL", "FDR_pVal_ASD_vs_CTL",
                          logFC_cut = asd_thresh$logFC, p_cut = asd_thresh$p,
                          remove_xy = remove_xy, xy_genes = xy_genes)
  
  SEX_DEGs <- filter_DEGs(sex_dge, "logFC", "adj.P.Val",
                          logFC_cut = ifelse(looser, 0, 0.2),
                          p_cut = ifelse(looser, 0.1, 0.05),
                          remove_xy = remove_xy, xy_genes = xy_genes)
  
  ASD_split <- split(asd_dge, asd_dge$ASD_Broad_Celltype)
  SEX_split <- split(sex_dge, sex_dge$Sex_Broad_Celltype_Br)
  ASD_DEG_split <- split(ASD_DEGs, ASD_DEGs$ASD_Broad_Celltype)
  SEX_DEG_split <- split(SEX_DEGs, SEX_DEGs$Sex_Broad_Celltype_Br)
  
  results <- expand.grid(ASD = names(ASD_DEG_split), Sex = names(SEX_DEG_split)) %>%
    pmap_dfr(function(ASD, Sex) {
      ASD_genes <- unique(ASD_DEG_split[[ASD]]$gene)
      SEX_genes <- unique(SEX_DEG_split[[Sex]]$gene)
      background <- intersect(unique(ASD_split[[ASD]]$gene), unique(SEX_split[[Sex]]$gene))
      
      if (length(ASD_genes) == 0 || length(SEX_genes) == 0 || length(background) == 0) {
        return(tibble(
          ASD_celltype = ASD, Sex_celltype = Sex, overlap_genes = 0,
          total_ASD_DEGs = length(ASD_genes), total_Sex_DEGs = length(SEX_genes),
          total_background = length(background), odds_ratio = NA, p_value = NA
        ))
      }
      
      df <- tibble(Gene = background,
                   DE_ASD = background %in% ASD_genes,
                   DE_Sex = background %in% SEX_genes)
      
      tbl <- table(df$DE_Sex, df$DE_ASD)
      
      if (nrow(tbl) < 2 || ncol(tbl) < 2) {
        return(tibble(
          ASD_celltype = ASD, 
          Sex_celltype = Sex,
          overlap_genes = length(intersect(ASD_genes, SEX_genes)),
          total_ASD_DEGs = length(ASD_genes), 
          total_Sex_DEGs = length(SEX_genes),
          total_background = length(background), 
          odds_ratio = NA, 
          p_value = NA
        ))
      }
      
      test <- fisher.test(tbl)
      
      tibble(
        ASD_celltype = ASD, 
        Sex_celltype = Sex,
        overlap_genes = length(intersect(ASD_genes, SEX_genes)),
        total_ASD_DEGs = length(ASD_genes), 
        total_Sex_DEGs = length(SEX_genes),
        total_background = length(background),
        odds_ratio = test$estimate, 
        p_value = test$p.value
      )
    })
  
  results %>% mutate(FDR = p.adjust(p_value, "fdr"), Sig = FDR < 0.05)
}

#W permutations
fisher_list = list()
run_perm_fishers <- function(asd_dge, perm_list, ...) {
  
  for (i in seq_along(perm_list)) {
    message("Processing permutation ", i, " of ", length(perm_list), "...")
    
    df <- perm_list[[i]]
    
    df$brainregion_celltype <- df$Celltype_BrainRegion
    df$Sex_Brainregion <- sub(".*\\.", "", df$brainregion_celltype)
    df$Sex_Celltype    <- sub("\\..*", "", df$brainregion_celltype)
    
    df$brainregion_shorten <- dplyr::recode(
      df$Sex_Brainregion,
      "Angular_Gyrus"="AngGyr",
      "Caudal_Insula"= "CauIns",
      "Fusiform_Gyrus"="FusGyr",
      "Inferior_Lateral_TemporalOccipital"="InfLTO",
      "Retrosplenial_cortex"= "RetCor",
      "Intraparietal_sulcus" = "IntParSul")
    
    df <- df %>%
      dplyr::mutate(Sex_Broad_Celltype = case_when(
        Sex_Celltype %in% c("L23IT","L4IT","L5IT","L56NP","L6IT","L6b","L6CT",
                            "L5ET/L6IT","L6CT/L5IT","L23IT/L5IT") ~ "ExcNeu",
        Sex_Celltype %in% c("VIP","PVALB","SST","LAMP5","ADARB2","PAX6",
                            "Chandelier","Sncg") ~ "InhNeu",
        Sex_Celltype == "Astro" ~ "Astro",
        Sex_Celltype == "Micro" ~ "Micro",
        Sex_Celltype %in% c("Oligo","OPC") ~ "Oligo/OPC",
        Sex_Celltype == "BEC" ~ "BEC_Endo",
        TRUE ~ "Other"
      ))
    
    # Run fisher's test as above, but this produces Inf sometimes, since there are some perfect overlaps of low sex DEGs with autism DEGs 
    fisher_res <- run_fishers_conconcate(asd_dge, df, looser = TRUE)
    fisher_res$permutation_number <- i
    fisher_list[[i]] <- fisher_res
  }
  names(fisher_list) <- paste0("Permutation_", seq_along(perm_list))
  return(fisher_list)
}

add_columns_perm_fishers <- function(fisher_list) {
  lapply(fisher_list, function(df) {
    df$Sex_Brainregion <- sub(".*\\.", "", df$Sex_celltype)
    df$Sex_Celltype <- sub("\\..*", "", df$Sex_celltype)
    
    df$Sex_Brainregion_short <- dplyr::recode(
      df$Sex_Brainregion,
      "Angular_Gyrus"="AngGyr",
      "Caudal_Insula"= "CauIns",
      "Fusiform_Gyrus"="FusGyr",
      "Inferior_Lateral_TemporalOccipital"="InfLTO",
      "Retrosplenial_cortex"= "RetCor",
      "Intraparietal_sulcus" = "IntParSul"
    )
    df
  })
}

run_perm_fishers_NOY <- function(asd_dge, perm_list, gene_list_toexclude, ...) {
  
  for (i in seq_along(perm_list)) {
    message("Processing permutation ", i, " of ", length(perm_list), "...")
    
    df <- perm_list[[i]]
    
    df$brainregion_celltype <- df$Celltype_BrainRegion
    df$Sex_Brainregion <- sub(".*\\.", "", df$brainregion_celltype)
    df$Sex_Celltype    <- sub("\\..*", "", df$brainregion_celltype)
    
    df$brainregion_shorten <- dplyr::recode(
      df$Sex_Brainregion,
      "Angular_Gyrus"="AngGyr",
      "Caudal_Insula"= "CauIns",
      "Fusiform_Gyrus"="FusGyr",
      "Inferior_Lateral_TemporalOccipital"="InfLTO",
      "Retrosplenial_cortex"= "RetCor",
      "Intraparietal_sulcus" = "IntParSul")
    
    df <- df %>%
      dplyr::mutate(Sex_Broad_Celltype = case_when(
        Sex_Celltype %in% c("L23IT","L4IT","L5IT","L56NP","L6IT","L6b","L6CT",
                            "L5ET/L6IT","L6CT/L5IT","L23IT/L5IT") ~ "ExcNeu",
        Sex_Celltype %in% c("VIP","PVALB","SST","LAMP5","ADARB2","PAX6",
                            "Chandelier","Sncg") ~ "InhNeu",
        Sex_Celltype == "Astro" ~ "Astro",
        Sex_Celltype == "Micro" ~ "Micro",
        Sex_Celltype %in% c("Oligo","OPC") ~ "Oligo/OPC",
        Sex_Celltype == "BEC" ~ "BEC_Endo",
        TRUE ~ "Other"
      ))
    
    df_filt = df %>%
      dplyr::filter(!gene %in% gene_list_toexclude)
    
    asd_dge_filt = asd_dge %>%
      dplyr::filter(!gene %in% gene_list_toexclude)
    
    # Run fisher's test as above, but this produces Inf sometimes, since there are some perfect overlaps of low sex DEGs with autism DEGs 
    fisher_res <- run_fishers_conconcate(asd_dge_filt, df_filt, looser = TRUE)
    
    fisher_res$permutation_number <- i
    
    fisher_list[[i]] <- fisher_res
  }
  
  names(fisher_list) <- paste0("Permutation_", seq_along(perm_list))
  
  return(fisher_list)
  
}

add_columns_perm_fishers <- function(fisher_list) {
  lapply(fisher_list, function(df) {
    df$Sex_Brainregion <- sub(".*\\.", "", df$Sex_celltype)
    df$Sex_Celltype <- sub("\\..*", "", df$Sex_celltype)
    
    df$Sex_Brainregion_short <- dplyr::recode(
      df$Sex_Brainregion,
      "Angular_Gyrus"="AngGyr",
      "Caudal_Insula"= "CauIns",
      "Fusiform_Gyrus"="FusGyr",
      "Inferior_Lateral_TemporalOccipital"="InfLTO",
      "Retrosplenial_cortex"= "RetCor",
      "Intraparietal_sulcus" = "IntParSul"
    )
    df
  })
}


# ============================================================
# DEG Directionality Pipeline — Self-contained functions
# ============================================================
map_broad_celltype <- function(df, celltype_col) {
  df %>%
    dplyr::mutate(Broad_Celltype = case_when(
      !!sym(celltype_col) %in% c("EXT_6_L23", "EXT_3_L23", "EXT_11_L34", "EXT_10_L34", 
                                 "EXT_1_L4", "EXT_13_L45", "EXT_14_L5B", "EXT_4_L56", 
                                 "EXT_7_L2356", "EXT_5_L56", "EXT_2_L56", "EXT_12_L56", 
                                 "EXT_8_L6", "EXT_9_L6", "Ex-L2", "Ex-L23", "Ex-L3", 
                                 "Ex-L4_MYLK", "Ex-L45_MET", "Ex-L45_LRRK1", "Ex-L5b_HTR2C", 
                                 "Ex-L56", "Ex-L56_CC_NTNG2", "Ex-L6_CC_SEMA3A", "Ex-L6b_SEMA3D", 
                                 "Ex-L6b_SEMA3E", "L23IT", "L4IT", "L5IT", "L56NP", "L6IT", 
                                 "L6b", "L6CT", "L5ET/ L6IT", "L6CT / L5IT", "L23IT/ L5IT") ~ "ExcNeu",
      !!sym(celltype_col) %in% c("INT_5_SST", "INT_2_PV_SULF", "INT_6_PV", "INT_1_KIT", 
                                 "INT_3_RELN_ND", "INT_4_VIP_CR", "In-Rosehip_CHST9", 
                                 "In-Rosehip_TRPC3", "In-Reelin", "In-VIP", "In-PV_Chandelier", 
                                 "In-PV_Basket", "In-SST", "VIP", "PVALB", "SST", "LAMP5", 
                                 "ADARB2", "PAX6", "Chandelier", "Sncg") ~ "InhNeu",
      !!sym(celltype_col) %in% c("ASTRO_1", "ASTRO_2", "ASTRO_3", "Ast", "Astro") ~ "Astro",
      !!sym(celltype_col) %in% c("MG1", "MG2", "BBB_MG3", "Mic", "Micro") ~ "Micro",
      !!sym(celltype_col) %in% c("ODC_1", "ODC_2", "ODC_3", "OPCS_1", "OPCS_2", "OPCS_3", 
                                 "OPCS_4", "Oli", "OPC", "Oligo") ~ "Oligo/OPC",
      !!sym(celltype_col) %in% c("BBB_Endo_1", "BBB_Pericytes", "BBB_SM(Mural)", "Endo", 
                                 "Pericytes", "BEC") ~ "NeuVasc",
      TRUE ~ "Other"
    ))
}

.filter_sex_degs <- function(sex_dge,
                             logFC_cut = 0,
                             p_cut     = 0.1) {
  stopifnot(
    all(c("gene", "logFC", "adj.P.Val",
          "Sex_Broad_Celltype", "brainregion_shorten") %in% names(sex_dge))
  )
  sex_dge %>%
    filter(abs(logFC) > logFC_cut, adj.P.Val < p_cut)
}

.filter_asd_degs <- function(asd_dge,
                             logFC_cut = 0.2,
                             p_cut     = 0.05) {
  stopifnot(
    all(c("gene", "logFC_ASD_vs_CTL", "FDR_pVal_ASD_vs_CTL",
          "ASD_Broad_Celltype") %in% names(asd_dge))
  )
  asd_dge %>%
    filter(abs(logFC_ASD_vs_CTL) > logFC_cut, FDR_pVal_ASD_vs_CTL < p_cut)
}

.get_gene_overlaps <- function(sex_degs, asd_degs) {
  
  asd_split <- split(asd_degs, asd_degs$ASD_Broad_Celltype)
  
  sex_degs  <- sex_degs %>%
    mutate(Sex_CT_BR = paste(Sex_Broad_Celltype, brainregion_shorten, sep = "."))
  sex_split <- split(sex_degs, sex_degs$Sex_CT_BR)
  
  results <- expand.grid(
    ASD = names(asd_split),
    Sex = names(sex_split),
    stringsAsFactors = FALSE
  )
  
  overlaps <- Map(function(a, s) {
    genes <- intersect(
      unique(sex_split[[s]]$gene),
      unique(asd_split[[a]]$gene)
    )
    list(genes = genes, n = length(genes))
  }, results$ASD, results$Sex)
  
  results$IntersectingGenes <- lapply(overlaps, `[[`, "genes")
  results$NumberOverlap     <- vapply(overlaps, `[[`, integer(1), "n")
  results$Joint_Celltype    <- paste(results$ASD, results$Sex, sep = "_")
  
  results
}

.compute_mean_logfc <- function(sex_dge, asd_dge) {
  
  mean_sex <- sex_dge %>%
    dplyr::group_by(Sex_Broad_Celltype, brainregion_shorten, gene) %>%
    dplyr::summarise(mean_logFC_sex = mean(logFC), .groups = "drop") %>%
    dplyr::mutate(celltype = paste(Sex_Broad_Celltype, brainregion_shorten, sep = "."))
  
  mean_asd <- asd_dge %>%
    dplyr::group_by(ASD_Broad_Celltype, gene) %>%
    dplyr::summarise(mean_logFC_autism = mean(logFC_ASD_vs_CTL), .groups = "drop")
  
  mean_sex %>%
    inner_join(mean_asd, by = "gene") %>%
    mutate(Joint_Celltype = paste(ASD_Broad_Celltype, celltype, sep = "_"))
}

#   Assign directional quadrant per gene 
#   MB.AU = male-biased & upregulated in ASD  (concordant with FPE)
#   MB.AD = male-biased & downregulated in ASD
#   FB.AU = female-biased & upregulated in ASD
#   FB.AD = female-biased & downregulated in ASD  (concordant with FPE)

.assign_directions <- function(gene_overlaps, mean_logfc) {
  
  gene_overlaps %>%
    dplyr::filter(NumberOverlap > 0) %>%
    rowwise() %>%
    reframe(
      mean_logfc %>%
        dplyr::filter(
          gene %in% IntersectingGenes,
          Joint_Celltype == paste(ASD, Sex, sep = "_")
        )
    ) %>%
    dplyr::mutate(
      Sex_Dir    = if_else(mean_logFC_sex    > 0, "MB", "FB"),
      Autism_Dir = if_else(mean_logFC_autism > 0, "AU", "AD"),
      Joint_Dir  = paste(Sex_Dir, Autism_Dir, sep = ".")
    )
}


.build_counts <- function(directed_degs, empirical_pvals,
                          celltype_levels, brainregion_levels) {
  
  sig_combos <- empirical_pvals %>%
    dplyr::filter(p_empirical < 0.05) %>%
    dplyr::mutate(Joint_Celltype = paste(ASD_celltype, Sex_celltype, sep = "_"))
  
  directed_degs %>%
    dplyr::count(celltype, ASD_Broad_Celltype, Joint_Celltype, Joint_Dir) %>%
    dplyr::mutate(
      Align = if_else(Joint_Dir %in% c("FB.AU", "MB.AD"), "Neg", "Pos"),
      Joint_Dir = factor(
        Joint_Dir,
        levels = rev(c("MB.AU", "FB.AD", "FB.AU", "MB.AD"))),
      ASD_Broad_Celltype = factor(
        ASD_Broad_Celltype,
        levels = celltype_levels
      ),
      Sex_CT = factor(
        sub("\\..*", "", celltype),
        levels = celltype_levels
      ),
      Sex_BR = factor(
        sub(".*\\.", "", celltype),
        levels = brainregion_levels
      ),
      significant = if_else(Joint_Celltype %in% sig_combos$Joint_Celltype, "*", "")
    )
}

.run_directional_fisher <- function(df) {
  
  counts <- df %>%
    dplyr::mutate(
      sex_bias = if_else(grepl("^MB", Joint_Dir), "MB", "FB"),
      asd_dir  = if_else(grepl("AU$",  Joint_Dir), "AU", "AD")
    ) %>%
    dplyr::group_by(sex_bias, asd_dir) %>%
    dplyr::summarise(n = sum(n), .groups = "drop")
  
  keys <- paste(counts$sex_bias, counts$asd_dir, sep = ".")
  vals <- counts$n
  
  mat <- matrix(
    c(
      ifelse("MB.AU" %in% keys, vals[keys == "MB.AU"], 0),
      ifelse("FB.AU" %in% keys, vals[keys == "FB.AU"], 0),
      ifelse("MB.AD" %in% keys, vals[keys == "MB.AD"], 0),
      ifelse("FB.AD" %in% keys, vals[keys == "FB.AD"], 0)
    ),
    nrow = 2,
    dimnames = list(c("MB", "FB"), c("AU", "AD"))
  )
  
  ft <- fisher.test(mat)
  
  tibble(
    Joint_Celltype = unique(df$Joint_Celltype),
    MB.AU          = mat["MB", "AU"],
    MB.AD          = mat["MB", "AD"],
    FB.AU          = mat["FB", "AU"],
    FB.AD          = mat["FB", "AD"],
    total_overlap  = sum(mat),
    odds_ratio     = ft$estimate,
    OR_lower       = ft$conf.int[1],
    OR_upper       = ft$conf.int[2],
    p_value        = ft$p.value
  )
}

.run_all_fisher <- function(joint_dir_counts) {
  joint_dir_counts %>%
    filter(significant == "*") %>%
    group_by(Joint_Celltype) %>%
    group_modify(~ .run_directional_fisher(.x)) %>%
    ungroup() %>%
    mutate(
      p_adj           = p.adjust(p_value, method = "BH"),
      dir_significant = p_adj < 0.05
    )
}


.make_directionality_plot <- function(joint_dir_counts) {
  
  plot_data <- joint_dir_counts %>%
    mutate(n = if_else(significant == "*", n, 0L))
  
  if (nrow(plot_data) == 0)
    stop("No significant cell-type combinations to plot. Check empirical_pvals.")
  
  ggplot(plot_data, aes(y = n, x = "", fill = Joint_Dir, alpha = Align)) +
    geom_bar(position = "fill", stat = "identity") +
    coord_polar(theta = "y") +
    scale_fill_manual(
      values = c(
        "MB.AU" = "#8C510A",
        "MB.AD" = "#DFC27D",
        "FB.AD" = "#01665E",
        "FB.AU" = "#80CDC1"
      ),
      drop = FALSE
    ) +
    scale_alpha_manual(values = c("Pos" = 1, "Neg" = 0.5)) +
    facet_grid( Sex_CT~ Sex_BR + ASD_Broad_Celltype,switch = "y") +
    labs(
      x     = "Sex Broad Cell Type",
      y     = "Autism Broad Cell Type",
      fill  = "Joint Direction",
      alpha = "Alignment",
      title = "Directionality of Shared DEGs"
    ) +
    theme_void()+
    theme(
      legend.position = "top",
      strip.text      = element_text(size = 15),
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 25)
    )
}


run_deg_directionality <- function(
    sex_dge,
    asd_dge,
    empirical_pvals,
    sex_logFC_cut     = 0,
    sex_p_cut         = 0.1,
    asd_logFC_cut     = 0.2,
    asd_p_cut         = 0.05,
    celltype_levels   = c("ExcNeu", "InhNeu", "Astro", "Oligo/OPC", "Micro", "NeuVasc"),
    brainregion_levels = c("FusGyr", "AngGyr", "CauIns", "InfLTO", "IntParSul", "RetCor")
) {
  
  message("Step 1/6: Filtering DEGs...")
  sex_degs <- .filter_sex_degs(sex_dge, sex_logFC_cut, sex_p_cut)
  asd_degs <- .filter_asd_degs(asd_dge, asd_logFC_cut, asd_p_cut)
  message(sprintf("  %d sex DEGs | %d ASD DEGs", nrow(sex_degs), nrow(asd_degs)))
  
  message("Step 2/6: Computing gene overlaps...")
  gene_overlaps <- .get_gene_overlaps(sex_degs, asd_degs)
  n_with_overlap <- sum(gene_overlaps$NumberOverlap > 0)
  message(sprintf("  %d / %d cell-type combos have overlapping DEGs",
                  n_with_overlap, nrow(gene_overlaps)))
  
  message("Step 3/6: Computing mean logFC per gene × cell-type...")
  mean_logfc <- .compute_mean_logfc(sex_dge, asd_dge)
  
  message("Step 4/6: Assigning directional quadrants...")
  directed_degs <- .assign_directions(gene_overlaps, mean_logfc)
  
  message("Step 5/6: Building quadrant counts & flagging significant combos...")
  counts <- .build_counts(directed_degs, empirical_pvals,
                          celltype_levels, brainregion_levels)
  n_sig <- length(unique(counts$Joint_Celltype[counts$significant == "*"]))
  message(sprintf("  %d empirically significant combos flagged", n_sig))
  
  message("Step 6/6: Running directional Fisher tests & building plot...")
  fisher_results <- .run_all_fisher(counts)
  n_dir_sig <- sum(fisher_results$dir_significant, na.rm = TRUE)
  message(sprintf("  %d combos show significant directional bias (BH p < 0.05)",
                  n_dir_sig))
  
  plot <- .make_directionality_plot(counts)
  
  message("Done.")
  list(
    plot           = plot,
    fisher_results = fisher_results,
    counts         = counts,
    directed_degs  = directed_degs
  )
}


summarize_dge <- function(df, logfc_col, dataset_name, group_cols = c("Broad_Celltype", "gene")) {
  df %>%
    dplyr::group_by(across(all_of(group_cols))) %>%
    dplyr::summarise(
      mean_logFC = mean(!!sym(logfc_col), na.rm = TRUE),
      Dataset = dataset_name,
      .groups = "drop"
    )
}
#Functions
filter_DEGs <- function(df, logFC_col, p_col, logFC_cut = 0.2, p_cut = 0.05, remove_xy = FALSE, xy_genes = NULL) {
  out <- df %>%
    dplyr::filter(abs(.data[[logFC_col]]) > logFC_cut, .data[[p_col]] < p_cut) %>%
    dplyr::filter(!grepl("^(RPL|RPS|MRPL|MRPS)", gene), !grepl("^MT-", gene))
  if (remove_xy && !is.null(xy_genes)) out <- out %>% dplyr::filter(!gene %in% xy_genes)
  out
}

run_sc_fisher <- function(deg1, deg2, dge1, dge2, genes) {
  bg      <- intersect(dge1$gene, dge2$gene)
  
  shared  <- intersect(deg1$gene, deg2$gene)
  
  tbl     <- table(
    DE_shared = bg %in% shared,
    SC        = bg %in% genes
  )
  
  fisher.test(tbl)
  
}


run_fisher <- function(dge1, dge2, deg1, deg2) {
  bg  <- intersect(dge1$gene, dge2$gene)
  tbl <- table(
    DE_1 = bg %in% deg1$gene,
    DE_2 = bg %in% deg2$gene
  )
  fisher.test(tbl, alternative = "two.sided")
}


run_fishers_conconcate = function(asd_dge, sex_dge, looser = FALSE, remove_xy = FALSE, xy_genes = NULL, asd_thresh = list(logFC = 0.2, p = 0.05)) {
  
  sex_dge$Sex_Broad_Celltype_Br = paste(sex_dge$Sex_Broad_Celltype, sex_dge$brainregion_shorten, sep=".")
  
  ASD_DEGs <- filter_DEGs(asd_dge, "logFC_ASD_vs_CTL", "FDR_pVal_ASD_vs_CTL",
                          logFC_cut = asd_thresh$logFC, p_cut = asd_thresh$p,
                          remove_xy = remove_xy, xy_genes = xy_genes)
  
  SEX_DEGs <- filter_DEGs(sex_dge, "logFC", "adj.P.Val",
                          logFC_cut = ifelse(looser, 0, 0.2),
                          p_cut = ifelse(looser, 0.1, 0.05),
                          remove_xy = remove_xy, xy_genes = xy_genes)
  
  ASD_split <- split(asd_dge, asd_dge$ASD_Broad_Celltype)
  SEX_split <- split(sex_dge, sex_dge$Sex_Broad_Celltype_Br)
  ASD_DEG_split <- split(ASD_DEGs, ASD_DEGs$ASD_Broad_Celltype)
  SEX_DEG_split <- split(SEX_DEGs, SEX_DEGs$Sex_Broad_Celltype_Br)
  
  results <- expand.grid(ASD = names(ASD_DEG_split), Sex = names(SEX_DEG_split)) %>%
    pmap_dfr(function(ASD, Sex) {
      ASD_genes <- unique(ASD_DEG_split[[ASD]]$gene)
      SEX_genes <- unique(SEX_DEG_split[[Sex]]$gene)
      background <- intersect(unique(ASD_split[[ASD]]$gene), unique(SEX_split[[Sex]]$gene))
      
      if (length(ASD_genes) == 0 || length(SEX_genes) == 0 || length(background) == 0) {
        return(tibble(
          ASD_celltype = ASD, Sex_celltype = Sex, overlap_genes = 0,
          total_ASD_DEGs = length(ASD_genes), total_Sex_DEGs = length(SEX_genes),
          total_background = length(background), odds_ratio = NA, p_value = NA
        ))
      }
      
      df <- tibble(Gene = background,
                   DE_ASD = background %in% ASD_genes,
                   DE_Sex = background %in% SEX_genes)
      
      tbl <- table(df$DE_Sex, df$DE_ASD)
      
      if (nrow(tbl) < 2 || ncol(tbl) < 2) {
        return(tibble(
          ASD_celltype = ASD, 
          Sex_celltype = Sex,
          overlap_genes = length(intersect(ASD_genes, SEX_genes)),
          total_ASD_DEGs = length(ASD_genes), 
          total_Sex_DEGs = length(SEX_genes),
          total_background = length(background), 
          odds_ratio = NA, 
          p_value = NA
        ))
      }
      
      test <- fisher.test(tbl)
      
      tibble(
        ASD_celltype = ASD, 
        Sex_celltype = Sex,
        overlap_genes = length(intersect(ASD_genes, SEX_genes)),
        total_ASD_DEGs = length(ASD_genes), 
        total_Sex_DEGs = length(SEX_genes),
        total_background = length(background),
        odds_ratio = test$estimate, 
        p_value = test$p.value
      )
    })
  
  results %>% mutate(FDR = p.adjust(p_value, "fdr"), Sig = FDR < 0.05)
}

run_fisher_enrichment <- function(query_genes, background_genes, pathways) {
  
  results <- lapply(names(pathways), function(pathway_name) {
    
    pathway_genes <- pathways[[pathway_name]]
    
    in_query <- background_genes %in% query_genes
    in_pathway <- background_genes %in% pathway_genes
    
    contingency_table <- table(in_query, in_pathway)
    
    ft <- fisher.test(contingency_table, alternative = "greater")
    
    return(data.frame(
      Pathway = pathway_name,
      P_Value = ft$p.value,
      Odds_Ratio = ft$estimate
    ))
  })
  
  do.call(rbind, results)
}


#W permutations
fisher_list = list()
run_perm_fishers <- function(asd_dge, perm_list, ...) {
  
  for (i in seq_along(perm_list)) {
    message("Processing permutation ", i, " of ", length(perm_list), "...")
    
    df <- perm_list[[i]]
    
    df$brainregion_celltype <- df$Celltype_BrainRegion
    df$Sex_Brainregion <- sub(".*\\.", "", df$brainregion_celltype)
    df$Sex_Celltype    <- sub("\\..*", "", df$brainregion_celltype)
    
    df$brainregion_shorten <- dplyr::recode(
      df$Sex_Brainregion,
      "Angular_Gyrus"="AngGyr",
      "Caudal_Insula"= "CauIns",
      "Fusiform_Gyrus"="FusGyr",
      "Inferior_Lateral_TemporalOccipital"="InfLTO",
      "Retrosplenial_cortex"= "RetCor",
      "Intraparietal_sulcus" = "IntParSul")
    
    df <- df %>%
      dplyr::mutate(Sex_Broad_Celltype = case_when(
        Sex_Celltype %in% c("L23IT","L4IT","L5IT","L56NP","L6IT","L6b","L6CT",
                            "L5ET/L6IT","L6CT/L5IT","L23IT/L5IT") ~ "ExcNeu",
        Sex_Celltype %in% c("VIP","PVALB","SST","LAMP5","ADARB2","PAX6",
                            "Chandelier","Sncg") ~ "InhNeu",
        Sex_Celltype == "Astro" ~ "Astro",
        Sex_Celltype == "Micro" ~ "Micro",
        Sex_Celltype %in% c("Oligo","OPC") ~ "Oligo/OPC",
        Sex_Celltype == "BEC" ~ "BEC_Endo",
        TRUE ~ "Other"
      ))
    
    # Run fisher's test as above, but this produces Inf sometimes, since there are some perfect overlaps of low sex DEGs with autism DEGs 
    fisher_res <- run_fishers_conconcate(asd_dge, df, looser = TRUE)
    fisher_res$permutation_number <- i
    fisher_list[[i]] <- fisher_res
  }
  names(fisher_list) <- paste0("Permutation_", seq_along(perm_list))
  return(fisher_list)
}

add_columns_perm_fishers <- function(fisher_list) {
  lapply(fisher_list, function(df) {
    df$Sex_Brainregion <- sub(".*\\.", "", df$Sex_celltype)
    df$Sex_Celltype <- sub("\\..*", "", df$Sex_celltype)
    
    df$Sex_Brainregion_short <- dplyr::recode(
      df$Sex_Brainregion,
      "Angular_Gyrus"="AngGyr",
      "Caudal_Insula"= "CauIns",
      "Fusiform_Gyrus"="FusGyr",
      "Inferior_Lateral_TemporalOccipital"="InfLTO",
      "Retrosplenial_cortex"= "RetCor",
      "Intraparietal_sulcus" = "IntParSul"
    )
    df
  })
}

run_perm_fishers_NOY <- function(asd_dge, perm_list, gene_list_toexclude, ...) {
  
  for (i in seq_along(perm_list)) {
    message("Processing permutation ", i, " of ", length(perm_list), "...")
    
    df <- perm_list[[i]]
    
    df$brainregion_celltype <- df$Celltype_BrainRegion
    df$Sex_Brainregion <- sub(".*\\.", "", df$brainregion_celltype)
    df$Sex_Celltype    <- sub("\\..*", "", df$brainregion_celltype)
    
    df$brainregion_shorten <- dplyr::recode(
      df$Sex_Brainregion,
      "Angular_Gyrus"="AngGyr",
      "Caudal_Insula"= "CauIns",
      "Fusiform_Gyrus"="FusGyr",
      "Inferior_Lateral_TemporalOccipital"="InfLTO",
      "Retrosplenial_cortex"= "RetCor",
      "Intraparietal_sulcus" = "IntParSul")
    
    df <- df %>%
      dplyr::mutate(Sex_Broad_Celltype = case_when(
        Sex_Celltype %in% c("L23IT","L4IT","L5IT","L56NP","L6IT","L6b","L6CT",
                            "L5ET/L6IT","L6CT/L5IT","L23IT/L5IT") ~ "ExcNeu",
        Sex_Celltype %in% c("VIP","PVALB","SST","LAMP5","ADARB2","PAX6",
                            "Chandelier","Sncg") ~ "InhNeu",
        Sex_Celltype == "Astro" ~ "Astro",
        Sex_Celltype == "Micro" ~ "Micro",
        Sex_Celltype %in% c("Oligo","OPC") ~ "Oligo/OPC",
        Sex_Celltype == "BEC" ~ "BEC_Endo",
        TRUE ~ "Other"
      ))
    
    df_filt = df %>%
      dplyr::filter(!gene %in% gene_list_toexclude)
    
    asd_dge_filt = asd_dge %>%
      dplyr::filter(!gene %in% gene_list_toexclude)
    
    # Run fisher's test as above, but this produces Inf sometimes, since there are some perfect overlaps of low sex DEGs with autism DEGs 
    fisher_res <- run_fishers_conconcate(asd_dge_filt, df_filt, looser = TRUE)
    
    fisher_res$permutation_number <- i
    
    fisher_list[[i]] <- fisher_res
  }
  
  names(fisher_list) <- paste0("Permutation_", seq_along(perm_list))
  
  return(fisher_list)
  
}

add_columns_perm_fishers <- function(fisher_list) {
  lapply(fisher_list, function(df) {
    df$Sex_Brainregion <- sub(".*\\.", "", df$Sex_celltype)
    df$Sex_Celltype <- sub("\\..*", "", df$Sex_celltype)
    
    df$Sex_Brainregion_short <- dplyr::recode(
      df$Sex_Brainregion,
      "Angular_Gyrus"="AngGyr",
      "Caudal_Insula"= "CauIns",
      "Fusiform_Gyrus"="FusGyr",
      "Inferior_Lateral_TemporalOccipital"="InfLTO",
      "Retrosplenial_cortex"= "RetCor",
      "Intraparietal_sulcus" = "IntParSul"
    )
    df
  })
}

# ============================================================
# GO enrichment helpers
# ============================================================

run_enrichGO <- function(deg1, deg2, dge1, dge2, remove_genes = NULL) {
  shared <- intersect(deg1$gene, deg2$gene)
  
  bg     <- intersect(dge1$gene, dge2$gene)
  
  print("before removal")
  print(length(shared))
  
  if (!is.null(remove_genes)) {
    shared <- shared[!shared %in% remove_genes]
    bg     <- bg[!bg %in% remove_genes]
  }
  print("after removal")
  print(length(shared))
  enrch <- enrichGO(get_entrez(shared), ont = "BP",OrgDb = "org.Hs.eg.db",pvalueCutoff  = 0.05, pAdjustMethod = "BH",universe      = get_entrez(bg),minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2
  )
  list(
    Full_Enrich = enrch,
    Sig_Enrich  = enrch@result[enrch@result$p.adjust < 0.05, ]
  )
}

plot_enrichGO = function(enrichGO_res, title, fill_color) { 
  
  sig_enrichGO_res  = enrichGO_res@result[enrichGO_res@result$p.adjust < 0.05, ]
  
  sig_enrichGO_res_top = sig_enrichGO_res[order(sig_enrichGO_res$zScore, decreasing=T),] [1:10,]
  
  sig_enrichGO_res_top = as.data.frame(sig_enrichGO_res_top)
  
  sig_enrichGO_res_top = sig_enrichGO_res_top[order(sig_enrichGO_res_top$zScore, decreasing=F),]
  
  sig_enrichGO_res_top$Description = factor(sig_enrichGO_res_top$Description,levels = unique(sig_enrichGO_res_top$Description))
  
  sig_enrichGO_res_top$`Gene List` = str_wrap(sig_enrichGO_res_top$Description, width =30)
  
  sig_enrichGO_res_top$`Gene List` = factor(sig_enrichGO_res_top$`Gene List`,levels = rev(unique(sig_enrichGO_res_top$`Gene List`) ))
  
  sig_enrichGO_res_top$significant = "Yes"
  
  sig_enrichGO_res_top$significant[sig_enrichGO_res_top$p.adjust > 0.05] = "No"
  
  p = ggplot(sig_enrichGO_res_top, aes(x = `Gene List` , y =zScore)) +
    geom_segment( aes(x=`Gene List`, xend=`Gene List`, y=0, yend=zScore, color = significant)) +
    geom_point(position = position_dodge(width = 0.6), aes(size = -log10(p.adjust), color = significant), shape =21, stroke = 1.5,fill = fill_color) +
    scale_color_manual(values = c("Yes" = "black", "No" = "white")) +
    labs(title = title,x = "Gene List",y = "zScore") +
    theme_cowplot() +
    scale_size_continuous(name   = "-log10(FDR q)",range  = c(5, 30),limits=c(0,100)) +
    theme_cowplot() +
    theme(
      axis.text.x        = element_text(angle = 40, hjust = 1, size = 16),
      axis.text.y        = element_text(size = 18),
      axis.title         = element_text(size = 20),
      plot.margin        = margin(t = 10, r = 10, b = 20, l = 40, unit = "pt"),
      panel.grid.major.x = element_line(color = "grey92"),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(hjust = 0.5, size = 20, face = "bold"),
      legend.position    = "right"
    )+
    ylim(0,15)
  
  return(p)
  
}


get_entrez = function(genes) {
  mapIds(org.Hs.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
}


run_enrichGO <- function(deg1, deg2, dge1, dge2, remove_genes = NULL) {
  
  shared <- intersect(deg1$gene, deg2$gene)
  
  bg     <- intersect(dge1$gene, dge2$gene)
  
  print("before removal")
  print(length(shared))
  
  if (!is.null(remove_genes)) {
    shared <- shared[!shared %in% remove_genes]
    bg     <- bg[!bg %in% remove_genes]
  }
  
  print("after removal")
  print(length(shared))
  enrch <- enrichGO(get_entrez(shared), ont = "BP",OrgDb = "org.Hs.eg.db",pvalueCutoff  = 0.05, pAdjustMethod = "BH",universe      = get_entrez(bg),minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2
  )
  
  list(
    Full_Enrich = enrch,
    Sig_Enrich  = enrch@result[enrch@result$p.adjust < 0.05, ]
  )
  
}

plot_enrichGO = function(enrichGO_res, title, fill_color) { 
  
  sig_enrichGO_res  = enrichGO_res@result[enrichGO_res@result$p.adjust < 0.05, ]
  
  sig_enrichGO_res_top = sig_enrichGO_res[order(sig_enrichGO_res$zScore, decreasing=T),] [1:10,]
  
  sig_enrichGO_res_top = as.data.frame(sig_enrichGO_res_top)
  
  sig_enrichGO_res_top = sig_enrichGO_res_top[order(sig_enrichGO_res_top$zScore, decreasing=F),]
  
  sig_enrichGO_res_top$Description = factor(sig_enrichGO_res_top$Description,levels = unique(sig_enrichGO_res_top$Description))
  
  sig_enrichGO_res_top$`Gene List` = str_wrap(sig_enrichGO_res_top$Description, width =30)
  
  sig_enrichGO_res_top$`Gene List` = factor(sig_enrichGO_res_top$`Gene List`,levels = rev(unique(sig_enrichGO_res_top$`Gene List`) ))
  
  sig_enrichGO_res_top$significant = "Yes"
  
  sig_enrichGO_res_top$significant[sig_enrichGO_res_top$p.adjust > 0.05] = "No"
  
  p = ggplot(sig_enrichGO_res_top, aes(x = `Gene List` , y =zScore)) +
    geom_point(position = position_dodge(width = 0.6), aes(size = -log10(p.adjust), color = significant), shape =21, stroke = 1.5,fill = fill_color) +
    scale_color_manual(values = c("Yes" = "black", "No" = "white")) +
    labs(title = title,x = "Gene List",y = "zScore") +
    theme_cowplot() +
    scale_size_continuous(name   = "-log10(FDR q)",range  = c(5, 30),limits=c(0,100)) +
    theme_cowplot() +
    theme(
      axis.text.x        = element_text(angle = 40, hjust = 1, size = 16),
      axis.text.y        = element_text(size = 18),
      axis.title         = element_text(size = 20),
      plot.margin        = margin(t = 10, r = 10, b = 20, l = 40, unit = "pt"),
      panel.grid.major.x = element_line(color = "grey92"),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(hjust = 0.5, size = 20, face = "bold"),
      legend.position    = "right"
    )+
    ylim(0,15)
  
  return(p)
  
}


#Risk Genes
SFARI_genes <- read_xlsx("/vf/users/HGBomics/GabiDugan/OxCam/Human.Transcriptomics/SFARI-Gene_genes_07-08-2025release_09-26-2025export.xlsx",skip =1)

SFARI_genes= list("Autism - SFARI Score 1"=SFARI_genes$`gene-symbol`[SFARI_genes$`gene-score`=="1"],"Autism - SFARI Score 2"=SFARI_genes$`gene-symbol`[SFARI_genes$`gene-score`=="2"])

Fu_genes_list = list("Autism - Fu et al." = c("MED13", "NIPBL", "SMC3", "CNOT1", "BRD4", "FBN1", "XPO1", "PRR12", "SYNGAP1", "HDAC4", "SMARCA4", "AGO1", "ARHGAP35", "TRIP12", "KMT2E", "UBTF", "KDM3B", "TRRAP", "GRIN2B", "USP7", "ATP1A3",  "ASH1L", "MED13L", "KMT2A", "CREBBP", "KAT6A", "KMT2B", "SIN3A", "CHD2", "ARID1A", "SCN1A", "SPEN", "EHMT1", "TLE3", "UBR5", "PTCH1", "CHD7", "CNOT3", "AHDC1", "SUV420H1", "PHF21A", "HIVEP2", "DYNC1H1",  "CHD8", "WAC", "ZNF462", "WDR26", "STXBP1", "MN1", "MYT1L", "SPTBN1","SATB2", "EFTUD2", "TCF20", "DLL1", "NSD1", "SRCAP", "HNRNPK", "ARID2", "GATAD2B", "FBXO11", "NOTCH1", "TAB2", "EP300", "HDAC2", "SPTAN1", "CACNA1C", "ARID1B", "KMT2D", "ZNF148", "RFX3", "CLTC", "ANK2", "TBL1XR1", "ANKRD11", "PAFAH1B1", "ZEB2", "SETBP1", "NTRK2", "RORB", "PTPRD", "PHIP", "HNRNPU", "TLK2", "DHX30", "RAI1", "PIK3CA", "KIF5C", "SMARCC2", "SRRM2", "TSC1", "SON", "WHSC1", "STAG1","POGZ", "ZMYND11", "KIF11", "KMT2C", "ADNP", "SHANK3", "CACNA1E", "KCNB1", "NLGN2", "HECTD4", "SCN2A", "VCP", "COL4A1", "KAT6B", "CTNNB1", "SCN8A", "ZMYND8", "RFX7", "NR4A2", "U2AF2", "ITPR1",  "COL2A1", "CACNA1A", "PHF12", "CTR9", "PTPN11", "SETD1A", "KIF1A", "MBD5", "ZFHX3", "SHOC2", "TRIO", "KDM6B", "SOS1", "BAI1", "ZNF292", "CSNK1E", "SPOP", "ZFHX4", "GNB1", "ASXL2", "CHD3", "PPP1CB", "CTCF", "ATP2B2", "NFIX", "SLC6A1", "EBF3", "MAPK8IP3", "SETD1B", "CASZ1", "MARK2", "CHD4", "KCNQ2", "TANC2", "PPP3CA", "DPF2", "KIF2A", "UPF1", "ZNF865", "CASKIN1", "PAX5", "ARIH1", "GABBR2", "SOX9", "NACC1", "ACVR2A", "CAMTA1", "TFAP4", "TNRC6B", "CELF2","PPP2R5D", "FOXP1", "ABL1", "BRPF1", "PRPF8", "ZBTB18", "SET", "CIC", "MEIS2", "MTOR", "PACS1", "CTNND1", "EEF1A2", "SOX5", "GRIN2A", "DSCAM", "SMARCE1", "SKI", "TBR1", "GLI3", "GLYR1", "YY1", "ASXL3", "NFIA", "SMARCA2", "MORC2", "PDE4D", "UBE3A", "SETD2", "SYNCRIP", "SOX10", "BRAF", "CACNA1G", "MAP3K7", "RAB14", "DYRK1A", "FGFR1", "WDR48", "TCF4", "FOXP2", "SRP54", "SMARCB1", "AFF3", "SMAD4", "SF3B2", "SPAST", "NFIB", "SETD5", "HNRNPD", "CUL3", "PCBP2", "PIK3R1", "SNAP25", "FBXW7", "DLG4", "KANSL1","VEZF1", "CAMK2G", "SLC2A1", "SRSF1", "LONP1", "KLF7", "PUF60", "GRIA2", "TOP2B", "RYR2", "CAMK2A", "ZMIZ1", "DNM1", "AUTS2", "ADCY5", "NRXN1", "PBX1", "HECW2", "TRIM8", "EXT1", "ERF", "GNAO1", "CAMTA2", "PLXNA1", "PPP2CA", "MSI1", "GNB2", "ATP6V0A1", "RALA", "FGFR2", "CHAMP1", "TCF7L2", "TAOK1", "QRICH1", "KCNMA1", "BAP1",  "CTBP1", "PACS2", "KCND3", "NAA15", "BCL11B", "MFN2", "ATF2", "YWHAG", "ZBTB10", "NF1", "SATB1", "CSNK2A1", "PPP2R1A", "GLTSCR1",  "BICD2", "GRIN1", "CAPRIN1", "MAGEL2", "PIP5K1C", "BRD3", "ZBTB20",  "BCL11A", "CDK13", "ZBTB7A", "BRSK2", "FOXG1", "ZMYM2", "GSK3B", "HK1", "PUR A", "PSMC3", "GABRB3", "AP3B2", "KIDINS220", "CBX5", "PSMC5", "KCNQ3", "ATP6V1A", "DDX23", "COL4A3BP", "CAMK2B", "GABRA1", "PHLDB1", "TCF12", "GNAS", "KCNT2", "OTX2", "MAP2K1", "CTNNA2", "GNAI1", "KCNA2", "CSNK2B", "GABRB2", "KCNH1", "GATA3", "SRPR","GABRG2", "ARF1", "RORA", "IRF2BPL", "MYCN", "CDK8", "KIAA0100", "VAMP2", "POU3F3", "SOX11", "SYT1", "TRPM3", "TRAF7", "SHANK2", "SEMA6B", "RAC1", "KCNT1", "MSL2", "ABCC9", "GIGYF1", "PIK3R2","ATP1A2", "RHOBTB2", "TCF3", "MED23", "CDC42", "AP2S1", "DNM1L", "PTEN", "DHDDS", "KLHL20", "TNPO3", "NALCN", "PRKAR1B", "ODC1", "EIF5A", "SOX2", "KDM5B", "FGFR3", "MEF2C", "RAB11B", "ACTL6B", "RIT1", "DEAF1", "ASXL1", "POLR3B", "KBTBD7", "HRAS", "GFAP", "HIST1H2AC", "PPM1D", "H3F3A", "KRAS", "DNMT3A", "HIST1H1E", "HIST1H4J", "HIST1H2AE", "MIB1", "SMAD6", "LZTR1"))

Satterstrom_genes_list = list("Autism - Sattestrom et al." = c("CHD8", "SCN2A", "SYNGAP1", "ADNP", "FOXP1", "POGZ", "ARID1B", "SUV420H1", "DYRK1A", "SLC6A1", "GRIN2B", "PTEN", "SHANK3", "MED13L",  "GIGYF1", "CHD2", "ANKRD11", "ANK2", "ASH1L", "TLK2", "DNMT3A", "DEAF1", "CTNNB1", "KDM6B", "DSCAM", "SETD5", "KCNQ3", "SRPR", "KDM5B", "WAC", "SHANK2", "NRXN1", "TBL1XR1", "MYT1L", "BCL11A", "RORB", "RAI1", "DYNC1H1", "DPYSL2", "AP2S1", "KMT2C", "PAX5", "MKX", "GABRB3", "SIN3A", "MBD5", "MAP1A", "STXBP1", "CELF4", "PHF12", "TBR1", "PPP2R5D", "TM9SF4", "PHF21A", "PRR12", "SKI", "ASXL3", "SPAST", "SMARCC2", "TRIP12", "CREBBP", "TCF4", "CACNA1E",  "GNAI1", "TCF20", "FOXP2", "NSD1", "TCF7L2", "LDB1", "EIF3G", "PHF2", "KIAA0232", "VEZF1", "GFAP", "IRF2BPL", "ZMYND8", "SATB1", "RFX3", "SCN1A", "PPP5C", "TRIM23", "TRAF7", "ELAVL3", "GRIA2",  "LRRC4C", "CACNA2D3", "NUP155", "KMT2E", "NR3C2", "NACC1", "PTK7", "PPP1R9B", "GABRB2", "HDLBP", "TAOK1", "UBR1", "TEK", "KCNMA1",  "CORO1A", "HECTD4", "NCOA1", "DIP2A"))

Trubetskoy_Schizophrenia_riskgenes = list("Schizophrenia - Trubetskoy et al."= c("AC068490.2", "ACE", "ACTR1B", "AKT3", "ALMS1P", "ASPHD1", "ATG13", "ATP2A2", "BCL11B", "BCL2L12", "BNIP3L", "C12orf43", "CACNA1C", "CALN1", "CISD2", "CLCN3", "CNTN4", "CREB3L4", "CRHR1", "CSMD1", "CTD-2008L17.2", "CUL9", "CYP7B1", "DARS2", "DCC", "DCLK3", "DDHD2", "DLGAP2", "DNAJA3", "DPYD", "EMB", "EMX1", "ENOX1", "ENSG00000262319", "EPN2", "EYS", "FAM114A2", "FAM216A", "FAM83D", "FTSJ2", "FURIN", "GABBR2", "GATAD2A", "GPM6A", "GPR98", "GRAMD1B", "GRIN2A", "GRM1", "HYI", "IL1RAPL1", "IMMP2L", "INO80E", "IRF3", "KANSL1", "KIAA1549", "KLF6", "LINC00320", "LINC01068", "LINC01088", "LRRC4B", "LSM1", "MAD1L1", "MAN2A1", "MAPK3", "MAPT", "MLXIP", "MOB4", "MSI2", "MYO19", "NAB2", "NDFIP2", "NEBL", "NEGR1", "NLGN4X", "NRIP1", "NXPH1", "OPCML", "PAK6", "PCCB", "PCDHA2", "PCDHA8", "PCGF3", "PCNXL3", "PDE4B", "PDIA3", "PJA1", "PLCH2", "PSMA4", "PTPRD", "R3HDM2", "RERE", "RNASEH2C", "RP11-10L12.4", "RP11-165J3.6", "RP11-399D6.2", "RP11-490G2.2", "RP11-507B12.2", "RP11-53O19.3", "RP11-586K2.1", "RP11-73M18.2", "SETD6", "SF3B1", "SGCD", "SLC39A8", "SLC4A10", "SLC9B1", "SNAP91", "SP4", "THAP8", "THOC7", "TMTC1", "TRPC4", "TSNARE1", "TXNRD1", "WBSCR17", "WDR76", "WSCD2", "ZNF804A", "ZNF823", "ZNF835"))

Singh_Schizophrenia_riskgenes = list("Schizophrenia - Singh et al."= c("TRIO","SP4","CUL1","RB1CC1","HERC1","GRIN2A","CACNA1G","GRIA3","XPO7","SETD1A"))

combined_autism_genes <-list("Combined Autism"= unique(c(SFARI_genes$`Autism - SFARI Score 1`,SFARI_genes$`Autism - SFARI Score 2`, Satterstrom_genes_list$`Autism - Sattestrom et al.`, Fu_genes_list$`Autism - Fu et al.`)))

combined_schizophrenia_genes <-list("Combined Schizophrenia"= unique(c(Singh_Schizophrenia_riskgenes$`Schizophrenia - Singh et al.`, Trubetskoy_Schizophrenia_riskgenes$`Schizophrenia - Trubetskoy et al.`)))

riskgene_pathways = c(SFARI_genes, Fu_genes_list, Satterstrom_genes_list,Trubetskoy_Schizophrenia_riskgenes,combined_autism_genes,Singh_Schizophrenia_riskgenes,combined_schizophrenia_genes )

run_fisher_enrichment <- function(query_genes, background_genes, pathways) {
  
  results <- lapply(names(pathways), function(pathway_name) {
    
    pathway_genes <- pathways[[pathway_name]]
    
    in_query <- background_genes %in% query_genes
    in_pathway <- background_genes %in% pathway_genes
    
    contingency_table <- table(in_query, in_pathway)
    
    ft <- fisher.test(contingency_table, alternative = "greater")
    
    return(data.frame(
      Pathway = pathway_name,
      P_Value = ft$p.value,
      Odds_Ratio = ft$estimate
    ))
  })
  
  do.call(rbind, results)
}



#hennick et al. early development

####################################################
#Functions
####################################################
map_broad_celltype <- function(df, celltype_col) {
  df %>%
    dplyr::mutate(Broad_Celltype = case_when(
      !!sym(celltype_col) %in% c("EXT_6_L23", "EXT_3_L23", "EXT_11_L34", "EXT_10_L34", 
                                 "EXT_1_L4", "EXT_13_L45", "EXT_14_L5B", "EXT_4_L56", 
                                 "EXT_7_L2356", "EXT_5_L56", "EXT_2_L56", "EXT_12_L56", 
                                 "EXT_8_L6", "EXT_9_L6", "Ex-L2", "Ex-L23", "Ex-L3", 
                                 "Ex-L4_MYLK", "Ex-L45_MET", "Ex-L45_LRRK1", "Ex-L5b_HTR2C", 
                                 "Ex-L56", "Ex-L56_CC_NTNG2", "Ex-L6_CC_SEMA3A", "Ex-L6b_SEMA3D", 
                                 "Ex-L6b_SEMA3E", "L23IT", "L4IT", "L5IT", "L56NP", "L6IT", 
                                 "L6b", "L6CT", "L5ET/ L6IT", "L6CT / L5IT", "L23IT/ L5IT") ~ "ExcNeu",
      !!sym(celltype_col) %in% c("INT_5_SST", "INT_2_PV_SULF", "INT_6_PV", "INT_1_KIT", 
                                 "INT_3_RELN_ND", "INT_4_VIP_CR", "In-Rosehip_CHST9", 
                                 "In-Rosehip_TRPC3", "In-Reelin", "In-VIP", "In-PV_Chandelier", 
                                 "In-PV_Basket", "In-SST", "VIP", "PVALB", "SST", "LAMP5", 
                                 "ADARB2", "PAX6", "Chandelier", "Sncg") ~ "InhNeu",
      !!sym(celltype_col) %in% c("ASTRO_1", "ASTRO_2", "ASTRO_3", "Ast", "Astro") ~ "Astro",
      !!sym(celltype_col) %in% c("MG1", "MG2", "BBB_MG3", "Mic", "Micro") ~ "Micro",
      !!sym(celltype_col) %in% c("ODC_1", "ODC_2", "ODC_3", "OPCS_1", "OPCS_2", "OPCS_3", 
                                 "OPCS_4", "Oli", "OPC", "Oligo") ~ "Oligo/OPC",
      !!sym(celltype_col) %in% c("BBB_Endo_1", "BBB_Pericytes", "BBB_SM(Mural)", "Endo", 
                                 "Pericytes", "BEC") ~ "NeuVasc",
      TRUE ~ "Other"
    ))
}


summarize_dge <- function(df, logfc_col, dataset_name, group_cols = c("Broad_Celltype", "gene")) {
  df %>%
    dplyr::group_by(across(all_of(group_cols))) %>%
    dplyr::summarise(
      mean_logFC = mean(!!sym(logfc_col), na.rm = TRUE),
      Dataset = dataset_name,
      .groups = "drop"
    )
}


prep_gsea_list <- function(summarized_df, split_col) {
  split_data <- split(summarized_df, summarized_df[[split_col]])
  lapply(split_data, function(x) {
    ranks <- setNames(x$mean_logFC, x$gene)
    sort(na.omit(ranks), decreasing = TRUE)
  })
}


plot_enrichment <- function(df, title_text, alpha_col = "significant") {
  ggplot(df, aes(x = pathway, y = Broad_Celltype, size = abs(NES), fill = NES)) +
    geom_point(aes(alpha = !!sym(alpha_col)), shape = 21, stroke = 1.5) +
    facet_grid(rows = vars(BrainRegion), cols = vars(PrenatalSex), scales="free") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3)) +
    scale_fill_gradient2(low = "#2b83ba", mid = "#f7f7f7", high = "#b35806", midpoint = 0) +
    scale_size_continuous(limits = c(0, 3)) +
    theme_cowplot() +
    labs(x = "Development - Gene List", y = "Adulthood - Broad Celltype", title = title_text) +
    theme(axis.text.x = element_text(hjust = 1, angle = 45), plot.title = element_text(size=25, hjust =0.5), strip.text=element_text(size=18))
}

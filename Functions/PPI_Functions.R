# ===============================
# Read in data
# ===============================
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

#DGE
#Autism
asd_dge <- read_xlsx("~/Data_S1_CELLTYPE-SPECIFIC_RAW_DEG.xlsx", col_names = TRUE) %>%
  filter(!grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", gene)) %>%
  map_broad_celltype("CELL TYPE")

asd_mean <- summarize_dge(asd_dge, "logFC_ASD_vs_CTL", "Autism") %>%
  mutate(BrainRegion = "PFC")

#Sex
sex_dge <-readRDS("~/sex_dge") %>%
  map_broad_celltype("celltype_only")

sex_mean <- summarize_dge(sex_dge, "logFC", "Sex", 
                          group_cols = c("Broad_Celltype", "brainregion_shorten", "gene")) %>%
  mutate(Broad_Celltype_BR = paste(Broad_Celltype, brainregion_shorten, sep = "."))


#Sex and autism matrices that CCA was run on 
matrix1 =readRDS("~/matrix1")
matrix2= readRDS("~/matrix2")

#rCCA with and without XY genes (including PAR)
rcca_sex_schizophrenia =readRDS("~/rcca_sex_schizophrenia")
rcca_sex_schizophrenia_auto = readRDS("~/rcca_sex_schizophrenia_auto")
rcc.model_filtered_ASD =readRDS("~/rcca_sex_autism")
rcc.model_filtered_ASD_auto = readRDS("~/rcca_sex_autism_auto")
CCA1_weights <-read_xlsx("~/CCA1_weights.xlsx")

CCA1_weights = CCA1_weights %>%
  dplyr::mutate(mean = (CCA1_weights$CCA1_Weight_Sex+CCA1_weights$CCA1_Weight_ASD)/2)

CCA1_weights_auto = read_xlsx("~/CCA1_weights_auto.xlsx") 

#SFARI genes 
SFARI_genes <- read_xlsx("~/SFARI-Gene_genes_07-08-2025release_09-26-2025export.xlsx",skip =1)
SFARI_genes= list("Autism - SFARI Score 1"=SFARI_genes$`gene-symbol`[SFARI_genes$`gene-score`=="1"],"Autism - SFARI Score 2"=SFARI_genes$`gene-symbol`[SFARI_genes$`gene-score`=="2"])
Fu_genes_list = list("Autism - Fu et al." = c("MED13", "NIPBL", "SMC3", "CNOT1", "BRD4", "FBN1", "XPO1", "PRR12", "SYNGAP1", "HDAC4", "SMARCA4", "AGO1", "ARHGAP35", "TRIP12", "KMT2E", "UBTF", "KDM3B", "TRRAP", "GRIN2B", "USP7", "ATP1A3",  "ASH1L", "MED13L", "KMT2A", "CREBBP", "KAT6A", "KMT2B", "SIN3A", "CHD2", "ARID1A", "SCN1A", "SPEN", "EHMT1", "TLE3", "UBR5", "PTCH1", "CHD7", "CNOT3", "AHDC1", "SUV420H1", "PHF21A", "HIVEP2", "DYNC1H1",  "CHD8", "WAC", "ZNF462", "WDR26", "STXBP1", "MN1", "MYT1L", "SPTBN1","SATB2", "EFTUD2", "TCF20", "DLL1", "NSD1", "SRCAP", "HNRNPK", "ARID2", "GATAD2B", "FBXO11", "NOTCH1", "TAB2", "EP300", "HDAC2", "SPTAN1", "CACNA1C", "ARID1B", "KMT2D", "ZNF148", "RFX3", "CLTC", "ANK2", "TBL1XR1", "ANKRD11", "PAFAH1B1", "ZEB2", "SETBP1", "NTRK2", "RORB", "PTPRD", "PHIP", "HNRNPU", "TLK2", "DHX30", "RAI1", "PIK3CA", "KIF5C", "SMARCC2", "SRRM2", "TSC1", "SON", "WHSC1", "STAG1","POGZ", "ZMYND11", "KIF11", "KMT2C", "ADNP", "SHANK3", "CACNA1E", "KCNB1", "NLGN2", "HECTD4", "SCN2A", "VCP", "COL4A1", "KAT6B", "CTNNB1", "SCN8A", "ZMYND8", "RFX7", "NR4A2", "U2AF2", "ITPR1",  "COL2A1", "CACNA1A", "PHF12", "CTR9", "PTPN11", "SETD1A", "KIF1A", "MBD5", "ZFHX3", "SHOC2", "TRIO", "KDM6B", "SOS1", "BAI1", "ZNF292", "CSNK1E", "SPOP", "ZFHX4", "GNB1", "ASXL2", "CHD3", "PPP1CB", "CTCF", "ATP2B2", "NFIX", "SLC6A1", "EBF3", "MAPK8IP3", "SETD1B", "CASZ1", "MARK2", "CHD4", "KCNQ2", "TANC2", "PPP3CA", "DPF2", "KIF2A", "UPF1", "ZNF865", "CASKIN1", "PAX5", "ARIH1", "GABBR2", "SOX9", "NACC1", "ACVR2A", "CAMTA1", "TFAP4", "TNRC6B", "CELF2","PPP2R5D", "FOXP1", "ABL1", "BRPF1", "PRPF8", "ZBTB18", "SET", "CIC", "MEIS2", "MTOR", "PACS1", "CTNND1", "EEF1A2", "SOX5", "GRIN2A", "DSCAM", "SMARCE1", "SKI", "TBR1", "GLI3", "GLYR1", "YY1", "ASXL3", "NFIA", "SMARCA2", "MORC2", "PDE4D", "UBE3A", "SETD2", "SYNCRIP", "SOX10", "BRAF", "CACNA1G", "MAP3K7", "RAB14", "DYRK1A", "FGFR1", "WDR48", "TCF4", "FOXP2", "SRP54", "SMARCB1", "AFF3", "SMAD4", "SF3B2", "SPAST", "NFIB", "SETD5", "HNRNPD", "CUL3", "PCBP2", "PIK3R1", "SNAP25", "FBXW7", "DLG4", "KANSL1","VEZF1", "CAMK2G", "SLC2A1", "SRSF1", "LONP1", "KLF7", "PUF60", "GRIA2", "TOP2B", "RYR2", "CAMK2A", "ZMIZ1", "DNM1", "AUTS2", "ADCY5", "NRXN1", "PBX1", "HECW2", "TRIM8", "EXT1", "ERF", "GNAO1", "CAMTA2", "PLXNA1", "PPP2CA", "MSI1", "GNB2", "ATP6V0A1", "RALA", "FGFR2", "CHAMP1", "TCF7L2", "TAOK1", "QRICH1", "KCNMA1", "BAP1",  "CTBP1", "PACS2", "KCND3", "NAA15", "BCL11B", "MFN2", "ATF2", "YWHAG", "ZBTB10", "NF1", "SATB1", "CSNK2A1", "PPP2R1A", "GLTSCR1",  "BICD2", "GRIN1", "CAPRIN1", "MAGEL2", "PIP5K1C", "BRD3", "ZBTB20",  "BCL11A", "CDK13", "ZBTB7A", "BRSK2", "FOXG1", "ZMYM2", "GSK3B", "HK1", "PUR A", "PSMC3", "GABRB3", "AP3B2", "KIDINS220", "CBX5", "PSMC5", "KCNQ3", "ATP6V1A", "DDX23", "COL4A3BP", "CAMK2B", "GABRA1", "PHLDB1", "TCF12", "GNAS", "KCNT2", "OTX2", "MAP2K1", "CTNNA2", "GNAI1", "KCNA2", "CSNK2B", "GABRB2", "KCNH1", "GATA3", "SRPR","GABRG2", "ARF1", "RORA", "IRF2BPL", "MYCN", "CDK8", "KIAA0100", "VAMP2", "POU3F3", "SOX11", "SYT1", "TRPM3", "TRAF7", "SHANK2", "SEMA6B", "RAC1", "KCNT1", "MSL2", "ABCC9", "GIGYF1", "PIK3R2","ATP1A2", "RHOBTB2", "TCF3", "MED23", "CDC42", "AP2S1", "DNM1L", "PTEN", "DHDDS", "KLHL20", "TNPO3", "NALCN", "PRKAR1B", "ODC1", "EIF5A", "SOX2", "KDM5B", "FGFR3", "MEF2C", "RAB11B", "ACTL6B", "RIT1", "DEAF1", "ASXL1", "POLR3B", "KBTBD7", "HRAS", "GFAP", "HIST1H2AC", "PPM1D", "H3F3A", "KRAS", "DNMT3A", "HIST1H1E", "HIST1H4J", "HIST1H2AE", "MIB1", "SMAD6", "LZTR1"))
Satterstrom_genes_list = list("Autism - Sattestrom et al." = c("CHD8", "SCN2A", "SYNGAP1", "ADNP", "FOXP1", "POGZ", "ARID1B", "SUV420H1", "DYRK1A", "SLC6A1", "GRIN2B", "PTEN", "SHANK3", "MED13L",  "GIGYF1", "CHD2", "ANKRD11", "ANK2", "ASH1L", "TLK2", "DNMT3A", "DEAF1", "CTNNB1", "KDM6B", "DSCAM", "SETD5", "KCNQ3", "SRPR", "KDM5B", "WAC", "SHANK2", "NRXN1", "TBL1XR1", "MYT1L", "BCL11A", "RORB", "RAI1", "DYNC1H1", "DPYSL2", "AP2S1", "KMT2C", "PAX5", "MKX", "GABRB3", "SIN3A", "MBD5", "MAP1A", "STXBP1", "CELF4", "PHF12", "TBR1", "PPP2R5D", "TM9SF4", "PHF21A", "PRR12", "SKI", "ASXL3", "SPAST", "SMARCC2", "TRIP12", "CREBBP", "TCF4", "CACNA1E",  "GNAI1", "TCF20", "FOXP2", "NSD1", "TCF7L2", "LDB1", "EIF3G", "PHF2", "KIAA0232", "VEZF1", "GFAP", "IRF2BPL", "ZMYND8", "SATB1", "RFX3", "SCN1A", "PPP5C", "TRIM23", "TRAF7", "ELAVL3", "GRIA2",  "LRRC4C", "CACNA2D3", "NUP155", "KMT2E", "NR3C2", "NACC1", "PTK7", "PPP1R9B", "GABRB2", "HDLBP", "TAOK1", "UBR1", "TEK", "KCNMA1",  "CORO1A", "HECTD4", "NCOA1", "DIP2A"))
combined_autism_genes <-list("Combined Autism"= unique(c(SFARI_genes$`Autism - SFARI Score 1`,SFARI_genes$`Autism - SFARI Score 2`, Satterstrom_genes_list$`Autism - Sattestrom et al.`, Fu_genes_list$`Autism - Fu et al.`)))
Trubetskoy_Schizophrenia_riskgenes = list("Schizophrenia - Trubetskoy et al."= c("AC068490.2", "ACE", "ACTR1B", "AKT3", "ALMS1P", "ASPHD1", "ATG13", "ATP2A2", "BCL11B", "BCL2L12", "BNIP3L", "C12orf43", "CACNA1C", "CALN1", "CISD2", "CLCN3", "CNTN4", "CREB3L4", "CRHR1", "CSMD1", "CTD-2008L17.2", "CUL9", "CYP7B1", "DARS2", "DCC", "DCLK3", "DDHD2", "DLGAP2", "DNAJA3", "DPYD", "EMB", "EMX1", "ENOX1", "ENSG00000262319", "EPN2", "EYS", "FAM114A2", "FAM216A", "FAM83D", "FTSJ2", "FURIN", "GABBR2", "GATAD2A", "GPM6A", "GPR98", "GRAMD1B", "GRIN2A", "GRM1", "HYI", "IL1RAPL1", "IMMP2L", "INO80E", "IRF3", "KANSL1", "KIAA1549", "KLF6", "LINC00320", "LINC01068", "LINC01088", "LRRC4B", "LSM1", "MAD1L1", "MAN2A1", "MAPK3", "MAPT", "MLXIP", "MOB4", "MSI2", "MYO19", "NAB2", "NDFIP2", "NEBL", "NEGR1", "NLGN4X", "NRIP1", "NXPH1", "OPCML", "PAK6", "PCCB", "PCDHA2", "PCDHA8", "PCGF3", "PCNXL3", "PDE4B", "PDIA3", "PJA1", "PLCH2", "PSMA4", "PTPRD", "R3HDM2", "RERE", "RNASEH2C", "RP11-10L12.4", "RP11-165J3.6", "RP11-399D6.2", "RP11-490G2.2", "RP11-507B12.2", "RP11-53O19.3", "RP11-586K2.1", "RP11-73M18.2", "SETD6", "SF3B1", "SGCD", "SLC39A8", "SLC4A10", "SLC9B1", "SNAP91", "SP4", "THAP8", "THOC7", "TMTC1", "TRPC4", "TSNARE1", "TXNRD1", "WBSCR17", "WDR76", "WSCD2", "ZNF804A", "ZNF823", "ZNF835"))

Singh_Schizophrenia_riskgenes = list("Schizophrenia - Singh et al."= c("TRIO","SP4","CUL1","RB1CC1","HERC1","GRIN2A","CACNA1G","GRIA3","XPO7","SETD1A"))

combined_schizophrenia_genes <-list("Combined Schizophrenia"= unique(c(Singh_Schizophrenia_riskgenes$`Schizophrenia - Singh et al.`, Trubetskoy_Schizophrenia_riskgenes$`Schizophrenia - Trubetskoy et al.`)))


#Sex-chromosome genes 
xy_genes = readRDS("~/Redraft_February2026/xy_genes")
x_genes = readRDS("~/Redraft_February2026/x_genes")
y_genes = readRDS("~/Redraft_February2026/y_genes")
par_1_genes <- c("AKAP17A", "ASMTL", "ASMTL-AS1", "CD99", "CD99P1", "CRLF2", "CSF2RA", "DHRSX", "DHRSX-IT1", "FABP5P13", "GTPBP6", "IL3RA", "LINC00102", "LINC00106", "LINC00685", "MIR3690", "MIR6089", "PLCXD1", "PPP2R3B", "P2RY8", "SHOX", "SLC25A6", "XG", "ZBED1")
par_2_genes <- c("AMD1P2", "DDX11L16", "DPH3P2", "IL9R", "SPRY3", "ELOCP24", "TRPC6P1", "VAMP7", "WASH6P", "WASIR1")
par_genes = c(par_1_genes,par_2_genes)
xy_genes_nopar = xy_genes[!xy_genes %in% c(par_1_genes,par_2_genes)]
x_genes_nopar = x_genes[!x_genes %in% c(par_1_genes,par_2_genes)]
y_genes_nopar = y_genes[!y_genes %in% c(par_1_genes,par_2_genes)]

# Psychencode developmental transcription
df_psychencode_metadata <- read_excel("~/PsychEncode_DevTraj_Metadata.xlsx")
mrna_exp <- read.table("~/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.count.txt", head = T, row.names = 1)

mrna_exp_log <- log10(mrna_exp+1) 
mrna_exp_centered <- scale(mrna_exp_log, center = TRUE, scale = FALSE)+1
split_names <- strsplit(colnames(mrna_exp_centered), "\\.")

mrna_psychencode_metadata <- do.call(rbind, lapply(split_names, function(x) {
  data.frame(sample_id = x[1], brainregion = x[2], stringsAsFactors = FALSE)
}))
mrna_psychencode_metadata$column_names <- paste(mrna_psychencode_metadata$sample_id,mrna_psychencode_metadata$brainregion,sep=".") 
mrna_psychencode_metadata$window <- df_psychencode_metadata$Window[match(mrna_psychencode_metadata$sample_id,df_psychencode_metadata$Braincode)]

dput(unique(mrna_psychencode_metadata$brainregion))
allregions <- c("A1C", "AMY", "CBC", "DFC", "HIP", "IPC", "ITC", "M1C", "MFC", 
                "OFC", "S1C", "STR", "V1C", "VFC", "STC", "MD", "CGE", "DTH", 
                "LGE", "M1CS1C", "MGE", "OC", "PC", "URL", "TC", "MSC")
cortical_regions <- c("A1C","CBC","DFC","IPC","ITC","M1C","MFC","OFC","S1C","V1C","VFC","M1CS1C","OC","PC","TC","MSC")

df_sample_window <- mrna_psychencode_metadata %>% 
  dplyr::filter(brainregion %in% cortical_regions)

pattern <- paste(cortical_regions, collapse = "|")
mrna_cortex <- mrna_exp_centered[, grepl(pattern, colnames(mrna_exp_centered))]
rownames(mrna_cortex) =  sub(".*\\|", "", rownames(mrna_cortex))


# ===============================
# Functions, dataframe clean up & network generation
# ===============================

annotate_broad_celltype_sex <- function(celltype) {
  dplyr::case_when(
    celltype %in% c("L23IT","L4IT","L5IT","L56NP","L6IT","L6b","L6CT",
                    "L5ET/ L6IT","L6CT / L5IT","L23IT/ L5IT") ~ "ExcNeu",
    celltype %in% c("VIP","PVALB","SST","LAMP5","ADARB2","PAX6",
                    "Chandelier","Sncg") ~ "InhNeu",
    celltype %in% c("Astro") ~ "Astro",
    celltype %in% c("Micro") ~ "Micro",
    celltype %in% c("Oligo","OPC") ~ "Oligo/OPC",
    celltype %in% c("BEC") ~ "NeuVasc",
    TRUE ~ NA_character_
  )
}

annotate_broad_celltype_autism <- function(celltype) {
  dplyr::case_when(
    celltype %in% c("EXT_6_L23","EXT_3_L23","EXT_11_L34","EXT_10_L34",
                    "EXT_1_L4","EXT_13_L45","EXT_14_L5B","EXT_4_L56",
                    "EXT_7_L2356","EXT_5_L56","EXT_2_L56","EXT_12_L56",
                    "EXT_8_L6","EXT_9_L6") ~ "ExcNeu",
    celltype %in% c("INT_5_SST","INT_2_PV_SULF","INT_6_PV",
                    "INT_1_KIT","INT_3_RELN_ND","INT_4_VIP_CR") ~ "InhNeu",
    celltype %in% c("ASTRO_1","ASTRO_2","ASTRO_3") ~ "Astro",
    celltype %in% c("MG1","MG2","BBB_MG3") ~ "Micro",
    celltype %in% c("ODC_1","ODC_2","ODC_3",
                    "OPCS_1","OPCS_2","OPCS_3","OPCS_4") ~ "Oligo/OPC",
    celltype %in% c("BBB_Endo_1","BBB_Pericytes","BBB_SM(Mural)") ~ "NeuVasc",
    TRUE ~ NA_character_
  )
}

annotate_brainregion_short <- function(brainregion) {
  dplyr::recode(
    brainregion,
    "Angular_Gyrus" = "AngGyr",
    "Caudal_Insula" = "CauIns",
    "Fusiform_Gyrus" = "FusGyr",
    "Inferior_Lateral_TemporalOccipital" = "InfLTO",
    "Intraparietal_sulcus" = "IntParSul",
    "Retrosplenial_cortex" = "RetCor",
    "Frontal Cortex" = "Frontal Cortex"
  )
}

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


#Same function as for dev. traj
get_all_top_gene_subsets <- function(gene_scores, score_col = "mean_score", top_n) {
  extract_top_genes_expr <- function(gene_scores, score_col, direction = "up", top_n) {
    if (direction == "up") {
      top_genes <- gene_scores[order(-gene_scores[[score_col]]), ][1:top_n, ]
    } else if (direction == "down") {
      top_genes <- gene_scores[order(gene_scores[[score_col]]), ][1:top_n, ]
    } else {
      stop("Invalid direction. Use 'up' or 'down'.")
    }
    return(top_genes)
  }
  list(
    CCA1_up = extract_top_genes_expr(gene_scores, score_col, "up", top_n),
    CCA1_down = extract_top_genes_expr(gene_scores, score_col, "down", top_n)
  )
}


STRINGdb_networks = function(dataframe, gene_column, seed = 124){
  
  set.seed(seed)
  
  gene_data <- as.data.frame(dataframe)
  
  mapped <- string_db$map(gene_data, "gene", removeUnmappedRows = TRUE)
  
  interactions <- string_db$get_interactions(mapped$STRING_id)
  
  network <- interactions
  
  network$from = mapped$gene[match(network$from,mapped$STRING_id)]
  
  network$to = mapped$gene[match(network$to,mapped$STRING_id)]
  
  #Graph PPI networks
  graph_from_df_network <- graph_from_data_frame(network[, c("from", "to")], directed = FALSE)
  
  degree_df <- data.frame(Degree= degree(graph_from_df_network, mode = "all"))
  
  degree_df =rownames_to_column(degree_df, "Gene")
  
  return(list("Mapped"= mapped,
              "Network"= network,
              "Graph_df"= graph_from_df_network,
              "Degree" = degree_df))
}

degree_matched_permutation <- function(graph_network,
                                       obs_genes,
                                       n_permutations = 1000,
                                       n_bins = 20,
                                       stat = c("n_edges", "largest_comp", "mean_degree"),
                                       allow_repeats = FALSE,
                                       seed = 24,
                                       k_nearest = 50) {
  
  stopifnot(inherits(graph_network, "igraph"))
  if (!"Gene" %in% names(obs_genes)) stop("obs_genes must have a 'Gene' column")
  
  stat <- match.arg(stat)
  set.seed(seed)
  
  all_nodes <- V(graph_network)$name
  if (is.null(all_nodes)) stop("graph_network vertices must have a 'name' attribute")
  
  obs_gene_vec <- unique(intersect(obs_genes$Gene, all_nodes))
  
  if (length(obs_gene_vec) == 0) stop("No observed genes present in graph vertices")
  if (length(obs_gene_vec) < length(unique(obs_genes$Gene))) {
    warning(sprintf(
      "%d observed gene(s) not found in graph and were dropped",
      length(unique(obs_genes$Gene)) - length(obs_gene_vec)
    ))
  }
  
  # --- degree table & bins -----------------------------------------------
  deg_df <- data.frame(
    Gene   = all_nodes,
    Degree = degree(graph_network, mode = "all"),
    stringsAsFactors = FALSE
  )
  
  probs  <- seq(0, 1, length.out = n_bins + 1)
  breaks <- unique(quantile(deg_df$Degree, probs = probs, type = 1))
  
  if (length(breaks) < 2) {
    deg_df$bin <- "1"
  } else {
    deg_df$bin <- as.character(cut(deg_df$Degree, breaks = breaks, include.lowest = TRUE))
  }
  
  candidates_by_bin <- split(deg_df$Gene, deg_df$bin)
  
  obs_info <- deg_df[deg_df$Gene %in% obs_gene_vec, c("Gene", "Degree", "bin")]
 
   # preserve the order of obs_gene_vec
  obs_info <- obs_info[match(obs_gene_vec, obs_info$Gene), ]
  
  n_obs <- nrow(obs_info)
  if (length(all_nodes) < n_obs) stop("Graph has fewer nodes than observed genes")
  
  # genes that must never appear in a permuted set, regardless of branch
  real_genes <- obs_info$Gene
  
  # k_nearest fallback pool for a single gene, given a set of genes to exclude
  fallback_pool <- function(gene, gene_degree, exclude_set) {
    deg_diff <- abs(deg_df$Degree - gene_degree)
    ord <- deg_df$Gene[order(deg_diff, deg_df$Gene)]
    ord <- setdiff(ord, c(gene, exclude_set))
    if (length(ord) == 0) {
      stop(sprintf(
        "Could not find any eligible background gene for '%s' (all candidates excluded).",
        gene
      ))
    }
    ord[seq_len(min(k_nearest, length(ord)))]
  }
  
  
  draw_one_set <- function() {
    sampled <- character(n_obs)
    names(sampled) <- obs_info$Gene
    
    for (i in seq_len(n_obs)) {
      g <- obs_info$Gene[i]
      d <- obs_info$Degree[i]
      b <- obs_info$bin[i]
      
      #can't have duplicate genes within the same perm network
      already_chosen <- if (allow_repeats) character(0) else sampled[seq_len(i - 1)]
      exclude_set <- union(real_genes, already_chosen)
      
      pool <- setdiff(candidates_by_bin[[b]], exclude_set)
      
      if (length(pool) == 0) {
        pool <- fallback_pool(g, d, exclude_set)
      }
      
      sampled[i] <- sample(pool, 1)
    }
    
    sampled
  }
  
  # --- statistic on a gene set ---------------------------------------------
  compute_stat <- function(gene_set) {
    gene_set <- unique(gene_set)
    gene_set_in_sub <- intersect(gene_set, all_nodes)
    if (length(gene_set_in_sub) < 2) return(0)
    
    g_sub <- induced_subgraph(graph_network, vids = gene_set_in_sub)
    
    if (stat == "n_edges")      return(gsize(g_sub))
    if (stat == "largest_comp") return(max(components(g_sub)$csize))
    if (stat == "mean_degree") {
      d <- degree(g_sub)
      return(if (length(d) == 0) 0 else mean(d))
    }
  }
  
  obs_stat <- compute_stat(obs_gene_vec)
  
  perm_stats <- numeric(n_permutations)
  for (p in seq_len(n_permutations)) {
    perm_set <- draw_one_set()
    perm_stats[p] <- compute_stat(perm_set)
  }
  
  p_emp <- (sum(perm_stats >= obs_stat) + 1) / (n_permutations + 1)
  
  list(
    obs_stat       = obs_stat,
    perm_stats     = perm_stats,
    p_val          = p_emp,
    n_permutations = n_permutations,
    stat           = stat
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


run_fishers_conconcate = function(hub_genes, sex_dge, looser = FALSE, remove_xy = FALSE, xy_genes = NULL) {
  
  sex_dge$Sex_Broad_Celltype_Br = paste(sex_dge$Sex_Broad_Celltype, sex_dge$brainregion_shorten, sep=".")
  
  SEX_DEGs <- filter_DEGs(sex_dge, "logFC", "adj.P.Val",
                          logFC_cut = ifelse(looser, 0, 0.2),
                          p_cut = ifelse(looser, 0.1, 0.05),
                          remove_xy = remove_xy, xy_genes = xy_genes)
  
  SEX_split <- split(sex_dge, sex_dge$Sex_Broad_Celltype_Br)
  
  SEX_DEG_split <- split(SEX_DEGs, SEX_DEGs$Sex_Broad_Celltype_Br)
  
  results <- expand.grid(Sex = names(SEX_DEG_split), stringsAsFactors = FALSE) %>%
    
    pmap_dfr(function(Sex) {
      
      SEX_genes <- unique(SEX_DEG_split[[Sex]]$gene)
      
      background <- unique(SEX_split[[Sex]]$gene)
      
      if (length(SEX_genes) == 0 || length(background) == 0) {
        return(tibble(
          Sex_celltype = Sex, 
          overlap_genes = 0,
          total_hub = length(hub_genes),
          total_Sex_DEGs = length(SEX_genes),
          total_background = length(background), odds_ratio = NA, p_value = NA
        ))
      }
      
      df <- tibble(Gene = background,
                   hub = background %in% hub_genes,
                   DE_Sex = background %in% SEX_genes)
      
      tbl <- table(df$DE_Sex, df$hub)
      
      if (nrow(tbl) < 2 || ncol(tbl) < 2) {
        return(tibble(
          Sex_celltype = Sex,
          overlap_genes = length(intersect(hub_genes, SEX_genes)),
          total_hub = length(hub_genes), 
          total_Sex_DEGs = length(SEX_genes),
          total_background = length(background), 
          odds_ratio = NA, 
          p_value = NA
        ))
      }
      
      test <- fisher.test(tbl)
      
      tibble(
        Sex_celltype = Sex,
        overlap_genes = length(intersect(hub_genes, SEX_genes)),
        total_hub = length(hub_genes),
        total_Sex_DEGs = length(SEX_genes),
        total_background = length(background),
        odds_ratio = test$estimate, 
        p_value = test$p.value
      )
    })
  
  results %>% mutate(FDR = p.adjust(p_value, "fdr"), Sig = FDR < 0.05)
}

#Fishers exact 
fishers_exact_test = function(background_genes, gene_list, gene_set){
  background = data.frame("bg"= background_genes)
  background$Genes2Test <- background$bg  %in% gene_list
  
  CCA1_ppi <- data.frame(background$bg, background$Genes2Test)
  colnames(CCA1_ppi) <- c("gene","Genes2Test")
  
  CCA1_ppi$GeneSet = CCA1_ppi$gene %in% gene_set
  
  # Run Fisher's test for Sex top 1 percent
  tab <- table(CCA1_ppi$GeneSet, CCA1_ppi$Genes2Test)
  fishers <- fisher.test(tab)
  return(fishers)
}

summarise_hub_dge <- function(gene_col, dir_col, ppi_label) {
  
  asd_sub <- asd_dge %>%
    dplyr::filter(gene %in% hub_genes[[gene_col]]) %>%
    dplyr::mutate(
      Direction          = hub_genes[[dir_col]][match(gene, hub_genes[[gene_col]])],
      brainregion_shorten = "Frontal Cortex"
    ) %>%
    dplyr::group_by(Broad_Celltype, brainregion_shorten, Direction) %>%
    dplyr::summarise(mean_logFC = mean(logFC_ASD_vs_CTL, na.rm = TRUE), .groups = "drop")
  
  sex_sub <- sex_dge %>%
    dplyr::filter(gene %in% hub_genes[[gene_col]]) %>%
    dplyr::mutate(Direction = hub_genes[[dir_col]][match(gene, hub_genes[[gene_col]])]) %>%
    dplyr::group_by(Broad_Celltype, brainregion_shorten, Direction) %>%
    dplyr::summarise(mean_logFC = mean(logFC, na.rm = TRUE), .groups = "drop")
  
  dplyr::bind_rows(asd_sub, sex_sub) %>%
    dplyr::mutate(PPI = ppi_label)
}


enrichGO_function_single <- function(select_genes, bg) {
  entrez_select<- mapIds(org.Hs.eg.db, keys = select_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_bg     <- mapIds(org.Hs.eg.db, keys = bg,     column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  
  enrch     <- enrichGO(entrez_select, OrgDb = "org.Hs.eg.db", pvalueCutoff = 0.05,
                        pAdjustMethod = "fdr", universe = entrez_bg,
                        minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,ont = "BP")
  
  sig_enrch <- enrch@result[enrch@result$p.adjust < 0.05, ]
  
  list(Full_Enrich = enrch, Sig_Enrich = sig_enrch)
}

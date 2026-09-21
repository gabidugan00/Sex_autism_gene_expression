###-------------------------------------
#Running permutations on sex DGE across all celltypes and brain regions. 
#Shuffling metadata first across permutations
###------------------------------------

#This code calculations 1000 permutations of limma-voom results for sex-differences.

#Install packages
suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(ggplot2)
  library(corrplot)
  library(pheatmap)
  library(edgeR)
  library(limma)
  library(tidyverse)
  library(foreach)
  library(doParallel)
  library(edgeR)
  library(limma)
  library(dplyr)
})

####################################################
#Read in data to set up inputs for permutations
####################################################

#Read in pseudobulked object (Cellclass.subjectID)
alex_pb_cellclass <- readRDS("/data/DNU/alex/snRNAseq/alevin-all/h5ad/all_class_pseudo_corrected.rds")

#Extract cell type nams 
celltype <- sapply(colnames(alex_pb_cellclass), function(df) {
  parts <- unlist(strsplit(df, "\\."))
  if (length(parts) > 1) parts[2]   
})

#Extract sampleIDs
sample.id <- sapply(colnames(alex_pb_cellclass), function(df) {
  parts <- unlist(strsplit(df, "\\."))
  if (length(parts) > 1) parts[4]   
})
sample.id <- paste("sample",sample.id,sep ="_")

#Rename the column names to the celltype & sampleID
colnames(alex_pb_cellclass) <- paste(celltype,sample.id, sep =".") 

###############################################
#Set up inputs for limmma-voom permutations
###############################################

#Metadata
meta <- readRDS('/data/DNU/alex/snRNAseq/alevin-all/h5ad/meta_merged.rds')
colnames(meta)[which(colnames(meta) == 'BEC.arterial')] <- 'arterial'
colnames(meta)[which(colnames(meta) == 'BEC.venous')] <- 'venous'
colnames(meta)[which(colnames(meta) == 'BEC.cap')] <- 'cap'
meta$new.id <- paste("sample", meta$new.id, sep="_")
meta$add.sample.origident <- paste("sample", meta$orig.ident, sep="_")
meta <- meta[meta$repeat. != "repeat",]

# Get unique cell types to iterate through
celltypes <- sub("\\..*", "",colnames(alex_pb_cellclass))
celltypes_unique <- unique(celltypes)
celltypes_unique <- celltypes_unique[!grepl("SM|Tcell|Fibro|L5ET|Peri",celltypes_unique)]

sampleID <- sub(".*\\.", "",colnames(alex_pb_cellclass))
Brainregion <- meta$region[match(sampleID,meta$new.id)] 
brainregions <- unique(Brainregion)
meta_br <- split(meta,meta$region)

# pb_counts <- list()
# pb_counts_br <- list()
# d0<- list()
# for(celltype in unique(celltypes_unique)){for(brainregion in unique(brainregions)){
#   # Pseudobulk processing for each cell type and brain region
#   pb_counts[[celltype]] <- alex_pb_cellclass[, grepl(paste(celltype), colnames(alex_pb_cellclass))]
#   colnames(pb_counts[[celltype]]) <- sub(paste(celltype, ".", sep = ""), "", colnames(pb_counts[[celltype]]))
#   pb_counts_br[[celltype]][[brainregion]] <- pb_counts[[celltype]][, colnames(pb_counts[[celltype]]) %in% meta_br[[brainregion]]$new.id]
# 
#   # Filter samples with <30 cells
#   c <- meta[, c('new.id', celltype)]
#   colnames(c)[2] <- 'count'
#   keep_samples <- subset(c, count >= 30)$new.id
#   pb_counts_br[[celltype]][[brainregion]] <- pb_counts_br[[celltype]][[brainregion]][, colnames(pb_counts_br[[celltype]][[brainregion]]) %in% keep_samples]
# 
#   if(length(pb_counts_br[[celltype]][[brainregion]])>1){
#     # Remove mitochondrial and ribosomal genes
#     ribo <- rownames(pb_counts_br[[celltype]][[brainregion]])[grep("^(RPL|RPS|MRPL|MRPS)", rownames(pb_counts_br[[celltype]][[brainregion]]))]
#     ribo <- setdiff(ribo, c("RPS4X", "RPS4Y1", "RPS4Y2"))
#     mt <- rownames(pb_counts_br[[celltype]][[brainregion]])[grep("^MT-", rownames(pb_counts_br[[celltype]][[brainregion]]))]
#     pb_counts_br[[celltype]][[brainregion]] <- pb_counts_br[[celltype]][[brainregion]][!rownames(pb_counts_br[[celltype]][[brainregion]]) %in% c(ribo, mt), ]
# 
#     # Filter for highly expressed genes
#     cpm_cutoff <- 5
#     cpm_means <- rowMeans(edgeR::cpm(pb_counts_br[[celltype]][[brainregion]]))
#     keep_genes <- names(cpm_means[cpm_means > cpm_cutoff])
#     d0[[celltype]][[brainregion]] <- DGEList(pb_counts_br[[celltype]][[brainregion]][keep_genes, ])
#     d0[[celltype]][[brainregion]] <- calcNormFactors(d0[[celltype]][[brainregion]], method = 'TMM')
#   }
# }}
# 
# saveRDS(d0,"/gpfs/gsfs12/users/HGBomics/GabiDugan/OxCam/Human.Transcriptomics/Sex_fulldimensions_limma_voom/d0")

#################################################
#Limma-voom permutation pipeline 
#################################################
d0<- readRDS("/gpfs/gsfs12/users/HGBomics/GabiDugan/OxCam/Human.Transcriptomics/Sex_fulldimensions_limma_voom/d0")

shuffled_meta_perm<- readRDS("/data/HGBomics/GabiDugan/OxCam/Human.Transcriptomics/Sex_fulldimensions_limma_voom/ShuffleIndividuals_Permutation/permuted_metadata")

setwd("/vf/users/DNU/gabidugan/Sex_Permutations_DGE_114_Redo/")

# Define the number of permutations
num_permutations <- 1000

# Set up the cluster and register it as the parallel backend
n_cores <- min(35, parallel::detectCores() - 1)

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Verify backend is registered
getDoParWorkers() 
getDoParName()

all_files <- list.files(pattern = "_perm_")
file_paths <- file.path("~/Shuffle_114_Permutations/", all_files)

file_info <- file.info(file_paths)

valid_files <- all_files[file_info$size > 0]

standardize_filename <- function(fname) {
  base <- sub(".rds$", "", fname)
  parts <- strsplit(base, "_perm_")[[1]]
  name_part <- parts[1]
  perm_part <- parts[2]
  name_split <- strsplit(name_part, "_")[[1]]
  if (grepl("^[A-Z0-9]+$", name_split[1])) {
    celltype <- name_split[1]
    brainregion <- paste(name_split[-1], collapse = "_")
  } else {
    brainregion <- paste(name_split[-1], collapse = "_")
    celltype <- name_split[1]
  }
  
  paste0(brainregion, "_", celltype, "_perm_", perm_part, ".rds")
}

standardized_files <- sapply(valid_files, standardize_filename)

# Define output directories
output_dir <- "/vf/users/DNU/gabidugan/Sex_Permutations_DGE_114_Redo"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

setwd("/vf/users/DNU/gabidugan/Sex_Permutations_DGE_114_Redo/")

# Create list for results 
permres <- vector("list", 114)

# Generate all combinations to name the results with 
combinations <- unique(expand.grid(celltypes_unique,brainregions))

# Paste each combination
result <- apply(combinations, 1, function(row) paste(row[1], row[2], sep ="_"))
names(permres) <- result
top_results <- list()
all_perm_results<- list()
# Parallelize across cell types and brain regions
foreach(celltype = celltypes_unique, .packages = c("edgeR", "limma", "dplyr")) %:%
  foreach(brainregion = brainregions, .combine = 'c') %dopar% {
    
    tryCatch({
      print(paste("Processing:", celltype, "-", brainregion))
      missing_perms <- which(!(paste0(brainregion, "_", celltype, "_perm_", 1:1000, ".rds") %in% standardized_files))
      
      for (k in missing_perms) {
        #Filter out samples with fewer 30 cells 
        meta_br_filtered_shuffle <- shuffled_meta_perm[[k]][shuffled_meta_perm[[k]]$add.sample.origident %in% colnames(d0[[celltype]][[brainregion]]),]
        
        #Filter out permutations that do not have 5 individuals per sex
        f <- length(meta_br_filtered_shuffle$Shuffled_Sex[meta_br_filtered_shuffle$Shuffled_Sex=="Female"])
        
        m <- length(meta_br_filtered_shuffle$Shuffled_Sex[meta_br_filtered_shuffle$Shuffled_Sex=="Male"])
        
        
        if (min(c(f,m)) > 5) {
          
          design <- model.matrix(~0 + Shuffled_Sex + age + PMI + RIN + nreads + total_cells, data=meta_br_filtered_shuffle, na.action="na.pass")
          
          if(dim(design)[1] <= dim(design)[2]) {print('not enough individuals')} else {
            
            #Limma-voom
            y <- voomWithQualityWeights(d0[[celltype]][[brainregion]], design)
            
            fit <- lmFit(y, design)
            
            contrast <- makeContrasts(Shuffled_SexMale - Shuffled_SexFemale, levels=colnames(coef(fit)))
            fit <- contrasts.fit(fit, contrast)
            fit <- eBayes(fit)
            
            # Extract results and add permutation number
            top_results[[celltype]][[brainregion]] <- topTable(fit, sort.by = "P", n = Inf)
            
            top_results[[celltype]][[brainregion]] <- data.frame(top_results[[celltype]][[brainregion]])
            
            top_results[[celltype]][[brainregion]]$gene <- rownames(top_results[[celltype]][[brainregion]])
            
            top_results[[celltype]][[brainregion]]$permutation_number <- k
            
            top_results[[celltype]][[brainregion]]$Celltype_BrainRegion <- paste(celltype, brainregion, sep = ".")
            
            # Append to the list for this combination
            permres[[paste(celltype,brainregion, sep ="_")]][[paste0("perm_",k)]] <- top_results[[celltype]][[brainregion]]
            
            # Saave out each permutation result to read in
            saveRDS(permres[[paste(celltype,brainregion, sep ="_")]][[paste0("perm_",k)]],
                    file = file.path(output_dir, paste0(celltype, "_", brainregion, "_perm_", k, ".rds")))
          }
        }
      }
    })
  }

# Stop the parallel backend
stopCluster(cl)
registerDoSEQ()





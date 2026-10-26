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

# Read in pseudobulked object (Cellclass.subjectID) from sex-differences dataset
alex_pb_cellclass <- readRDS("~/all_class_pseudo_corrected.rds") 

#Extract cell type nams 
celltype <- sapply(colnames(alex_pb_cellclass), function(df) {
  parts <- unlist(strsplit(df, "\\."))
  if (length(parts) > 1) parts[2]   # Extract the part between the first and second periods (celltype)
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
meta <- readRDS('~/meta_merged.rds') #metadata from sex-differences dataset
colnames(meta)[which(colnames(meta) == 'BEC.arterial')] <- 'arterial'
colnames(meta)[which(colnames(meta) == 'BEC.venous')] <- 'venous'
colnames(meta)[which(colnames(meta) == 'BEC.cap')] <- 'cap'
meta$new.id <- paste("sample", meta$new.id, sep="_")
meta$add.sample.origident <- paste("sample", meta$orig.ident, sep="_")
meta <- meta[meta$repeat. != "repeat",]

# Get unique cell types to iterate through
celltypes <- sub("\\..*", "",colnames(alex_pb_cellclass))
celltypes_unique <- unique(celltypes)
celltypes_unique <- celltypes_unique[!grepl("SM|Tcell|Fibro|L5ET",celltypes_unique)]

sampleID <- sub(".*\\.", "",colnames(alex_pb_cellclass))
Brainregion <- meta$region[match(sampleID,meta$new.id)] #get brain regions for respective sample IDs from metadata
brainregions <- unique(Brainregion)
meta_br <- split(meta,meta$region)

# #Shuffle the sex label within each brain region and celltype for 1000 permutations
shuffled_meta_perm<- lapply(1:1000, function(i) {
   shuffled_meta_perm <-meta %>%
     mutate(
       RIN = scale(RIN),
       PMI = scale(PMI),
       nreads = scale(nreads),
       total_cells = scale(total_cells),
       age = scale(age),
       sex = sample(sex) #permute sex label separately for each permutation
     )
   return(shuffled_meta_perm)
 })
names(shuffled_meta_perm) <- seq_along(1:1000)
 
# shuffled_meta_perm<- readRDS("~/ShuffleFirst_Permutation/shuffled_meta_perm") #save out copy of shuffled meta for future runs, if needed

shuffled_meta_perm <- lapply(seq_along(shuffled_meta_perm), function(i) {
  shuffled_meta_perm[[i]]$Permutation <- i
  shuffled_meta_perm[[i]]
})

combine_shuffle_meta_perm <- Reduce(rbind,shuffled_meta_perm)

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
# saveRDS(d0,"~/d0")

#################################################
#Limma-voom permutation pipeline 
#################################################
d0<- readRDS("~/d0")

# Define the number of permutations
num_permutations <- 1000

# Set up the cluster and register it as the parallel backend
n_cores <- min(35, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Verify backend is registered
getDoParWorkers() 
getDoParName()

setwd("~/ShuffleFirst_Permutations_114BRCT/")

# Define output directories
output_dir <- "~/ShuffleFirst_Permutations_114BRCT/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Create list for results 
permres <- vector("list", 120)

#Remove sm and tcells from celltypes
celltypes <- celltypes[!grepl("SM|Tcell|Fibro|L5ET",celltypes)]

# Generate all combinations to name the results with 
combinations <- unique(expand.grid(brainregions, celltypes))

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
      for (k in seq_len(num_permutations)) {
        #Reduce  metadata to only kept samples for brainregion and for samples with >30 cells for celltype (already done for d0)
        meta_br_filtered_shuffle <- shuffled_meta_perm[[k]][shuffled_meta_perm[[k]]$add.sample.origident %in% colnames(d0[[celltype]][[brainregion]]),]
        
        #Get dimensions of each sex
        f <- length(meta_br_filtered_shuffle$sex[meta_br_filtered_shuffle$sex=="Female"])
        m <- length(meta_br_filtered_shuffle$sex[meta_br_filtered_shuffle$sex=="Male"])
        
        # Check if there are enough individuals for both sexes
        if (min(c(f,m)) > 5) {
          design <- model.matrix(~0 + sex + age + PMI + RIN + nreads + total_cells, data=meta_br_filtered_shuffle, na.action="na.pass")
          # Check dimensions degrees of freedom before proceeding
          if(dim(design)[1] <= dim(design)[2]) {print('not enough individuals')} else {
            #Limma-voom
            y <- voomWithQualityWeights(d0[[celltype]][[brainregion]], design)
            fit <- lmFit(y, design)
            contrast <- makeContrasts(sexMale - sexFemale, levels=colnames(coef(fit)))
            fit <- contrasts.fit(fit, contrast)
            fit <- eBayes(fit)
            
            # Extract results and add permutation number
            top_results[[celltype]][[brainregion]] <- topTable(fit, sort.by = "P", n = Inf)
            top_results[[celltype]][[brainregion]] <- data.frame(top_results[[celltype]][[brainregion]])
            top_results[[celltype]][[brainregion]]$gene <- rownames(top_results[[celltype]][[brainregion]])
            top_results[[celltype]][[brainregion]]$permutation_number <- k
            top_results[[celltype]][[brainregion]]$Celltype_BrainRegion <- paste(celltype, brainregion, sep = ".")
            
            # Append to the list for this combination
            permres[[paste(celltype,brainregion, sep ="_")]][[k]] <- top_results[[celltype]][[brainregion]]
            
            # Save each permutation result individually (optional)
            saveRDS(permres[[paste(celltype,brainregion, sep ="_")]][[k]],
                    file = file.path(output_dir, paste0(celltype, "_", brainregion, "_perm_", k, ".rds")))
          }
        }
      }
    })
  }


# Stop the parallel backend
stopCluster(cl)
registerDoSEQ()


library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(tidyverse)

# Load the merged Seurat object after initial filtering
merged_seurat_filtered = readRDS("../Desktop/Rprojects/scRNA-seq/toturials project/my project/preprocessing_1/merged_seurat_filtered")

# View metadata to check loaded object info
view(merged_seurat_filtered@meta.data)

# Split the merged object by sample to process each sample independently
sample_list <- SplitObject(merged_seurat_filtered, split.by = "sample_id")


# Loop through each sample for DoubletFinder workflow
for (i in 1:length(sample_list)) {
  seurat_obj <- sample_list[[i]]
  
  # Step 1: Normalize data and find variable features per sample
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  
  # Step 2: Scale data and run PCA for dimensionality reduction
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  
  # Step 3: Run neighborhood graph and clustering to prepare for doublet detection
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
  
  # Step 4: Run parameter sweep to find optimal pK for doublet detection
  sweep.res.list <- paramSweep(seurat_obj, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # Optional: Inspect and select best pK value (peak of BCmetric)
  best.pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))
  
  # Step 5: Estimate expected number of doublets (7.5% of cells)
  nExp_poi <- round(0.075 * ncol(seurat_obj))
  
  # Run DoubletFinder with optimized parameters
  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:20, pN = 0.25, 
                              pK = best.pK , nExp =  nExp_poi, sct = FALSE)
  
  # Step 6: Save the updated Seurat object back into the list
  sample_list[[i]] <- seurat_obj
}

# Check metadata of the second sample after DoubletFinder step
head(sample_list[[2]]@meta.data)

###########################################################
# Filter out predicted doublets from each sample

for (i in 1:length(sample_list)) {
  seurat_obj <- sample_list[[i]]
  
  # Identify the column containing DoubletFinder classifications
  df_col <- grep("DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
  
  # Assign doublet classification to a standard column 'doublet_status'
  seurat_obj$doublet_status <- seurat_obj@meta.data[[df_col]]
  
  # Keep only singlet cells
  seurat_obj <- subset(seurat_obj, subset = doublet_status == "Singlet")
  
  # Save filtered object back into the list
  sample_list[[i]] <- seurat_obj
}

# Merge filtered samples back into one Seurat object
merged_filtered <- merge(sample_list[[1]], y = sample_list[2:length(sample_list)])

# Save the final merged and filtered Seurat object for downstream analysis
saveRDS(merged_filtered, 
        file = "../Desktop/Rprojects/scRNA-seq/toturials project/my project/DoubletFinder/merged_filtered.rds")

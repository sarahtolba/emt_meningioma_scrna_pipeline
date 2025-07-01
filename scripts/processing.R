# Load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# Set base path for your sample folders
base_path = "../Desktop/Rprojects/scRNA-seq/toturials project/tutorials/interation/files/"

# List all subdirectories (samples)
loc <- list.dirs(path = base_path , recursive = FALSE, full.names = FALSE)

# Loop through each sample folder and read data
for (i in loc) {
  folder_path <- file.path(base_path, i)
  
  matrix_file <- list.files(folder_path, pattern = "matrix.mtx$", full.names = TRUE)
  features_file <- list.files(folder_path, pattern = "features.tsv$", full.names = TRUE)
  barcodes_file <- list.files(folder_path, pattern = "barcodes.tsv$", full.names = TRUE)
  
  sample_name <- i
  
  # Read count matrix (genes x cells)
  cts <- ReadMtx(mtx = matrix_file,
                 features = features_file,
                 cells = barcodes_file)
  
  # Create Seurat object for each sample with the sample name
  assign(sample_name, CreateSeuratObject(counts = cts))
  
}

# Merge all Seurat objects into one combined object
merged_seurat <- merge(
  x = MSC1,
  y = c(MSC2,MSC3,MSC4,MSC4_Dura,MSC5,MSC5_BTI,MSC5_Dura,MSC6,MSC6_BTI),
  add.cell.ids = c("MSC1_tumor","MSC2_tumor","MSC3_tumor","MSC4_tumor","MSC4_Dura",
                   "MSC5_tumor","MSC5_BTI", "MSC5_Dura","MSC6_tumor","MSC6_BTI")
)

# Quick QC: View metadata table
view(merged_seurat@meta.data)

# Extract sample info by splitting cell barcodes metadata
merged_seurat@meta.data$sample <- rownames(merged_seurat@meta.data)
merged_seurat@meta.data <- separate(
  merged_seurat@meta.data,
  col = "sample",
  into = c("patient","type", "Barcode"),
  sep = "_" 
)

# Check unique sample types
unique(merged_seurat@meta.data$type)

# Calculate mitochondrial gene percentage per cell
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')

# Visualize QC metrics: nFeature_RNA, nCount_RNA, mitoPercent
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3) 

# Scatter plot of counts vs features with smoothing to check linearity
FeatureScatter(merged_seurat , feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# Filter cells based on quality metrics
merged_seurat_filtered <- subset(merged_seurat,
                                 subset = nFeature_RNA > 500 & nFeature_RNA < 7500 &
                                   nCount_RNA > 1000 &
                                   mitoPercent < 10
)

# Check filtered object summary
merged_seurat_filtered
merged_seurat

# Scatter plot for filtered data
FeatureScatter(merged_seurat_filtered , feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

# Standard workflow to assess batch effect before integration

# Normalize data
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)

# Identify variable features (genes)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)

# Scale data
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)

# PCA dimensionality reduction
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)

# View top PCA features and heatmap for first component
print(merged_seurat_filtered[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(merged_seurat_filtered, dims = 1, cells = 500, balanced = TRUE)

# Determine optimal number of PCs with elbow plot
ElbowPlot(merged_seurat_filtered)

# Construct nearest neighbor graph
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)

# Cluster cells
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)

# Run UMAP for visualization
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)

# Visualize batch effect or sample type on UMAP
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'type')
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'patient')

# Arrange plots side by side
grid.arrange(p1, p2, ncol = 2, nrow = 1)

# Additional preprocessing notes (to be done later):
# - Split merged object by sample
# - Normalize and find variable features individually
# - Run DoubletFinder per sample
# - Re-normalize and scale data after doublet removal
# - Integrate cleaned samples for downstream analysis

# Add sample_id column based on barcode pattern
head(merged_seurat_filtered@meta.data)
merged_seurat_filtered$sample_id <- gsub("_[^_]+$", "", rownames(merged_seurat_filtered@meta.data))
view(merged_seurat_filtered@meta.data)

# Save filtered and preprocessed Seurat object for later use
saveRDS(merged_seurat_filtered, 
        file = "../Desktop/Rprojects/scRNA-seq/toturials project/tutorials/integration/merged_seurat_filtered")

# Important notes about DoubletFinder & integration:
# 1. DoubletFinder works best on individually preprocessed samples.
# 2. Requires tuning parameters pN, pK, and nExp.
# 3. Run after normalization, variable feature selection, scaling, PCA.
# 4. Remove predicted doublets before integration.
# 5. Re-normalize and scale after doublet removal.
# 6. Then perform integration to minimize batch effects.


# Load required libraries
library(tidyverse)
library(Seurat)
library(multtest)
library(metap)
library(ggsci)

# Load the integrated Seurat object and precomputed markers table
seurat.integrated = readRDS("../Desktop/new_seu.rds")
markers = read.csv("../Desktop/markers")

# Optional: View metadata to understand the data structure
view(seurat.integrated@meta.data)


# ====== Data Scaling, PCA, UMAP and Clustering ======

# Normalize and scale data, then run PCA for dimensionality reduction
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)

# Visualize variance explained by PCs to choose number of PCs
ElbowPlot(seurat.integrated)

# Run UMAP using first 20 PCs for 2D visualization
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:20)

# Find neighbors and clusters at multiple resolutions to explore cluster granularity
seurat.integrated = FindNeighbors(seurat.integrated,reduction = "pca", dims = 1:20) 
seurat.integrated = FindClusters(seurat.integrated,resolution = c(0.1, 0.3, 0.5, 0.7, 1))

# View metadata after clustering
head(seurat.integrated@meta.data)


# ====== Visualize Clusters and Sample Types ======

p1 = DimPlot(seurat.integrated, group.by = "integrated_snn_res.0.1", label = TRUE) +
  scale_color_npg()

p2 = DimPlot(seurat.integrated, group.by = "type", label = TRUE) +
  scale_color_npg()

p1 + p2  # Display side-by-side for comparison


# ====== Set Cluster Identity for Downstream Analysis ======

Idents(seurat.integrated) <- "integrated_snn_res.0.1"  # Choose resolution 0.1 clusters
table(Idents(seurat.integrated), seurat.integrated$type)  # Table of clusters vs sample types


# ====== Prepare RNA Assay for Analysis ======

# Set default assay to RNA to work with expression data
DefaultAssay(seurat.integrated) <- "RNA"

# Join layers if working with multi-layer data (optional step)
seurat.integrated <- JoinLayers(seurat.integrated)


# ====== EMT Gene Expression Visualization ======

# Define canonical EMT marker genes
emt_genes <- c("VIM", "FN1", "ZEB1", "SNAI1", "SNAI2", "CDH2", "TWIST1")

# DotPlot for EMT genes across clusters or samples
DotPlot(seurat.integrated, features = emt_genes) + RotatedAxis()

# FeaturePlot to visualize spatial expression on UMAP, thresholding out low expression
FeaturePlot(seurat.integrated, features = emt_genes, min.cutoff = 'q10')


# ====== Calculate EMT Scores per Cell ======

# Wrap genes in a list as required by AddModuleScore()
emt_genes <- list(c("VIM", "FN1", "ZEB1", "SNAI1", "SNAI2", "CDH2", "TWIST1"))

# Add EMT module score (average expression of gene set adjusted against control gene sets)
seurat.integrated <- AddModuleScore(
  object = seurat.integrated,
  features = emt_genes,
  name = "EMT"
)

# View new EMT scores added as "EMT1" column in metadata
view(seurat.integrated@meta.data)
head(seurat.integrated@meta.data)

# Visualize EMT scores on UMAP and violin plots grouped by clusters and sample types
FeaturePlot(seurat.integrated, features = "EMT1", min.cutoff = "q10", max.cutoff = "q90")

VlnPlot(seurat.integrated, features = "EMT1", group.by = "integrated_snn_res.0.1", pt.size = 0) +
  scale_fill_npg()

VlnPlot(seurat.integrated, features = "EMT1", group.by = "type", pt.size = 0) + 
  scale_fill_npg()


# ====== Classify Cells as EMT_high or EMT_low ======

seurat.integrated$EMT_status <- ifelse(seurat.integrated$EMT1 > 0.5, "EMT_high", "EMT_low")


# ====== Visualize Top Marker Genes from Differential Expression ======

top_genes <- rownames(emt_markers)[1:10]  # First 10 top markers

# Violin plot and dot plot by EMT_status for top genes
VlnPlot(seurat.integrated, features = top_genes, group.by = "EMT_status", pt.size = 0)
DotPlot(seurat.integrated, features = top_genes, group.by = "EMT_status")



# ====== Differential Expression: EMT_high vs EMT_low ======

Idents(seurat.integrated) <- "EMT_status"

markers = FindMarkers(seurat.integrated, ident.1 = "EMT_high", ident.2 = "EMT_low")

# Plot EMT1 score and EMT_status UMAPs for visual check
VlnPlot(seurat.integrated, features = "EMT1", group.by = "EMT_status")
DimPlot(seurat.integrated, group.by = "EMT_status")

# Example FeaturePlot for one key gene (FN1)
FeaturePlot(seurat.integrated, features = c('FN1'), min.cutoff = 'q10')


# ====== Filter markers by significance and effect size ======

markers_filtered <- markers %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
row.names(markers_filtered) = markers_filtered$X


# ====== Visualize filtered top markers ======

top10 <- rownames(markers)[1:10]

VlnPlot(seurat.integrated, features = top10, group.by = "EMT_status", pt.size = 0) + scale_fill_npg()
DotPlot(seurat.integrated, features = top10, group.by = "EMT_status") + RotatedAxis() 
FeaturePlot(seurat.integrated, features = top10) +  scale_fill_npg()

# Heatmap of top markers
heatmap_plot <- DoHeatmap(seu, features = top10) + 
  ggtitle("Top Marker Genes Heatmap")


# ====== Save DE results ======

write.csv(markers_filtered,file = "results/markers_filtered.csv")
write.csv(markers, file = "results/markers.csv", row.names = TRUE)

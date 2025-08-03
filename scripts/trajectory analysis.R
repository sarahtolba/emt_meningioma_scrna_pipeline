# Load required libraries
library(monocle3)
library(Seurat)
library(tidyverse)
library(SeuratWrappers)
library(ggplot2)

# Load data
seurat_obj <- readRDS("../Desktop/emt_meningioma_scrna_pipeline/seurat_objects/singler.rds")
markers <- read.csv("../Desktop/emt_meningioma_scrna_pipeline/results/markers.csv")
row.names(markers) <- markers$X
markers <- markers[, -1]  # gene metadata

# Convert Seurat object to Monocle CellDataSet
cds <- as.cell_data_set(seurat_obj)

# Preprocess CDS (PCA)
cds <- preprocess_cds(cds, num_dim = 50)

# Align Seurat UMAP to Monocle
cds@int_colData@listData[["reducedDims"]]$UMAP <- seurat_obj@reductions$umap@cell.embeddings

# Cluster cells
cds <- cluster_cells(cds, reduction_method = "UMAP")

# Learn principal graph
cds <- learn_graph(cds)

# Define root cell from EMT-low cells
root_cell <- colnames(seurat_obj)[seurat_obj$EMT_status == "EMT_low"][1]

# Order cells in pseudotime
cds <- order_cells(cds, root_cells = root_cell)

# Plot pseudotime trajectory
pdf("../Desktop/emt_meningioma_scrna_pipeline/results/PseudotimeTrajectory.pdf", width = 8, height = 6)
plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)
dev.off()

# Add pseudotime to Seurat object
seurat_obj$pseudotime <- pseudotime(cds)

# Save updated Seurat object
saveRDS(seurat_obj, "../Desktop/emt_meningioma_scrna_pipeline/seurat_objects/pseudo.rds")

# Add gene names
rowData(cds)$gene_short_name <- rownames(cds)

# Plot gene expression over pseudotime (example genes)
pdf("../Desktop/emt_meningioma_scrna_pipeline/results/PseudotimeTrajectory1.pdf", width = 10, height = 5)
plot_genes_in_pseudotime(
  cds[c("VIM", "CDH1", "ZEB1", "SNAI2"), ],
  color_cells_by = "pseudotime"
)
dev.off()

# Identify genes changing along pseudotime
deg_pseudo <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
write.csv(deg_pseudo, "../Desktop/emt_meningioma_scrna_pipeline/results/dynamic_genes_along_pseudotime.csv", row.names = TRUE)

# Filter significant genes
deg_sig <- deg_pseudo %>% filter(q_value < 0.05)
write.csv(deg_sig, "../Desktop/emt_meningioma_scrna_pipeline/results/significant_dynamic_genes.csv", row.names = TRUE)

# Replace Inf pseudotime values with NA
seurat_obj$pseudotime <- pseudotime(cds)
seurat_obj$pseudotime[!is.finite(seurat_obj$pseudotime)] <- NA

# Summary of pseudotime values
head(seurat_obj$pseudotime)
summary(seurat_obj$pseudotime)
any(!is.finite(seurat_obj$pseudotime))

# Visualizations
FeaturePlot(seurat_obj, features = "pseudotime", reduction = "umap")
VlnPlot(seurat_obj, features = "pseudotime", group.by = "EMT_status")

# Top genes changing along pseudotime
top_genes <- deg_sig %>% arrange(q_value) %>% slice(1:50) %>% pull(gene_short_name)
plot_genes_in_pseudotime(cds[top_genes, ])

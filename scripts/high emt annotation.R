# Cell type annotation using SingleR and reference dataset
# SingleR annotates cells based on gene expression similarity to reference profiles
# Dimensionality reduction is NOT required for SingleR, but QC steps beforehand are recommended

library(Seurat)
library(SingleR)
library(celldex) # for reference datasets
library(ggsci)
library(ggplot2)

# Load Seurat object previously processed
seu = readRDS("../Desktop/Rprojects/scRNA-seq/toturials project/my project/find cluster marker/emt markers/new_seu.rds")
head(seu@meta.data)

# Set default assay to RNA counts (recommended for SingleR)
DefaultAssay(seu) = "RNA"

# Set cluster identities (e.g., resolution 0.1)
Idents(seu) = "integrated_snn_res.0.1"

# Visualize clusters on UMAP
DimPlot(seu, reduction = "umap") + scale_fill_npg()

# Load reference dataset (Human Primary Cell Atlas) for annotation
reference = celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(reference)))

# Extract raw counts matrix from Seurat object (required input for SingleR)
seu_counts <- GetAssayData(seu, slot = 'counts')
head(seu_counts)

# Annotate all cells by comparing test data to reference labels
annotate <- SingleR(test = seu_counts,
                    ref = reference,
                    labels = reference$label.main)

# Add SingleR labels to Seurat metadata, matching cell names
seu$singleR.labels <- annotate$labels[match(rownames(seu@meta.data), rownames(annotate))]
head(seu@meta.data)

# Visualize SingleR labels on UMAP
DimPlot(seu, reduction = 'umap', group.by = 'singleR.labels')

# Annotate only the subset of cells with high EMT score

# Subset high EMT cells from Seurat object
seu_high_emt <- subset(seu, subset = EMT_status == "EMT_high")

# Extract counts for high EMT cells
seu_high_counts <- GetAssayData(seu_high_emt, slot = 'counts')

# Run SingleR annotation on high EMT subset
annotate_high <- SingleR(test = seu_high_counts,
                         ref = reference,
                         labels = reference$label.main)

# Initialize a column in Seurat metadata for high EMT annotations (NA by default)
seu$SingleR_high_emt <- NA

# Match annotations back to full Seurat object cells (only high EMT cells get labels)
seu$SingleR_high_emt[colnames(seu_high_emt)] <- annotate_high$labels

# Visualize high EMT cell annotations on UMAP
DimPlot(seu, reduction = 'umap', group.by = 'SingleR_high_emt') +
  ggtitle("SingleR Annotation for High EMT Cells")

# Display contingency tables to compare annotations with EMT status
table(seu$SingleR_high_emt, seu$EMT_status)
table(seu$singleR.labels, seu$EMT_status)

# Save updated Seurat object with annotations for future use
library(readr)
write_rds(seu, "../Desktop/singler.rds", compress = "gz")

# Optionally reload saved Seurat object when needed
seu = readRDS("../Desktop/singler.rds")


# Diagnostic steps to evaluate SingleR predictions

setwd('/Users/SarahSamir/Desktop/')

library(Seurat)       # Seurat object handling
library(SingleR)      # Automated cell annotation
library(celldex)      # Reference datasets
library(pheatmap)     # Heatmaps for annotation evaluation

# Load previously annotated Seurat object
seu <- readRDS("../Desktop/singler.rds")

# Extract raw counts (required for diagnostic plots)
seu_counts <- GetAssayData(seu, slot = "counts")

# Load reference again
reference <- celldex::HumanPrimaryCellAtlasData()

# Run SingleR to obtain full diagnostic prediction object (includes confidence scores)
diagnostic_pred <- SingleR(test = seu_counts,
                           ref = reference,
                           labels = reference$label.main)

# Plot heatmap showing confidence scores for each label per cell
plotScoreHeatmap(diagnostic_pred)

# Plot delta distribution showing confidence gap between top and second-best labels
plotDeltaDistribution(diagnostic_pred)

# Generate confusion matrix comparing SingleR labels to Seurat clustering
tab <- table(Assigned = seu$singleR.labels, Clusters = Idents(seu))

# Plot confusion matrix as a heatmap (log-transformed for clarity)
pheatmap(log10(tab + 10),
         color = colorRampPalette(c("white", "blue"))(10),
         main = "SingleR vs Seurat Cluster Agreement")

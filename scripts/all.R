############################################################
# 1) LOAD LIBRARIES & INITIAL MERGE + QC
############################################################

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


############################################################
# 2) DOUBLETFINDER PER SAMPLE
############################################################

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


############################################################
# 3) INTEGRATION AFTER DOUBLETFINDER
############################################################

# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

#read our merged file
merged.filtered = readRDS("../Desktop/Rprojects/scRNA-seq/toturials project/my project/DoubletFinder/merged_filtered")
str(merged.filtered)

view(merged.filtered@meta.data) 

# perform integration to correct for batch effects ------
obj.list <- SplitObject( merged.filtered, split.by = 'patient')



for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

# Fix for Seurat v5: flatten the layers
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- JoinLayers(obj.list[[i]], assay = "RNA")
  DefaultAssay(obj.list[[i]]) <- "RNA"
}

# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)

view(seurat.integrated@meta.data)

# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
ElbowPlot(seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:20)

p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'patient')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'type')

p3+p4

saveRDS(seurat.integrated, file = "../Desktop/Rprojects/scRNA-seq/toturials project/my project/integration/seurat_integrated.rds")


############################################################
# 4) EMT SCORING, CLUSTERING & DE (EMT_high vs EMT_low)
############################################################

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


############################################################
# 5) CELL TYPE ANNOTATION WITH SINGLER + DIAGNOSTICS
############################################################

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


############################################################
# 6) PSEUDOTIME ANALYSIS WITH MONOCLE3
############################################################

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


############################################################
# GO ENRICHMENT SCRIPT (EMT_high vs EMT_low)
############################################################

# Load libraries
library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# ------------------------------------------
# 1. Load Seurat object
# ------------------------------------------

seurat.integrated <- readRDS("../seurat_objects/singler.rds")

table(seurat.integrated$EMT_status)

# ------------------------------------------
# 2. Differential Expression
# ------------------------------------------

Idents(seurat.integrated) <- "EMT_status"

# Fast DE using downsampling (recommended)
markers <- FindMarkers(
  seurat.integrated,
  ident.1 = "EMT_high",
  ident.2 = "EMT_low",
  logfc.threshold = 0,
  min.pct = 0,
  downsample = 2000      # makes DE very fast
)

# Save DE table
write.csv(markers, "../results/markers_for_GO.csv")

# ------------------------------------------
# 3. Extract upregulated genes in EMT_high
# ------------------------------------------

# Significant genes (adjust if needed)
markers_sig <- markers[markers$p_val_adj < 0.05, ]

# Upregulated only
up_genes <- rownames(markers_sig[markers_sig$avg_log2FC > 0, ])

cat("Number of upregulated genes:", length(up_genes), "\n")

# ------------------------------------------
# 4. GO Biological Process Enrichment
# ------------------------------------------

ego <- enrichGO(
  gene          = up_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",          # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# Convert to data frame
ego_df <- as.data.frame(ego)

# Save results
write.csv(ego_df, "../results/GO_EMT_high_enrichment.csv", row.names = FALSE)

# ------------------------------------------
# 5. Plots 
# ------------------------------------------

# Top 20 enriched processes
pdf("../results/GO_BP_dotplot_top20.pdf", width = 8, height = 6)
dotplot(ego, showCategory = 20) +
  ggtitle("Top 20 Enriched GO Biological Processes (EMT_high)")
dev.off()

pdf("../results/GO_BP_barplot_top20.pdf", width = 8, height = 6)
barplot(ego, showCategory = 20) +
  ggtitle("Top 20 GO Processes Upregulated in EMT_high Cells")
dev.off()




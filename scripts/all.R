# =========================================
# Single-Cell RNA-seq Analysis Pipeline
# EMT–Tumor Progression with Seurat + SingleR
# Author: You
# =========================================

# ========== Load Required Libraries ========== #
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(DoubletFinder)
library(SingleR)
library(celldex)
library(pheatmap)
library(ggsci)
library(multtest)
library(metap)

# ========== Step 1: Load Raw 10x Data ========== #
base_path <- "../Desktop/Rprojects/scRNA-seq/toturials project/tutorials/interation/files/"
loc <- list.dirs(path = base_path , recursive = FALSE, full.names = FALSE)

for (i in loc) {
  folder_path <- file.path(base_path, i)
  matrix_file <- list.files(folder_path, pattern = "matrix.mtx$", full.names = TRUE)
  features_file <- list.files(folder_path, pattern = "features.tsv$", full.names = TRUE)
  barcodes_file <- list.files(folder_path, pattern = "barcodes.tsv$", full.names = TRUE)
  sample_name <- i
  cts <- ReadMtx(mtx = matrix_file, features = features_file, cells = barcodes_file)
  assign(sample_name, CreateSeuratObject(counts = cts))
}

# ========== Step 2: Merge All Samples ========== #
merged_seurat <- merge(
  x = MSC1,
  y = c(MSC2,MSC3,MSC4,MSC4_Dura,MSC5,MSC5_BTI,MSC5_Dura,MSC6,MSC6_BTI),
  add.cell.ids = c("MSC1_tumor","MSC2_tumor","MSC3_tumor","MSC4_tumor","MSC4_Dura",
                   "MSC5_tumor","MSC5_BTI", "MSC5_Dura","MSC6_tumor","MSC6_BTI")
)

# ========== Step 3: Metadata Preparation & QC ========== #
merged_seurat@meta.data$sample <- rownames(merged_seurat@meta.data)
merged_seurat@meta.data <- separate(
  merged_seurat@meta.data,
  col = "sample",
  into = c("patient","type", "Barcode"),
  sep = "_" )

merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')

# Filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 &
                                   nCount_RNA > 1000 & mitoPercent < 10)

# Save filtered Seurat object
saveRDS(merged_seurat_filtered, file = "../Desktop/Rprojects/scRNA-seq/toturials project/tutorials/integration/merged_seurat_filtered")

# ========== Step 4: Run DoubletFinder on Each Sample ========== #
sample_list <- SplitObject(merged_seurat_filtered, split.by = "sample_id")

for (i in 1:length(sample_list)) {
  seurat_obj <- sample_list[[i]]
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
  sweep.res.list <- paramSweep(seurat_obj, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  best.pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))
  nExp_poi <- round(0.075 * ncol(seurat_obj))
  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:20, pN = 0.25, pK = best.pK, nExp = nExp_poi, sct = FALSE)
  sample_list[[i]] <- seurat_obj
}

for (i in 1:length(sample_list)) {
  seurat_obj <- sample_list[[i]]
  df_col <- grep("DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
  seurat_obj$doublet_status <- seurat_obj@meta.data[[df_col]]
  seurat_obj <- subset(seurat_obj, subset = doublet_status == "Singlet")
  sample_list[[i]] <- seurat_obj
}

merged_filtered <- merge(sample_list[[1]], y = sample_list[2:length(sample_list)])
saveRDS(merged_filtered, file = "../Desktop/Rprojects/scRNA-seq/toturials project/my project/DoubletFinder/merged_filtered")

# ========== Step 5: Integration with Seurat ========== #
obj.list <- SplitObject(merged_filtered, split.by = 'patient')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]])
  obj.list[[i]] <- JoinLayers(obj.list[[i]], assay = "RNA")
  DefaultAssay(obj.list[[i]]) <- "RNA"
}
features <- SelectIntegrationFeatures(object.list = obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
seurat.integrated <- IntegrateData(anchorset = anchors)

seurat.integrated <- ScaleData(seurat.integrated)
seurat.integrated <- RunPCA(seurat.integrated)
seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:20)

p1 <- DimPlot(seurat.integrated, group.by = 'patient')
p2 <- DimPlot(seurat.integrated, group.by = 'type')
p1 + p2
saveRDS(seurat.integrated, file = "../Desktop/Rprojects/scRNA-seq/toturials project/my project/integration/seurat_integrated.rds")

# ========== Step 6: EMT Gene Scoring ========== #
emt_genes <- list(c("VIM", "FN1", "ZEB1", "SNAI1", "SNAI2", "CDH2", "TWIST1"))
seurat.integrated <- AddModuleScore(object = seurat.integrated, features = emt_genes, name = "EMT")
seurat.integrated$EMT_status <- ifelse(seurat.integrated$EMT1 > 0.5, "EMT_high", "EMT_low")

# ========== Step 7: Differential Expression ========== #
Idents(seurat.integrated) <- "EMT_status"
emt_markers <- FindMarkers(seurat.integrated, ident.1 = "EMT_high", ident.2 = "EMT_low")
markers_filtered <- emt_markers %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
row.names(markers_filtered) = markers_filtered$X

# ========== Step 8: Cell Type Annotation with SingleR ========== #
reference <- celldex::HumanPrimaryCellAtlasData()
seu_counts <- GetAssayData(seurat.integrated, slot = 'counts')
annotate <- SingleR(test = seu_counts, ref = reference, labels = reference$label.main)
seurat.integrated$singleR.labels <- annotate$labels[match(rownames(seurat.integrated@meta.data), rownames(annotate))]

# Subset and reannotate high EMT cells
seu_high_emt <- subset(seurat.integrated, subset = EMT_status == "EMT_high")
seu_high_counts <- GetAssayData(seu_high_emt, slot = 'counts')
annotate_high <- SingleR(test = seu_high_counts, ref = reference, labels = reference$label.main)
seurat.integrated$SingleR_high_emt <- NA
seurat.integrated$SingleR_high_emt[colnames(seu_high_emt)] <- annotate_high$labels

# Diagnostic checks
diagnostic_pred <- SingleR(test = seu_counts, ref = reference, labels = reference$label.main)
plotScoreHeatmap(diagnostic_pred)
plotDeltaDistribution(diagnostic_pred)
tab <- table(Assigned = seurat.integrated$singleR.labels, Clusters = Idents(seurat.integrated))
pheatmap(log10(tab + 10), color = colorRampPalette(c("white", "blue"))(10),
         main = "SingleR vs Seurat Cluster Agreement")
saveRDS(diagnostic_pred, "../Desktop/singler_diagnostics.rds")

# ========== Final Note ========== #
# All steps from raw data to EMT labeling and cell annotation are included.
# Next steps: trajectory inference, cell-cell communication, pathway enrichment.
# Let me know to continue those.

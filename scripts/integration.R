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







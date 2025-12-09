############################################################
# GSEA SCRIPT (EMT_high vs EMT_low)

# Load libraries
library(Seurat)
library(tidyverse)
library(msigdbr)
library(fgsea)
library(ggplot2)

# ------------------------------------------
# 1. Load Seurat object (EDIT THIS PATH)
# ------------------------------------------

seurat.integrated <- readRDS("../seurat_objects/singler.rds")

# Check metadata to ensure EMT_status exists
table(seurat.integrated$EMT_status)

# ------------------------------------------
# 2. Differential Expression
# ------------------------------------------

Idents(seurat.integrated) <- "EMT_status"

markers <- read.csv("../results/markers.csv")
head(markers)


# Keep everything (no logFC threshold → needed for GSEA)
markers <- FindMarkers(
  seurat.integrated,
  ident.1 = "EMT_high",
  ident.2 = "EMT_low",
  logfc.threshold = 0,
  min.pct = 0,
  downsample = 2000
)


# Save markers (optional)
write.csv(markers, "../results/markers_for_GSEA.csv")

# ------------------------------------------
# 3. Prepare ranked gene list for GSEA
# ------------------------------------------

geneList <- markers$avg_log2FC
names(geneList) <- rownames(markers)
geneList <- sort(geneList, decreasing = TRUE)

# ------------------------------------------
# 4. Load Hallmark gene sets (MSigDB)
# ------------------------------------------

m_df <- msigdbr(species = "Homo sapiens", category = "H")

hallmark_sets <- m_df %>%
  select(gs_name, gene_symbol) %>%
  split(.$gs_name)

# ------------------------------------------
# 5. Run GSEA using fgsea
# ------------------------------------------

gsea_res <- fgsea(
  pathways = hallmark_sets,
  stats = geneList,
  minSize = 10,
  maxSize = 500
)

gsea_res <- gsea_res[order(gsea_res$padj), ]

# Save results
write.csv(gsea_res, "GSEA_results_Hallmark.csv", row.names = FALSE)

# ------------------------------------------
# 6. Simple barplot of top 10 pathways
# ------------------------------------------

top10 <- gsea_res %>% head(10)

ggplot(top10, aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 Enriched Hallmark Gene Sets (EMT_high vs EMT_low)",
    x = "Pathway",
    y = "NES"
  ) +
  theme_bw()

ggsave("GSEA_Top10_barplot.png", width = 8, height = 5)

# ------------------------------------------
# 7. Plot EMT enrichment curve (if exists)
# ------------------------------------------

if ("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" %in% names(hallmark_sets)) {
  pdf("GSEA_EMT_enrichment.pdf", width = 7, height = 5)
  plotEnrichment(
    hallmark_sets[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
    geneList
  ) + ggtitle("Enrichment: EMT Hallmark")
  dev.off()
}

############################################################
# END OF SCRIPT
############################################################

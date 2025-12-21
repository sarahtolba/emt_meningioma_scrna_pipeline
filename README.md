# 🧬 Single-Cell Characterization of EMT Programs in Human Meningioma

This project implements an end to end scRNA seq analysis workflow in R using GEO dataset GSE183655. The pipeline includes quality control filtering, per sample doublet detection with DoubletFinder, Seurat integration and clustering, EMT module scoring and EMT high vs EMT low differential expression, SingleR cell type annotation, Monocle3 pseudotime trajectory analysis, and GO Biological Process enrichment with clusterProfiler.

Tools: R, Seurat, DoubletFinder, SingleR, celldex, Monocle3, clusterProfiler, org.Hs.eg.db
All analyses are performed in **R**, using well-established frameworks including **Seurat, DoubletFinder, SingleR, Monocle3, fgsea**, and **msigdbr**.

---

## 📘 Scientific Background

EMT is a fundamental cellular program that enables epithelial cells to acquire mesenchymal features, contributing to:

- tumor invasiveness  
- therapeutic resistance  
- cellular plasticity  
- metastatic potential  

Although EMT is known to play important roles in many cancers, its contribution to **meningioma heterogeneity and progression** remains poorly understood.

Single-cell RNA-seq provides the resolution needed to:

- Quantify EMT activation per cell  
- Identify EMT-high tumor subpopulations  
- Map EMT progression via pseudotime  
- Characterize biological pathways underlying EMT states  

This repository documents an end-to-end computational workflow designed to accomplish these goals.

---

## 📂 Dataset

- **GEO Accession:** [GSE183655](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183655)  
- **Tissue:** Human meningioma tumors  
- **Technology:** 10x Genomics Chromium  

Raw data are not included in this repository and must be obtained from GEO.

---

# 🔁 Complete Analysis Workflow (Summary)

### **1. Data Import & Preprocessing**
- Load raw 10X matrices  
- Create per-sample Seurat objects  
- Merge samples; extract patient, region, and barcode metadata  
- Compute QC metrics (UMIs, genes, mitochondrial %)  
- Apply filtering thresholds (nUMI, nGene, mito%)  

**Tools:** Seurat  

---

### **2. Doublet Detection & Removal**
DoubletFinder was applied **independently to each sample**:

- Parameter sweep to determine optimal **pK**  
- Estimated doublet rates per sample  
- Removal of predicted doublets  
- Reconstruction of a clean merged object  

**Tools:** DoubletFinder, Seurat  

---

### **3. Batch Correction & Sample Integration**
Anchor-based integration was performed to harmonize expression across patient samples:

- Variable feature identification  
- Anchor discovery  
- Integrated data construction  
- PCA and UMAP embedding  

**Tools:** Seurat (FindIntegrationAnchors, IntegrateData)  

---

### **4. Dimensionality Reduction & Clustering**
- Perform PCA  
- Construct shared nearest-neighbor graph  
- UMAP visualization  
- Clustering at multiple resolutions  
- Cluster–sample comparison  

**Tools:** Seurat  

---

### **5. EMT Scoring**
A curated set of EMT hallmark genes was used to compute EMT module scores:

**Genes:**  
`VIM, FN1, ZEB1, SNAI1, SNAI2, CDH2, TWIST1`

Cells were stratified into:

- **EMT_high**  
- **EMT_low**

This classification served as the foundation for downstream comparisons.

**Tools:** Seurat (AddModuleScore)

---

### **6. Differential Expression: EMT_high vs EMT_low**
Differentially expressed genes (DEGs) were computed to identify signatures associated with EMT activation.

Outputs include:

- log2 fold changes  
- adjusted p-values  
- top EMT markers  

**Tools:** Seurat (FindMarkers)

---

### **7. Gene Set Enrichment Analysis (GSEA)**

GSEA was performed on the ranked DEG list using GO terms from **org.Hs.eg.db**:

- Creation of ranked gene list  
- Running GSEA with clusterProfiler  
- Extraction of significantly enriched GO terms  
- Visualization of top enriched pathways  

**Enriched GO terms in EMT_high cells:**
- EMT-related biological processes  
- TGF-β–associated GO terms  
- Hypoxia-related GO terms  
- Inflammatory response pathways  
- IL6/JAK/STAT3-related terms  
- TNF/NF-κB-related signaling  

**Tools:** clusterProfiler, org.Hs.eg.db  

**Outputs:**
- `GSEA_results_GO.csv`  
- Top 10 enrichment barplot  
- EMT-related GO enrichment curve  
---

### **8. Cell Type Annotation**
SingleR was used for automated cell type identification using the **Human Primary Cell Atlas** reference dataset.

Diagnostic evaluations:

- score heatmaps  
- delta distributions  
- confusion matrices comparing cluster identity vs annotation  

**Tools:** SingleR, celldex, pheatmap  

---

### **9. Trajectory Inference (Monocle3)**
To understand EMT progression:

1. Convert Seurat object to Monocle CellDataSet  
2. Transfer Seurat UMAP embedding  
3. Define EMT_low cells as pseudotime root  
4. Learn principal graph  
5. Order cells along pseudotime  
6. Identify genes whose expression changes dynamically along EMT progression  

This reveals transcriptional transitions underlying EMT.

**Tools:** Monocle3, SeuratWrappers


---

# 🎯 Research Objectives

- Characterize transcriptional heterogeneity in meningioma  
- Identify EMT-associated tumor cell states  
- Quantify EMT activity at single-cell resolution  
- Predict functional pathways enriched in EMT-high populations  
- Reconstruct cellular trajectories representing EMT progression  
- Provide a reproducible computational pipeline for EMT research  

---

# 🧪 Planned Extensions

Although GSEA is completed, future analyses may include:

- KEGG & Reactome pathway enrichment  
- Gene Ontology (GO) biological process enrichment  
- Cell–cell communication using **CellChat** or **NicheNet**  
- Regulatory network inference (SCENIC)  

---

# ⚙️ Software Environment

- R ≥ 4.1  
- Seurat  
- DoubletFinder  
- SingleR + celldex  
- Monocle3  
- fgsea + msigdbr  
- tidyverse  
- pheatmap  
- ggsci  

---

# 🚀 Reproducibility

Each stage of the workflow is implemented as an independent script, enabling full end-to-end reproducibility. Starting from raw GEO data, the entire pipeline—from QC to DE, GSEA, annotation, and pseudotime—can be rerun using the scripts in `/scripts`.

---

# 👤 Author

**Sara Tolba**  
📧 sarahtolba842@gmail.com  
🔗 https://github.com/sarahtolba  

---

# 🙏 Acknowledgments

This work builds upon open-source tools including:

- **Seurat**  
- **DoubletFinder**  
- **SingleR** and **celldex**  
- **Monocle3**  
- **fgsea** and **msigdbr**  

Their contributions enable robust and reproducible single-cell research.

# EMT in Tumor Progression: Single-Cell RNA-seq Analysis Pipeline

## Project Overview

This project aims to investigate the role of **epithelial–mesenchymal transition (EMT)** in tumor progression using **single-cell RNA sequencing (scRNA-seq)** data.  
By analyzing cell-type-specific EMT activity across tumor and normal samples, we seek to identify molecular signatures and cell states associated with EMT, which may contribute to cancer metastasis and therapy resistance.

---

## Dataset Information

- **Source:** Publicly available scRNA-seq data  
- **Accession Number:** [GSE183655](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183655)  
- **Description:** Meningioma tumor samples with associated normal tissues, including multiple patient samples.

---

## Project Structure

/data/           # Raw sequencing data (not included; download from GEO accession)  
/seurat_objects/ # Serialized Seurat objects for different processing stages (excluded due to size)  
/results/        # Differential expression results, marker gene tables, metadata summaries  
/figures/        # Visualizations: UMAPs, violin plots, heatmaps, etc.  
/scripts/        # R scripts implementing each step of the pipeline  
.gitignore       # Specifies files and folders excluded from Git tracking  
README.md        # This file

---

## Analysis Pipeline

| Step                         | Description                                                            | Tools/Packages          |
|------------------------------|------------------------------------------------------------------------|------------------------|
| 1. Data Loading & QC          | Import raw data, filter cells based on quality metrics (mitochondrial %, gene counts) | Seurat                 |
| 2. Doublet Detection          | Identify and remove doublets to ensure data quality                    | DoubletFinder          |
| 3. Sample Integration         | Correct batch effects across samples using anchor-based integration   | Seurat                 |
| 4. EMT Scoring               | Quantify EMT gene expression signatures per cell                      | Seurat AddModuleScore  |
| 5. Differential Expression    | Identify genes differentially expressed between EMT-high vs EMT-low cells | Seurat FindMarkers     |
| 6. Cell Type Annotation       | Assign cell identities using reference datasets                        | SingleR                |
| 7. (Future) Trajectory Analysis | Infer EMT progression dynamics via pseudotime ordering               | Monocle3 / Slingshot   |
| 8. (Future) Cell-Cell Communication | Analyze intercellular signaling networks influencing EMT          | CellChat / NicheNet    |
| 9. (Future) Pathway Enrichment | Identify enriched biological pathways and processes                  | clusterProfiler / GSEA |

---

## Usage Instructions

1. **Download raw data** from [GEO Accession GSE183655](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183655) into the `/data` directory.

2. **Run analysis scripts** sequentially from the `/scripts` folder, or use provided R Markdown notebooks to reproduce the full workflow.

3. **Explore outputs** in `/results` and `/figures` for tables and visualizations relevant to EMT and cell type analyses.

---

## Software Dependencies

- R (≥ 4.1)  
- Seurat (≥ 4.0)  
- DoubletFinder  
- SingleR  
- tidyverse (ggplot2, dplyr, tidyr, etc.)  
- pheatmap  
- ggsci  
- clusterProfiler (for pathway enrichment, future)  
- Monocle3 or Slingshot (for trajectory inference, future)  
- CellChat or NicheNet (for cell-cell communication, future)

---

## Notes

- Large files such as raw sequencing data and Seurat objects are **not** included in this repository due to size limitations.  
- `.gitignore` is configured to exclude heavy files such as `.rds` Seurat objects and raw data.  
- Contact for access to processed data or collaboration inquiries.

---

## Contact Information

**Sarah Tolba**  
  
Email: [sarahtolba842@gmail.com](mailto:sarahtolba842@gmail.com)  
GitHub: [https://github.com/sarahtolba](https://github.com/sarahtolba)

---

## Acknowledgments

This project was inspired by methods and tools widely used in cutting-edge single-cell transcriptomics research.  
Special thanks to the developers of Seurat, SingleR, DoubletFinder, and other open-source tools used in this analysis.

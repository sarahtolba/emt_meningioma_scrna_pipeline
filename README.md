# 🧬 Single-Cell RNA-seq EMT Project

This project investigates the epithelial–mesenchymal transition (EMT) program in brain tumor samples using single-cell RNA sequencing (scRNA-seq). It focuses on identifying EMT-high subpopulations, characterizing their gene expression profiles, and mapping their developmental trajectories. Analyses were conducted using R and several state-of-the-art packages including Seurat, SingleR, DoubletFinder, and Monocle3.

---

## 📂 Dataset  
- **GSE183655**: scRNA-seq of meningioma tumor samples

---

## ✅ Completed Workflow

| Step                         | Description                                                                 | Tools/Packages          |
|------------------------------|-----------------------------------------------------------------------------|--------------------------|
| **1. Data Preprocessing**     | Imported raw data, filtered low-quality cells                              | Seurat                   |
| **2. Doublet Removal**        | Identified and removed artificial doublets                                  | DoubletFinder            |
| **3. Sample Integration**     | Integrated multiple samples and corrected for batch effects                 | Seurat (SCTransform)     |
| **4. EMT Scoring**            | Scored cells using curated EMT gene signature                              | Seurat (AddModuleScore)  |
| **5. Clustering & Annotation**| Clustered cells and assigned cell types based on reference data             | Seurat + SingleR         |
| **6. Differential Expression**| Identified DEGs between EMT-high vs EMT-low populations                     | Seurat (FindMarkers)     |
| **7. Trajectory Inference**   | Reconstructed pseudotime trajectories to model EMT dynamics                 | Monocle3                 |

---

## 📁 Project Structure

/data/ # Raw sequencing data (not included; download from GEO accession)
/seurat_objects/ # Serialized Seurat objects (excluded due to size)
/results/ # DEGs, marker tables, summaries
/figures/ # UMAPs, heatmaps, violin plots, etc.
/scripts/ # R scripts for each step
.gitignore # Files/folders excluded from Git tracking
README.md # This file

---

## 🧪 🔜 Planned Analysis

| Step                            | Description                                                  | Tools (to be used)       |
|---------------------------------|--------------------------------------------------------------|---------------------------|
| **1. Gene Set Enrichment**      | GO/KEGG/Reactome enrichment for DEGs and marker genes        | clusterProfiler / enrichR |
| **2. Cell–Cell Communication**  | Inference of signaling interactions affecting EMT            | CellChat / NicheNet       |

---

## 📌 Goals

- Identify EMT-driven cell subpopulations  
- Characterize EMT-related gene signatures  
- Understand EMT progression over pseudotime  
- Reveal molecular pathways and signaling mechanisms involved in EMT  

---

## ⚙️ Software Dependencies

- R (≥ 4.1)  
- Seurat (≥ 4.0)  
- DoubletFinder  
- SingleR  
- tidyverse  
- pheatmap  
- ggsci  
- clusterProfiler (for future enrichment)  
- Monocle3  
- CellChat / NicheNet (planned)

---

## 🚀 Usage Instructions

1. **Download raw data** from [GEO: GSE183655](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183655) into `/data`
2. **Run analysis scripts** in order from the `/scripts` directory
3. **View results** in `/results` and `/figures`

---

## 📝 Notes

- Raw data and `.rds` Seurat objects are excluded due to size.
- `.gitignore` excludes large files and intermediate outputs.
- Contact for access to processed objects or collaborative questions.

---

## 👤 Contact

**Sarah Tolba**  
📧 [sarahtolba842@gmail.com](mailto:sarahtolba842@gmail.com)  
🌐 [github.com/sarahtolba](https://github.com/sarahtolba)

---

## 🙏 Acknowledgments

Thanks to the developers of Seurat, SingleR, Monocle3, and other open-source packages used in this pipeline.


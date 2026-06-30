# scRNA-seq Analysis of Human Bone Marrow Immune Cells

End-to-end single-cell RNA-seq analysis pipeline in R (Seurat) covering quality control through cell-cell communication, applied to four human bone marrow samples (BMMC and CD34+ enriched progenitors) from Granja et al. (2019).

## Overview

This project processes ~21,000 cells across four donors to build an annotated immune cell atlas, characterize differentiation trajectories within the hematopoietic lineage, and map cell-cell signaling networks. It was completed as part of the Single Cell Bioinformatics course at Saarland University (ICBB, Müller Lab).

**Samples:** BMMC_D1T1, BMMC_D1T2 (bone marrow mononuclear cells) and CD34_D2T1, CD34_D3T1 (CD34+ enriched progenitors)

## Pipeline

1. **QC & Filtering** — per-sample filtering on `nFeature_RNA` and `nCount_RNA` thresholds, mitochondrial content checks
2. **Doublet removal** — DoubletFinder with per-sample pK optimization
3. **Normalization & feature selection** — LogNormalize, `FindVariableFeatures`
4. **Batch correction** — `FindIntegrationAnchors` / `IntegrateData` (compared against uncorrected merge)
5. **Dimensionality reduction** — PCA (20 PCs via elbow plot) → UMAP
6. **Clustering** — Louvain clustering, 19 clusters
7. **Cell type annotation**
   - Automatic: SingleR + Human Primary Cell Atlas reference
   - Manual: marker-based (CD3D, CD14, CD19, GATA2, etc. — see Table 2 in `docs/`)
8. **Differential expression** — B cells vs T cells, T cells vs Monocytes, BMMC vs CD34, volcano plots, dot plots
9. **Pathway enrichment** — EnrichR across ENCODE, GO, KEGG, MSigDB Hallmark, and Human Gene Atlas
10. **Trajectory analysis** — Monocle3, pseudotime along the HSC → LMPP → CLP differentiation path, manual vs automatic root node selection
11. **Cell-cell communication** — CellChat, comparing interaction number/strength between BMMC and CD34 samples, with a focused look at the MIF signaling pathway

## Key results

- Clear separation of batch effects in uncorrected UMAP, resolved after `IntegrateData` correction (see `figures/UMAP_batch_correction.png`)
- 18 manually annotated cell populations spanning stem/progenitor, B-lineage, T-lineage, NK, and myeloid compartments (`figures/UMAP_manual_annotation.png`)
- Differentiation trajectory recovered along HSC → LMPP → CLP with biologically consistent pseudotime ordering (`figures/Trajectory_pseudotime.png`)
- MIF signaling identified as a shared, active pathway across both BMMC and CD34 cell-cell communication networks (`figures/CellChat_MIF_pathway.png`)
- Top enriched pathways (BMMC vs CD34 comparison) included MYC target genes, mRNA splicing/spliceosome machinery, and neurodegeneration-associated KEGG pathways (mitochondrial/proteostasis genes also enriched in hematopoietic progenitors)

## Repository structure

```
scripts/    Assignment1_scRNA.R   — full annotated analysis script
figures/    Key output plots (UMAP, enrichment, trajectory, CellChat)
docs/       Full_Report.pdf       — complete write-up with all 40 figures, tables, and interpretation
```

## Tools & packages

R, Seurat, DoubletFinder, SingleR, celldex, EnrichR, Monocle3, CellChat, SeuratWrappers, tidyverse

## Reference

Granja, J. M. et al. (2019). Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia. *Nature Biotechnology*.

---


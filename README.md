# UPR-Net

Single-cell causal inference and clinical radiomics of endoplasmic-reticulum (ER) stress rewiring in drug-resistant tuberculosis.

## Research status

This repository contains the analysis code for an **ongoing study**. The manuscript is under preparation and **has not yet been published**. The scripts are provided to support reproducibility of the computational analyses.

## Overview

Tuberculosis drug resistance is shaped in part by host macrophage programs, yet the intracellular signaling that determines whether a macrophage clears or tolerates mycobacteria remains poorly resolved. This project combines three complementary layers of evidence:

- **Single-cell causal discovery** — pathway-constrained neural differential equations and structure learning over the unfolded protein response (UPR)
- **Independent validation** — counterfactual (biolord) and optimal-transport (moscot) temporal analyses
- **Clinical CT radiomics** — whole-lung radiomic features linked to single-cell molecular indices

## Scripts

### Causal discovery

| Script | Input | Output | Environment |
|--------|-------|--------|-------------|
| `p3_causal_compare2.R` | Seurat object (`p3_seurat.rds`) | PC and GES causal networks, edge lists | R ≥ 4.6, `pcalg`, `igraph`, `Seurat` |
| `p3_bootstrap.R` | Seurat object | 100× bootstrap edge-stability matrix | R ≥ 4.6, `pcalg` |
| `p3_harmony.R` / `p3_harmony2.R` | Seurat object | Harmony batch-corrected embeddings | R ≥ 4.6, `harmony`, `Seurat` |
| `fig1_single.R` | `p3_mac_clus_full.csv`, `p3_dotplot_data.csv` | Single-cell landscape figures | R ≥ 4.6, `ggplot2`, `viridis`, `ggsci` |
| `fig1_moscot.R` | moscot transition costs | Temporal transition figure | R ≥ 4.6, `ggplot2` |
| `fig2_single.R` | Frozen causal-effect values | Causal network, dumbbell, heatmap, lollipop | R ≥ 4.6, `ggraph`, `ggalt`, `ComplexHeatmap` |
| `fig3_single.R` | `causal_comparison.rds`, `bootstrap_edge_prob.rds` | Method-consensus and stability heatmaps | R ≥ 4.6, `pheatmap` |

### Network rewiring

| Script | Input | Output | Environment |
|--------|-------|--------|-------------|
| `p3_part3_clean.R` | `counts_T.mtx`, metadata | Pseudobulk DEGs, fGSEA enrichment, volcano/heatmap/barplot | R ≥ 4.6, `Seurat`, `fgsea`, `msigdbr`, `ggrepel` |

### Clinical radiomics

| Script | Input | Output | Environment |
|--------|-------|--------|-------------|
| `radiomics_extract.py` | Chest CT DICOM | Whole-lung PyRadiomics features (`radiomics_features.csv`) | Python 3.11, `PyRadiomics`, `SimpleITK`, `pydicom` |
| `radiomics_plot.R` | `radiomics_features.csv`, `molecular_scores_7.csv` | Radiomic–molecular correlation, scatter, boxplot, PCA | R ≥ 4.6, `ggplot2`, `pheatmap` |
| `render_ct_slices.py` | Chest CT DICOM | Lung-window CT slice images (PNG) | Python 3.11, `pydicom`, `matplotlib` |

### Dynamics

| Script | Input | Output | Environment |
|--------|-------|--------|-------------|
| `export_scdiffeq_velocity.py` | `tb_15samples.h5ad` | scDiffEq velocity vectors | Python 3.11, `scdiffeq`, `scanpy`, `torch` (CPU) |
| `install_scdiffeq.sh` | — | scDiffEq installation (CPU torch) | — |

## Data and acknowledgments

Clinical samples and data were provided with the support of:

- Xinjiang Second Medical College (新疆第二医学院)
- Shihezi University School of Medicine (石河子大学医学院)
- Xinjiang Medical University (新疆医科大学)
- The Sixth People's Hospital of Xinjiang Uygur Autonomous Region (新疆维吾尔自治区第六人民医院)

Public single-cell and bulk transcriptomic data were obtained from the Gene Expression Omnibus (GEO).

## Citation

The manuscript is not yet published. Once available, please cite the corresponding publication. For now, you may cite this repository.

## License

This project is licensed under the MIT License.

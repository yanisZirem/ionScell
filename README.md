<p align="center">
  <img src="logo.png" width="180" alt="IonScell logo"/>
</p>

<h1 align="center">IonScell</h1>
<p align="center"><em>Single-Cell Mass Spectrometry Imaging — End-to-End Analysis Pipeline</em></p>

<p align="center">
  <img src="https://img.shields.io/badge/Python-3.9%2B-blue?style=flat-square"/>
  <img src="https://img.shields.io/badge/License-MIT-green?style=flat-square"/>
</p>

---

## Overview

**IonScell** is a Python pipeline for **Single-Cell Mass Spectrometry Imaging (SCMSI)**. It processes raw imzML files acquired on MALDI, DESI, or other MSI platforms and performs the full analytical workflow, from data loading and cell segmentation, through soft clustering and quality control, to lipid/metabolite annotation and export.

IonScell is available in three complementary formats:

| Format | Best for |
|--------|----------|
| **Jupyter Notebook** (`ionScell_notebook.ipynb`) | Step-by-step exploratory analysis, publication figures |
| **Python class** (`ionscell_pipeline.py`) | Integration into existing pipelines, scripting |
| **Desktop application** (`ionScell.exe`) | Non-programmers, rapid interactive analysis |

---

## Key Features

- **Flexible m/z binning** — Fixed (0.1 Da) or adaptive (instrument-matched ppm resolution)
- **Adaptive TIC thresholding** — 8 algorithms (GMM, Otsu, Triangle, Li, Yen, …) with automatic selection
- **Cell segmentation** — Watershed and DFS algorithms with size filtering in pixels or µm
- **Soft clustering** — Gaussian Mixture Models (GMM) and Fuzzy C-Means (FCM) with probabilistic membership
- **Quality control** — SNR, sparsity, dynamic range, spectral entropy, number of detected peaks
- **Clone quality assessment** — Confidence scores, within-clone CV, between-clone separability
- **Differential analysis** — Kruskal-Wallis or ANOVA with BH FDR correction, volcano plots, Z-score heatmaps
- **Ion spatial mapping** — Single-ion and multi-ion overlay with cell contours
- **Lipid annotation** — Embedded LIPID MAPS / HMDB database, no internet required, adduct-aware
- **Export** — CSV (metadata + spectra), imzML

---

## Installation

### 1. Clone or download the repository

```bash
git clone https://github.com/your-org/ionscell.git
cd ionscell
```

### 2. Create a dedicated environment (recommended)

```bash
conda create -n ionscell python=3.10 -y
conda activate ionscell
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

**`requirements.txt`**
```
numpy>=1.24
pandas>=2.0
scipy>=1.11
scikit-image>=0.21
scikit-learn>=1.3
plotly>=5.18
matplotlib>=3.7
pyimzml>=1.5
tqdm>=4.66
umap-learn>=0.5
```

Optional (for Fuzzy C-Means):
```bash
pip install scikit-fuzzy
```

---

## Quick Start

### Option A — Jupyter Notebook

```bash
jupyter notebook ionScell_notebook.ipynb
```

Follow the numbered cells (Steps 1–12) to process your dataset interactively.

### Option B — Python script

```python
from ionscell_pipeline import SCMSIPipeline, LipidAnnotator

# 1. Load data
pipe = SCMSIPipeline()
pipe.load_imzml("path/to/your/data.imzML")
pipe.build_datacube(mz_min=350, mz_max=1100, mz_bin=0.1, normalization="TIC")

# 2. Visualise TIC
pipe.show_TIC_plotly()

# 3. Adaptive thresholding + segmentation
thresh, _ = pipe.auto_TIC_threshold(method='auto')
pipe.watershed_or_dfs_from_mask(TIC_threshold=thresh, method='watershed',
                                 min_cell_size=4, pixel_size_um=10)

# 4. Extract spectra + QC
pipe.extract_cell_spectra(agg='mean')
pipe.compute_spectral_quality_metrics()
pipe.filter_low_quality_cells(snr_threshold=2.0)

# 5. Soft clustering
pipe.compute_soft_clustering(n_clusters=None, max_clusters=8)
pipe.assess_clone_quality()

# 6. Visualise
pipe.show_umap_with_quality_overlay()
pipe.overlay_clusters_on_image_plotly()
pipe.plot_clone_spectra_with_uncertainty()

# 7. Differential analysis
top_markers, _ = pipe.run_differential_analysis(method='kruskal')

# 8. Annotation
ann_df, enrichment = pipe.annotate_ions(mode='neg', ppm_tolerance=10.0)

# 9. Export
pipe.export_cells_to_csv("results/cells.csv")
```

### Option C — Desktop GUI

Double-click `ionScell.exe`. See the **Interface User Manual** for full instructions.

---

## Pipeline Architecture

```
imzML file
    │
    ▼
load_imzml()
    │
    ▼
make_mz_axis() / build_adaptive_mz_axis()
    │
    ▼
build_datacube()          ← normalization: TIC / RMS / MEDIAN
    │
    ▼
auto_TIC_threshold()      ← 8 algorithms, automatic selection
    │
    ▼
watershed_or_dfs_from_mask()   ← Watershed or DFS
    │                            ← size filtering (px or µm)
    ├── remove_cells()           ← manual artefact removal
    │
    ▼
extract_cell_spectra()    ← mean / median / sum per cell
    │
    ▼
compute_spectral_quality_metrics()
filter_low_quality_cells()
    │
    ▼
compute_soft_clustering()      ← GMM (auto k via BIC/silhouette)
  or compute_fuzzy_cmeans()    ← FCM
    │
    ▼
assess_clone_quality()
    │
    ├── show_umap_with_quality_overlay()
    ├── overlay_clusters_on_image_plotly()
    ├── plot_clone_spectra_with_uncertainty()
    ├── show_mz_distribution_in_cells_plotly()
    ├── plot_multi_ion_overlay()
    ├── plot_distribution_by_clone()
    │
    ▼
run_differential_analysis()   ← Kruskal-Wallis / ANOVA + BH FDR
    │
    ▼
annotate_ions()               ← LipidAnnotator (embedded DB)
    │
    ▼
export_cells_to_csv()
export_cells_to_imzML()
```

---

## API Reference

### `SCMSIPipeline`

#### Data loading & cube building

| Method | Description |
|--------|-------------|
| `load_imzml(path)` | Load imzML file |
| `make_mz_axis(mz_min, mz_max, mz_bin)` | Fixed-grid m/z axis |
| `build_adaptive_mz_axis(...)` | Instrument-matched adaptive binning |
| `build_datacube(mz_min, mz_max, normalization)` | Build 3D data cube |

#### Thresholding & segmentation

| Method | Description |
|--------|-------------|
| `auto_TIC_threshold(method)` | Adaptive TIC threshold (8 methods + auto) |
| `watershed_or_dfs_from_mask(...)` | Cell segmentation |
| `remove_cells(remove_ids)` | Remove artefact cells |

#### Spectral extraction & QC

| Method | Description |
|--------|-------------|
| `extract_cell_spectra(agg)` | One spectrum per cell |
| `compute_spectral_quality_metrics()` | SNR, entropy, sparsity, … |
| `filter_low_quality_cells(...)` | Flag/remove low-quality cells |

#### Clustering

| Method | Description |
|--------|-------------|
| `compute_soft_clustering(n_clusters)` | GMM soft clustering |
| `compute_fuzzy_cmeans(max_clusters)` | Fuzzy C-Means |
| `assess_clone_quality()` | Quality report per clone |

#### Clone management

| Method | Description |
|--------|-------------|
| `remove_clone(clone_id)` | Remove entire clone |
| `combine_clones(clone_ids)` | Merge clones |
| `remove_clusters(clusters_to_remove)` | Remove by cluster ID |

#### Visualisation

| Method | Description |
|--------|-------------|
| `show_TIC_plotly()` | TIC heatmap |
| `show_umap_with_quality_overlay()` | UMAP coloured by clone + confidence |
| `overlay_clusters_on_image_plotly()` | Spatial clone map |
| `plot_cell_contours_by_clone()` | Contours on white/TIC background |
| `plot_clone_spectra_with_uncertainty()` | Mean spectra ± SD |
| `plot_global_spectrum()` | Global mean spectrum |
| `show_mz_distribution_in_cells_plotly(mz_value)` | Single-ion map |
| `plot_multi_ion_overlay(mz_list)` | Multi-ion RGB overlay |
| `plot_distribution_by_clone(mz_value)` | Violin/box plots |

#### Analysis & export

| Method | Description |
|--------|-------------|
| `run_differential_analysis(method)` | Statistical DA + volcano + heatmap |
| `annotate_ions(mode, ppm_tolerance)` | Lipid/metabolite annotation |
| `export_cells_to_csv(path)` | Export to CSV |
| `export_cells_to_imzML(path)` | Export to imzML |

---

### `LipidAnnotator`

```python
ann = LipidAnnotator(mode='pos', ppm_tolerance=5.0)
hits = ann.annotate_peak(800.56)          # single m/z
df   = ann.annotate_peak_list(mz_array)   # array of m/z
results, enrichment = ann.annotate_clones(pipeline)
df   = ann.search("ceramide")             # free-text search
```

**Supported adducts** (positive): `[M+H]+`, `[M+Na]+`, `[M+K]+`, `[M+NH4]+`, `[M+Li]+`, `[M+H-H2O]+`, `[M+2H]2+`  
**Supported adducts** (negative): `[M-H]-`, `[M+Cl]-`, `[M+FA-H]-`, `[M+Ac-H]-`, `[M-H-H2O]-`, `[M+Br]-`, `[M-2H]2-`

---

## Repository Structure

```
ionscell/
├── ionscell_pipeline.py        ← SCMSIPipeline + LipidAnnotator class
├── ionScell_notebook.ipynb     ← Clean Jupyter notebook (imports class)
├── requirements.txt
├── README.md
├── logo.png
├── docs/
│   ├── user_manual_jupyter.md  ← Jupyter + Google Colab user manual
│   └── user_manual_interface.md← Desktop interface user manual
└── _internal/                  ← (bundled app dependencies)
```

---

## Supported Instruments & Formats

| Instrument type | Format | Notes |
|----------------|--------|-------|
| MALDI-TOF/TOF   | imzML  | Continuous or processed mode |
| MALDI-Orbitrap  | imzML  | Use adaptive binning for best annotation |
| DESI-MS         | imzML  | |
| SIMS            | imzML  | |

---

## Citation

If you use IonScell in your research, please cite:

> **IonScell: A Python pipeline for single-cell mass spectrometry imaging**  
> *[Authors, Journal, Year — in preparation]*

---

## License

MIT License — see `LICENSE` for details.

---

## Contact

Issues and feature requests: [GitHub Issues](https://github.com/your-org/ionscell/issues)

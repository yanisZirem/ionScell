# IonScell — Jupyter Notebook & Google Colab User Manual

<p align="center"><em>Version 1.0 · Single-Cell Mass Spectrometry Imaging</em></p>

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Installation & Setup](#2-installation--setup)
3. [Google Colab Quick-Start](#3-google-colab-quick-start)
4. [Notebook Structure Overview](#4-notebook-structure-overview)
5. [Step-by-Step Walkthrough](#5-step-by-step-walkthrough)
   - 5.1 Loading your imzML file
   - 5.2 Building the data cube
   - 5.3 Visualising the TIC
   - 5.4 Adaptive TIC thresholding
   - 5.5 Cell segmentation
   - 5.6 Manual artefact removal
   - 5.7 Spectral extraction
   - 5.8 Quality control
   - 5.9 Soft clustering (GMM / FCM)
   - 5.10 Clone quality assessment
   - 5.11 Visualisations
   - 5.12 Differential analysis
   - 5.13 Lipid annotation
   - 5.14 Export
6. [Key Parameters Reference](#6-key-parameters-reference)
7. [Troubleshooting](#7-troubleshooting)
8. [FAQ](#8-faq)

---

## 1. Introduction

The IonScell Jupyter Notebook provides a **step-by-step, interactive** analysis of single-cell MSI datasets. All heavy-lifting is performed by the `SCMSIPipeline` and `LipidAnnotator` classes imported from `ionscell_pipeline.py` — the notebook itself contains **no class definitions**, only calls and visualisations.

This separation means you can:
- Reuse the class in your own scripts without the notebook
- Update the class independently from the notebook
- Run the notebook on Google Colab without a local Python installation

---

## 2. Installation & Setup

### Local Jupyter

```bash
# 1. Create environment
conda create -n ionscell python=3.10 -y
conda activate ionscell

# 2. Install dependencies
pip install numpy pandas scipy scikit-image scikit-learn plotly \
            matplotlib pyimzml tqdm umap-learn jupyterlab

# 3. Optional: Fuzzy C-Means
pip install scikit-fuzzy

# 4. Launch notebook
jupyter lab ionScell_notebook.ipynb
```

Place `ionscell_pipeline.py` in the **same folder** as the notebook, or add its path:

```python
import sys
sys.path.insert(0, "/path/to/ionscell_pipeline.py")
```

---

## 3. Google Colab Quick-Start

Click the badge below to open the notebook directly in Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/your-org/ionscell/blob/main/ionScell_notebook.ipynb)

### Setup cell (run once at the top of your Colab session)

```python
# ── Install dependencies ───────────────────────────────────────────────
!pip install pyimzml umap-learn scikit-fuzzy tqdm plotly -q

# ── Upload pipeline class ──────────────────────────────────────────────
from google.colab import files
uploaded = files.upload()          # upload ionscell_pipeline.py

# ── Upload your imzML + ibd files ─────────────────────────────────────
data = files.upload()              # upload .imzML AND .ibd together
# Note: both .imzML and .ibd must be in the same directory

# ── (Alternative) Mount Google Drive ──────────────────────────────────
from google.colab import drive
drive.mount('/content/drive')
imzml_path = '/content/drive/MyDrive/your_data/sample.imzML'
```

> **Tip:** For large files (> 500 MB), mounting Google Drive is strongly recommended over direct upload.

### GPU acceleration (optional)

For UMAP on large datasets (> 5 000 cells), switch to a GPU runtime:  
`Runtime → Change runtime type → T4 GPU`

Then install the GPU-accelerated UMAP:
```python
!pip install cuml-cu11  # RAPIDS cuML (CUDA 11)
```

---

## 4. Notebook Structure Overview

The notebook is organised into **12 numbered sections** that mirror the pipeline architecture:

| Section | Title | Key output |
|---------|-------|-----------|
| 0 | Imports & path setup | — |
| 1 | Load imzML | `pipe.parser`, `pipe.coords` |
| 2 | Build data cube | `pipe.data_cube`, `pipe.TIC` |
| 3 | TIC visualisation | Interactive heatmap |
| 4 | Adaptive thresholding | `thresh` value |
| 5 | Cell segmentation | `pipe.cell_labels` |
| 6 | Spectral extraction | `pipe.last_cell_spectra` |
| 7 | Quality control | QC metrics, flagged cell IDs |
| 8 | Soft clustering | Cluster assignments, UMAP |
| 9 | Clone assessment | Quality report |
| 10 | Differential analysis | `top_markers` DataFrame |
| 11 | Lipid annotation | `ann_df`, enrichment heatmap |
| 12 | Export | CSV / imzML files |

---

## 5. Step-by-Step Walkthrough

### 5.1 Loading your imzML file

```python
from ionscell_pipeline import SCMSIPipeline, LipidAnnotator

pipe = SCMSIPipeline(verbose=True)
pipe.load_imzml("path/to/sample.imzML")
```

**Expected output:**
```
imzML loaded: 12 450 spectral pixels; coord min = (1, 1)
```

> Both the `.imzML` and `.ibd` files must be in the **same directory** and share the same base name.

---

### 5.2 Building the data cube

**Fixed binning (recommended for most workflows):**

```python
pipe.build_datacube(
    mz_min=350,
    mz_max=1100,
    mz_bin=0.1,           # 0.1 Da bins
    normalization="TIC",  # or "RMS", "MEDIAN", None
)
```

**Adaptive binning (recommended for annotation workflows):**

```python
pipe.build_adaptive_mz_axis(
    mz_min=350, mz_max=1100,
    target_resolving_power=20000   # or ppm_bin=50.0
)
pipe.build_datacube()   # uses the axis just created
```

**Memory-efficient mode (large datasets > 8 GB):**

```python
pipe.build_datacube(
    mz_min=350, mz_max=1100,
    memmap_path="/tmp/ionscell_cube.dat"
)
```

---

### 5.3 Visualising the TIC

```python
pipe.show_TIC_plotly(
    sigma=0.5,       # Gaussian smoothing
    zoom_factor=1.0,
    cmap='magma'
)
```

The interactive Plotly figure lets you zoom, pan, and hover to inspect pixel intensities. Use this to:
- Identify tissue regions
- Spot artefacts (very bright/dark spots)
- Check spatial coverage

---

### 5.4 Adaptive TIC thresholding

```python
thresh, diagnostics = pipe.auto_TIC_threshold(
    method='auto',          # recommended
    show_diagnostics=True,  # shows candidate comparison plot
    min_cell_fraction=0.05, # minimum expected cell coverage
    max_cell_fraction=0.60  # maximum expected cell coverage
)
print(f"Selected threshold: {thresh:.1f}")
```

Available methods: `'auto'`, `'gmm'`, `'otsu'`, `'triangle'`, `'li'`, `'yen'`, `'isodata'`, `'percentile'`

The diagnostic plot shows all candidate thresholds with their scores. The selected method is highlighted.

> **Tip:** If the automatic selection misses many cells, try `method='otsu'` or lower `min_cell_fraction` to `0.02`.

---

### 5.5 Cell segmentation

```python
pipe.watershed_or_dfs_from_mask(
    TIC_threshold=thresh,
    method='watershed',      # or 'dfs'
    sigma=0.5,               # pre-segmentation smoothing
    footprint=(5, 5),        # local maxima neighbourhood
    min_cell_size=4,         # minimum 4 pixels
    max_cell_size=None,      # no upper limit
    pixel_size_um=10,        # physical pixel size
    size_unit='px',          # 'px' or 'um'
)
```

**Parameter tuning guide:**

| Parameter | Effect | Increase if… | Decrease if… |
|-----------|--------|-------------|-------------|
| `sigma` | Pre-smoothing | Too many fragments | Cells merged |
| `footprint` | Seed spacing | Over-segmentation | Under-segmentation |
| `min_cell_size` | Noise removal | Too many 1-pixel fragments | Real small cells removed |

**DFS segmentation** (alternative to watershed):
```python
pipe.watershed_or_dfs_from_mask(
    TIC_threshold=thresh,
    method='dfs',
    intensity_tolerance=0.15  # max intensity difference within a cell
)
```

---

### 5.6 Manual artefact removal

After visualising the segmentation, you may identify cells with suspicious IDs (artefacts, tissue folds, debris). Remove them:

```python
# Remove specific cells
pipe.remove_cells(
    remove_ids=[12, 45, 67],
    show_plotly=True
)
```

Hover over the segmentation figure to identify cell IDs before removal.

---

### 5.7 Spectral extraction

```python
df = pipe.extract_cell_spectra(agg='mean')  # 'mean', 'sum', or 'median'
print(df.shape)   # (n_cells, n_metadata + n_mz_bins)
df.head()
```

The returned DataFrame has:
- Metadata columns: `cell_id`, `area_px`, `area_um`, `centroid_row`, `centroid_col`
- One column per m/z bin (column name = m/z value as string)

---

### 5.8 Quality control

**Compute QC metrics:**

```python
qc = pipe.compute_spectral_quality_metrics()
# Adds columns: qc_SNR, qc_sparsity, qc_dynamic_range,
#               qc_total_intensity, qc_spectral_entropy, qc_n_peaks
```

**Filter low-quality cells:**

```python
flagged = pipe.filter_low_quality_cells(
    snr_threshold=2.0,
    intensity_percentile=10,     # remove bottom 10% by intensity
    sparsity_min=0.05,
    auto_remove=False,           # set True to remove automatically
    show_plot=True
)
print(f"Flagged cells: {flagged}")
```

The QC plot compares "Good" vs "Flagged" cells across all 6 metrics.

---

### 5.9 Soft clustering

**GMM (recommended):**

```python
memberships, labels, confidence = pipe.compute_soft_clustering(
    n_clusters=None,      # None = auto via BIC + silhouette
    max_clusters=8,       # upper bound for auto search
    covariance_type='tied',
    use_umap_space=True,
    n_neighbors=15,
    random_state=42
)
```

**Fuzzy C-Means (alternative):**

```python
memberships, labels, confidence = pipe.compute_fuzzy_cmeans(
    max_clusters=8,
    m=2.0,                # fuzziness exponent
    use_umap_space=True
)
```

After clustering, `pipe.last_cell_spectra` gains columns:
- `Class` — hard cluster assignment
- `soft_confidence` — max membership probability
- `prob_cluster_0`, `prob_cluster_1`, … — per-cluster probabilities
- `umap_0`, `umap_1` — 2D UMAP coordinates

---

### 5.10 Clone quality assessment

```python
quality = pipe.assess_clone_quality(
    min_confidence=0.6,
    min_cells_per_clone=5,
    show_report=True
)
```

**Status levels:**

| Status | Meaning |
|--------|---------|
| `GOOD` | High confidence, homogeneous, well-separated |
| `QUESTIONABLE` | 1–2 quality flags — review manually |
| `ARTIFACT` | ≥ 3 flags — likely noise or debris |
| `LOW_SIZE` | Too few cells — may merge with neighbour |

**Merge questionable clones:**

```python
pipe.combine_clones([2, 3])   # merge clones 2 and 3 into one
```

**Remove an artefact clone:**

```python
pipe.remove_clone(clone_id=4)
```

---

### 5.11 Visualisations

**UMAP coloured by clone (opacity = confidence):**
```python
pipe.show_umap_with_quality_overlay()
```

**Spatial TIC map with clone contours:**
```python
pipe.overlay_clusters_on_image_plotly(
    opacity_by_confidence=True,
    contour_width=1.0
)
```

**Cell contours on white background:**
```python
pipe.plot_cell_contours_by_clone(show_tic=False)
```

**Mean spectra per clone ± SD:**
```python
pipe.plot_clone_spectra_with_uncertainty(
    agg='mean',
    show_std=True,
    normalize=True,
    mode='overlay'    # or 'individual'
)
```

**Single-ion spatial distribution:**
```python
pipe.show_mz_distribution_in_cells_plotly(
    mz_value=760.59,   # e.g. PC(34:1) [M+H]+
    tolerance=0.1,
    smoothing_sigma=1.0
)
```

**Multi-ion overlay (up to 6 ions):**
```python
pipe.plot_multi_ion_overlay(
    mz_list=[760.59, 782.57, 810.60],
    colors=['#ff4d6d', '#00f5d4', '#fee440'],
    normalize_each=True
)
```

**Violin/box plots by clone:**
```python
pipe.plot_distribution_by_clone(
    mz_value=760.59,
    tolerance=0.1,
    plot_type='violin',   # or 'box', 'both'
    points='scatter'
)
```

---

### 5.12 Differential analysis

```python
top_markers, full_results = pipe.run_differential_analysis(
    top_n=30,
    pval_threshold=0.05,
    fc_threshold=1.5,
    method='kruskal',       # or 'anova'
    show_volcano=True,
    show_heatmap=True
)
print(top_markers[['mz','pval_adj','fold_change','best_clone']].head(10))
```

The method:
1. Runs Kruskal-Wallis (or ANOVA) across all m/z bins
2. Applies Benjamini-Hochberg FDR correction
3. Computes fold-change relative to global mean
4. Returns a ranked DataFrame of significant ions
5. Displays a volcano plot and Z-score heatmap

---

### 5.13 Lipid annotation

```python
annotator = LipidAnnotator(mode='neg', ppm_tolerance=10.0)

# Annotate top differential ions
ann_df, enrichment = annotator.annotate_clones(
    pipe,
    top_markers_df=top_markers,
    show_table=True,
    show_enrichment=True
)

# Or annotate a custom list
hits = annotator.annotate_peak_list(
    mz_array=[760.59, 782.57, 810.60, 885.55]
)
print(hits[['name','adduct','ppm_error','subclass','pathways']])

# Free-text search
annotator.search("ceramide")
```

**Mode selection:**

| Ionisation | mode |
|-----------|------|
| Positive | `'pos'` |
| Negative | `'neg'` |
| Both / unknown | `'both'` |

---

### 5.14 Export

**CSV export (full data):**
```python
pipe.export_cells_to_csv(
    "results/ionscell_cells.csv",
    include_spectra=True,
    round_intensities=4
)
```

**CSV export (metadata only, no spectra):**
```python
pipe.export_cells_to_csv(
    "results/ionscell_metadata.csv",
    include_spectra=False
)
```

**imzML export:**
```python
pipe.export_cells_to_imzML("results/ionscell_cells.imzML")
```

---

## 6. Key Parameters Reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mz_min` | — | Lower m/z bound (Da) |
| `mz_max` | — | Upper m/z bound (Da) |
| `mz_bin` | `0.1` | Fixed bin width (Da) |
| `normalization` | `None` | `'TIC'`, `'RMS'`, `'MEDIAN'`, or `None` |
| `TIC_threshold` | auto | Intensity threshold for binary mask |
| `sigma` | `0.5` | Gaussian smoothing σ before segmentation |
| `footprint` | `(3,3)` | Local maxima footprint |
| `min_cell_size` | `2` | Minimum cell area (px or µm depending on `size_unit`) |
| `pixel_size_um` | `10` | Physical pixel size in µm |
| `n_clusters` | `None` | Number of clusters (None = auto) |
| `max_clusters` | `10` | Upper bound for auto cluster search |
| `covariance_type` | `'tied'` | GMM covariance structure |
| `ppm_tolerance` | `5.0` | Mass accuracy for annotation |
| `pval_threshold` | `0.05` | FDR cutoff for differential analysis |
| `fc_threshold` | `1.5` | Fold-change cutoff |

---

## 7. Troubleshooting

**`RuntimeError: ImzML not loaded`**  
→ Run `pipe.load_imzml()` before any other step.

**`No cells detected after segmentation`**  
→ Lower `TIC_threshold` manually, or use `method='percentile'` with `percentile=50`.

**`UMAP failed: ...`**  
→ `umap-learn` is not installed. Run `pip install umap-learn`. The pipeline will fall back to PCA automatically.

**`GMM failed for all k values`**  
→ Try increasing `reg_covar` to `1e-3`, or switch to `compute_fuzzy_cmeans()`.

**Very slow datacube build**  
→ Use `mem_efficient=True` or `memmap_path` for datasets > 4 GB.

**Plotly figures not displaying in JupyterLab**  
```bash
pip install jupyterlab-plotly
jupyter labextension install jupyterlab-plotly
```

**Colab: imzML file not found**  
→ Ensure both `.imzML` AND `.ibd` files are uploaded to the same directory.

---

## 8. FAQ

**Q: Can I use IonScell with processed (centroided) imzML?**  
Yes. Set `mz_bin` to match your peak width (e.g. 0.01 Da for high-resolution data).

**Q: What is the difference between GMM and FCM clustering?**  
Both are soft clustering methods. GMM models each cluster as a Gaussian distribution — it works best when clusters have ellipsoidal shapes in UMAP space. FCM assigns partial memberships based on distance and works well for overlapping populations with less clear structure.

**Q: How do I choose the number of clusters?**  
Leave `n_clusters=None` — IonScell will automatically test 2 to `max_clusters` and select the best k by combining BIC (model fit) and silhouette score (cluster separation). You can verify by inspecting the UMAP.

**Q: Can I add my own lipids to the database?**  
Yes. Add entries to the `_LIPID_DB` list in `ionscell_pipeline.py` following the format:
`("Name", "Formula", neutral_mass_Da, "class", "subclass", "HMDB_ID", "LIPIDMAPS_ID", "pos/neg/both")`

**Q: My dataset has > 10 000 cells — is IonScell fast enough?**  
Yes. Use `memmap_path` for the data cube, and optionally subsample UMAP with `n_neighbors=10`. GMM clustering with 10 000 cells typically takes < 2 minutes on a standard laptop.

**Q: Can I process multiple samples together?**  
Run each sample through Steps 1–7 independently, then concatenate `pipe.last_cell_spectra` DataFrames (add a `sample_id` column first), and run Steps 8–14 on the combined DataFrame.

---

*IonScell — Single-Cell Mass Spectrometry Imaging Pipeline*  
*For questions and bug reports: GitHub Issues*

# IonScell Desktop Interface — User Manual

<p align="center"><em>Version 1.0 · Advanced GUI for Single-Cell Mass Spectrometry Imaging</em></p>

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Installation & First Launch](#2-installation--first-launch)
3. [Interface Overview](#3-interface-overview)
4. [Module-by-Module Guide](#4-module-by-module-guide)
   - 4.1 File & Project Panel
   - 4.2 Data Loading
   - 4.3 Data Cube Builder
   - 4.4 TIC Viewer
   - 4.5 Thresholding & Segmentation
   - 4.6 Cell Manager
   - 4.7 Spectral Extraction & QC
   - 4.8 Soft Clustering
   - 4.9 Clone Manager
   - 4.10 Ion Maps & Spatial Viewer
   - 4.11 Differential Analysis
   - 4.12 Lipid Annotation
   - 4.13 Export Manager
5. [Keyboard Shortcuts](#5-keyboard-shortcuts)
6. [Comparison: Interface vs Notebook](#6-comparison-interface-vs-notebook)
7. [Troubleshooting](#7-troubleshooting)

---

## 1. Introduction


The **IonScell Desktop Interface** is a standalone graphical application that exposes the full `SCMSIPipeline` engine through an interactive GUI — no Python knowledge required.

⚠️ The desktop application is **not distributed in this GitHub repository**.  
It must be downloaded separately from the official distribution link.

📥 Download IonScell GUI:  
https://nextcloud.univ-lille.fr/index.php/f/384061741

---

## 2. Installation & First Launch

### Download the Interface (Windows)

The IonScell GUI is distributed separately from the GitHub repository.

1. Download the latest release archive from:  
   https://nextcloud.univ-lille.fr/index.php/f/384061741
2. Extract the archive to a folder **without spaces** in the path  
   Example: `C:\ionscell\`
3. Inside the extracted folder, double-click **ionScell.exe**

> The GitHub repository only contains the pipeline and notebooks.  
> The executable is **not included in the repo**.

### First launch

On first launch, IonScell will:

1. Verify the `_internal/` folder is present alongside the executable  
2. Display the Welcome screen with a link to sample data  
3. Offer to create a default workspace folder (`Documents/IonScell_Projects/`)

⚠️ Do not move the executable away from the `_internal/` folder.

---

## 3. Interface Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│  IonScell  v6c           [File] [Analysis] [View] [Help]            │
├──────────────┬──────────────────────────────────────┬───────────────┤
│              │                                      │               │
│  NAVIGATOR   │         MAIN VIEWER                  │  INSPECTOR    │
│              │                                      │               │
│  ○ Project   │  ┌──────────────────────────────┐   │  Parameters   │
│  ○ Files     │  │                              │   │  ─────────    │
│  ○ Steps     │  │   Interactive image / plot   │   │  Results      │
│  ○ Results   │  │                              │   │  ─────────    │
│              │  └──────────────────────────────┘   │  Notes        │
│  PIPELINE    │                                      │               │
│  STATUS      │  ─────────── CONSOLE ────────────   │               │
│  ✅ Loaded   │  [14:32] imzML loaded (12 450 px)   │               │
│  ✅ Cube     │  [14:33] Datacube built 110×95×7500 │               │
│  … Segment   │  [14:33] Auto threshold: 8 432.1    │               │
└──────────────┴──────────────────────────────────────┴───────────────┘
│  [Step 1: Load] [Step 2: Cube] [Step 3: TIC] … [Step 12: Export]    │
└─────────────────────────────────────────────────────────────────────┘
```

The interface has three zones:

| Zone | Description |
|------|-------------|
| **Navigator (left)** | Project tree, pipeline status, step list |
| **Main Viewer (centre)** | All interactive figures, plots, tables |
| **Inspector (right)** | Context-sensitive parameters and results summary |
| **Step Bar (bottom)** | Quick-access buttons for each pipeline step |
| **Console (bottom centre)** | Live log output |

---

## 4. Module-by-Module Guide

### 4.1 File & Project Panel

**Creating a new project:**
`File → New Project` (or `Ctrl+N`)  
Enter a project name and choose a workspace folder. IonScell creates a structured directory:

```
MyProject/
├── data/           ← place .imzML + .ibd here
├── results/        ← CSV, imzML exports
├── figures/        ← saved plots
└── project.json    ← pipeline state (auto-saved)
```

**Opening an existing project:**
`File → Open Project` (or `Ctrl+O`) → select `project.json`

IonScell restores all pipeline steps performed in the previous session, including segmentation, clustering assignments, and annotation results.

**Auto-save:** The project state is saved automatically every 5 minutes. You can disable this in `Edit → Preferences → Auto-save`.

---

### 4.2 Data Loading

**Step 1 — Load imzML**

1. Click **[Step 1: Load]** in the Step Bar, or `Analysis → Load imzML`
2. In the Inspector, click **Browse** and navigate to your `.imzML` file
3. Ensure the matching `.ibd` file is in the **same folder**
4. Click **Load**

The Console will confirm loading:
```
[INFO] imzML loaded: 12 450 spectral pixels
[INFO] Coordinate range: x=[1,110], y=[1,95]
[INFO] Format: continuous / processed
```

**Inspector fields after loading:**
- Number of pixels
- Coordinate range
- Detected m/z range
- File size
- Instrument metadata (if present in imzML)

---

### 4.3 Data Cube Builder

**Step 2 — Build Cube**

1. Click **[Step 2: Cube]**
2. Configure parameters in the Inspector:

| Parameter | Description | Typical value |
|-----------|-------------|--------------|
| m/z min | Lower bound | 350 |
| m/z max | Upper bound | 1100 |
| Bin mode | Fixed / Adaptive | Fixed |
| Bin width | Fixed bin size (Da) | 0.1 |
| Resolving power | For adaptive mode | 20 000 |
| Normalization | TIC / RMS / MEDIAN / None | TIC |

3. Click **Build Cube**

A progress bar shows spectral binning progress. For large files, enable **Memory-efficient mode** (tick box) to use disk-mapped arrays.

**After build:**
- Main Viewer shows a colour-coded TIC preview
- Inspector shows cube dimensions and TIC statistics
- Step 2 status turns ✅ green

---

### 4.4 TIC Viewer

**Step 3 — Visualise TIC**

The TIC viewer is an interactive Plotly figure embedded in the Main Viewer.

**Controls:**
- **Scroll / pinch** — zoom
- **Click + drag** — pan
- **Hover** — shows pixel coordinates and intensity
- **Double-click** — reset zoom

**Inspector options:**
- `Gaussian σ` — smoothing level (0 = raw)
- `Colormap` — choose from magma, viridis, hot, gray, …
- `Zoom factor` — upscale display (does not affect analysis)

**Save figure:** Click the camera icon in the Plotly toolbar, or `View → Save Figure as PNG/SVG`.

---

### 4.5 Thresholding & Segmentation

**Step 4 — Threshold & Segment**

This step combines adaptive thresholding and cell segmentation.

**Sub-step A: Thresholding**

1. In the Inspector, select **Threshold method**: `Auto` (recommended), or choose manually from the dropdown
2. Adjust `Min cell fraction` / `Max cell fraction` if needed
3. Click **Compute Threshold**

The Main Viewer shows the TIC histogram with all candidate thresholds overlaid. The selected threshold is highlighted. A preview binary mask appears below.

**Sub-step B: Segmentation**

Configure segmentation parameters:

| Parameter | Description |
|-----------|-------------|
| Method | Watershed (default) or DFS |
| Smoothing σ | Gaussian pre-smoothing |
| Footprint | Local maxima neighbourhood (px) |
| Min cell size | Pixels or µm (toggle unit) |
| Max cell size | Upper limit (optional) |
| Pixel size (µm) | Physical pixel size |

Click **Run Segmentation**.

The Main Viewer shows:
- TIC heatmap in the background
- Red cell contours overlaid
- Each cell labelled with its ID (zoom in to read labels)

**Live preview:** Adjust `Smoothing σ` with the slider and the contour preview updates in real time (no recomputation needed).

---

### 4.6 Cell Manager

**Removing artefact cells:**

1. In the Main Viewer (segmentation figure), hover over a suspicious cell to see its ID in the tooltip
2. Right-click → **Mark as artefact**  
   OR  
   Type the cell ID directly into the Inspector field **Cell IDs to remove** and click **Remove**
3. The contour turns grey and the cell is excluded from all subsequent steps

**Undo removal:** `Edit → Undo` (Ctrl+Z), or click **Restore Cell** in the Inspector and enter the cell ID.

**Cell statistics table:**

Below the segmentation figure, a sortable table shows all cells with:
- Cell ID, area (px and µm²), centroid coordinates
- TIC intensity at centroid
- Status (Active / Removed)

Click any column header to sort. Click a row to highlight the corresponding cell in the image.

---

### 4.7 Spectral Extraction & QC

**Step 6 — Extract Spectra**

1. Choose aggregation: **Mean** (default), Median, or Sum
2. Click **Extract Spectra**

**Step 7 — Quality Control**

Click **Compute QC Metrics**. The Main Viewer switches to a 2×3 panel showing distributions of 6 QC metrics:

| Metric | Meaning |
|--------|---------|
| SNR | Signal-to-noise ratio |
| Sparsity | Fraction of non-zero m/z bins |
| Dynamic range | log₁₀(max/5th percentile) |
| Total intensity | Sum of all bins |
| Spectral entropy | Information content |
| N peaks | Bins above median intensity |

**Filter controls (Inspector):**
- `SNR threshold` — slider (default 2.0)
- `Intensity percentile` — remove bottom N% by intensity
- `Sparsity min` — minimum fill fraction
- `Auto-remove flagged` — toggle

Green bars = passing cells, red bars = flagged cells.

Click **Apply Filter** to confirm removal, or **Preview Only** to see the effect without removing.

---

### 4.8 Soft Clustering

**Step 8 — Clustering**

Choose clustering algorithm:
- **GMM** (Gaussian Mixture Models) — recommended
- **Fuzzy C-Means** — for heavily overlapping populations

**GMM parameters (Inspector):**

| Parameter | Description | Default |
|-----------|-------------|---------|
| N clusters | Auto / Manual | Auto |
| Max clusters | Upper bound for auto | 8 |
| Covariance | tied / full / diag | tied |
| UMAP neighbours | n_neighbors | 15 |
| UMAP min dist | min_dist | 0.1 |
| Random seed | Reproducibility | 42 |

Click **Run Clustering**.

**Main Viewer shows two panels:**
- Left: UMAP scatter coloured by clone (opacity = confidence)
- Right: Spatial image with clone-coloured contours

**Live re-clustering:** Change `N clusters` with the spinner and click **Re-run** — results update within seconds.

---

### 4.9 Clone Manager

**Clone overview table:**

After clustering, the Inspector shows a table with one row per clone:

| Column | Description |
|--------|-------------|
| Clone ID | 0-indexed integer |
| N cells | Cell count |
| Mean confidence | Average soft assignment probability |
| CV (variability) | Within-clone coefficient of variation |
| Status | GOOD / QUESTIONABLE / ARTIFACT |
| Actions | Remove / Merge / Rename |

**Merging clones:**
1. Hold `Ctrl` and click two or more clone rows in the table
2. Click **Merge Selected**

**Removing a clone:**
1. Click the clone row
2. Click **Remove Clone** — all its cells are excluded

**Renaming:**
Double-click the clone name cell and type a label (e.g. "Lipid-rich cells").

**Clone quality report:**  
Click **Generate Quality Report** to open a detailed HTML report in your browser summarising confidence, homogeneity, and separability metrics for all clones.

---

### 4.10 Ion Maps & Spatial Viewer

**Single-ion map:**

1. Navigate to `Analysis → Ion Maps → Single Ion`
2. Type or search an m/z value in the Inspector (e.g. `760.59`)
3. Adjust tolerance (Da) and smoothing
4. Click **Show**

The Main Viewer shows the ion spatial distribution masked to cell pixels, with cell contours overlaid.

**Multi-ion overlay:**

1. `Analysis → Ion Maps → Multi-Ion Overlay`
2. Drag m/z values from the **Differential Results** table into the overlay list  
   OR type values manually (one per line)
3. Click a colour swatch next to each ion to change its colour
4. Adjust opacity with the slider
5. Click **Render Overlay**

Up to 6 ions can be overlaid simultaneously. Each ion is mapped to a distinct colour, and the composite image uses additive blending.

**Clone distribution plots:**

1. `Analysis → Clone Distribution`
2. Select ion (m/z) or QC metric from the dropdown
3. Choose plot type: Violin, Box, or Both
4. Toggle scatter points on/off
5. Click **Plot**

---

### 4.11 Differential Analysis

**Step 10 — Differential Analysis**

1. Click **[Step 10: Diff]** or `Analysis → Differential Analysis`
2. Configure in Inspector:

| Parameter | Default | Description |
|-----------|---------|-------------|
| Statistical test | Kruskal-Wallis | or ANOVA |
| FDR threshold | 0.05 | Adjusted p-value cutoff |
| Fold-change threshold | 1.5 | Minimum FC |
| Top N ions | 30 | Ions to highlight |

3. Click **Run Analysis**

**Main Viewer shows:**
- **Volcano plot** — each point is one m/z bin; significant ions colour-coded by best clone
- **Z-score heatmap** — top N ions × cells (sorted by clone)

**Interactivity:**
- Hover on any point in the volcano plot to see m/z, adjusted p-value, fold-change, and best clone
- Click a point in the volcano to jump directly to that ion's spatial map
- Click a row in the heatmap to open the ion distribution for that m/z

**Export results:**  
Click **Export Table** to save the full differential results as CSV.

---

### 4.12 Lipid Annotation

**Step 11 — Annotation**

1. Click **[Step 11: Annotate]**
2. Configure:

| Parameter | Description | Default |
|-----------|-------------|---------|
| Ionisation mode | pos / neg / both | pos |
| ppm tolerance | Mass accuracy | 5.0 |
| Top N ions | Annotate top N from DA | 30 |
| Source | Top markers / Custom list | Top markers |

3. Click **Annotate**

**Results panel:**

A searchable, sortable table appears with columns:
- m/z query, lipid name, adduct, ppm error, lipid class, subclass, HMDB ID, confidence score, best clone, biological pathways

**Lipid class enrichment heatmap:**  
Below the table, a heatmap shows lipid class counts per clone (normalised by column maximum).

**Free-text search:**  
Type in the search box (e.g. "ceramide", "PC 34") to filter the annotation table in real time.

**Manual lookup:**  
Type any m/z in the Inspector search field and click **Lookup** to annotate a single value independently.

---

### 4.13 Export Manager

**Step 12 — Export**

`File → Export` or click **[Step 12: Export]**

**Available exports:**

| Format | Content | Use case |
|--------|---------|----------|
| CSV (full) | Metadata + m/z intensity matrix | Statistical analysis in R/Python |
| CSV (metadata) | Cell ID, area, centroid, clone, QC | Quick overview |
| CSV (ions only) | m/z intensity matrix only | Machine learning input |
| imzML | Segmented cell pixels | MSI downstream tools |
| HTML Report | All figures + tables | Sharing with collaborators |
| Figure pack | All plots as PNG + SVG | Publication submission |

**HTML Report:**  
Generates a self-contained HTML file with all visualisations (TIC, segmentation, UMAP, spectra, volcano, annotation). Open in any browser — no internet required.

**Figure export settings:**
- Resolution: 150 / 300 / 600 DPI
- Format: PNG, SVG, TIFF, PDF
- Background: White / Black / Transparent

---

## 5. Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| `Ctrl+N` | New project |
| `Ctrl+O` | Open project |
| `Ctrl+S` | Save project |
| `Ctrl+Z` | Undo last action |
| `Ctrl+Y` | Redo |
| `Ctrl+E` | Export |
| `F5` | Re-run current step |
| `F11` | Toggle fullscreen |
| `Ctrl+1` … `Ctrl+9` | Jump to step 1…9 |
| `Ctrl+0` | Jump to step 10 (Diff) |
| `Space` | Play/pause live preview |
| `Delete` | Remove selected cell(s) in Cell Manager |
| `Escape` | Cancel running computation |

---

## 6. Comparison: Interface vs Notebook

| Feature | Desktop Interface | Jupyter Notebook |
|---------|-----------------|-----------------|
| Python required | ❌ No | ✅ Yes |
| Click-to-remove cells | ✅ Yes | ❌ Manual ID entry |
| Real-time preview | ✅ Yes | ❌ No |
| Custom code | ❌ No | ✅ Full Python access |
| Google Colab support | ❌ No | ✅ Yes |
| Reproducible scripts | ✅ | ✅ |
| Pipeline state save | ✅ Auto | Manual |
| Batch processing | ❌ (v1.0) | ✅ Loop in Python |

**Recommendation:**  
Use the **interface** for rapid QC, interactive exploration, and sharing with collaborators. Use the **notebook** for publication-quality reproducible analyses and custom downstream processing.

---

## 7. Troubleshooting

**Application fails to start**  
→ Verify `_internal/` is in the same folder as `ionScell.exe`.  
→ Try running as Administrator (right-click → Run as administrator).

**Blank white window on startup**  
→ Update your graphics drivers.  
→ Try disabling hardware acceleration: edit `_internal/qt.conf`, add `HardwareAcceleration=false`.

**"imzML parser error"**  
→ Ensure the `.ibd` file is present alongside `.imzML` with the same base name.  
→ Check the imzML file is not corrupted (try opening in MSiReader as a sanity check).

**Segmentation produces thousands of tiny fragments**  
→ Increase `Smoothing σ` to 1.5–2.0 before segmentation.  
→ adapat `Min cell size` and `max cell size` 
→ Try lowering the TIC threshold manually.

**UMAP takes > 5 minutes**  
→ Reduce `n_neighbors` to 10.  
→ This is normal for > 5 000 cells without GPU. The interface will show a progress bar.

**Annotations show no results**  
→ Check that your m/z range overlaps the lipid database (350–1 100 Da).  
→ Increase ppm tolerance to 20 ppm for lower-resolution data.  
→ Switch ionisation mode (`pos` ↔ `neg`) to match your acquisition.

---

*IonScell Desktop Interface v1.0*  
*For bug reports and feature requests: GitHub Issues*

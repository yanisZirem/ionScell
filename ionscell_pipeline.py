"""
ionScell Pipeline — SCMSIPipeline
===================================
Single-Cell Mass Spectrometry Imaging (SCMSI) analysis pipeline.

This module provides the SCMSIPipeline class and the LipidAnnotator class.
Import this module in your Jupyter notebook or Python script:

    from ionscell_pipeline import SCMSIPipeline, LipidAnnotator

Author: Yanis Zirem, e-mail : yanis.zirem@univ-lille.fr/ yanis.zirem2016@gmail.com
Version: 1.0
"""

# ─── Standard library ────────────────────────────────────────────────────────
import os
import sys as _sys
import io as _io
import warnings

# Patch stdout/stderr before tqdm loads — required for PyInstaller console=False
if _sys.stdout is None:
    _sys.stdout = _io.StringIO()
if _sys.stderr is None:
    _sys.stderr = _io.StringIO()

# ─── Third-party ─────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
import plotly.colors as pc
from pyimzml.ImzMLParser import ImzMLParser
from pyimzml.ImzMLWriter import ImzMLWriter
from scipy.ndimage import gaussian_filter, zoom as scipy_zoom, label as ndi_label
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from skimage import measure
from skimage.filters import threshold_otsu
from tqdm import tqdm


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN PIPELINE CLASS
# ══════════════════════════════════════════════════════════════════════════════

class SCMSIPipeline:
    """
    Enhanced pipeline for Single-Cell Mass Spectrometry Imaging (SCMSI) with:
    - Fixed and adaptive m/z binning
    - Adaptive TIC thresholding (GMM, Otsu, Triangle, Li, Yen, …)
    - Watershed and DFS-based cell segmentation
    - Soft clustering (GMM / Fuzzy C-Means)
    - Quality control for spectral and clone assessment
    - Differential analysis & volcano plot
    - Multi-ion spatial overlay
    - Lipid annotation (embedded LIPID MAPS / HMDB database)
    - Export to CSV and imzML
    """

    def __init__(self, imzml_path=None, verbose=True):
        self.imzml_path = imzml_path
        self.parser = None
        self.coords = None
        self.mz_axis = None
        self.data_cube = None
        self.TIC = None
        self.cell_labels = None
        self.last_cell_spectra = None
        self.verbose = verbose
        self._coord_min = (0, 0)
        self.pixel_size_um = 10
        self.qc_results = None          # QC metrics storage
        self.soft_memberships = None    # Soft clustering probabilities

    # ──────────────────────────────────────────────────────────────────────────
    #  DATA LOADING
    # ──────────────────────────────────────────────────────────────────────────

    def load_imzml(self, path=None):
        """Load an imzML file and initialise coordinates."""
        if path is not None:
            self.imzml_path = path
        if self.imzml_path is None:
            raise ValueError("No imzML path provided.")
        self.parser = ImzMLParser(self.imzml_path)
        self.coords = np.array(self.parser.coordinates)
        xs = self.coords[:, 0].astype(int)
        ys = self.coords[:, 1].astype(int)
        xmin, ymin = xs.min(), ys.min()
        self._coord_min = (xmin, ymin)
        if self.verbose:
            print(f"imzML loaded: {len(self.coords)} spectral pixels; coord min = {self._coord_min}")
        return self.parser

    # ──────────────────────────────────────────────────────────────────────────
    #  m/z AXIS CONSTRUCTION
    # ──────────────────────────────────────────────────────────────────────────

    def make_mz_axis(self, mz_min=600, mz_max=1100, mz_bin=0.1):
        """
        Build a fixed-grid m/z axis (uniform bin width).

        Recommended for exploratory clustering. For annotation workflows
        use build_adaptive_mz_axis() instead.

        Parameters
        ----------
        mz_min, mz_max : float  — m/z range (Da)
        mz_bin : float          — bin width (Da), default 0.1
        """
        self.mz_axis = np.arange(mz_min, mz_max + mz_bin / 2, mz_bin)
        self._mz_binning_mode = 'fixed'
        self._mz_bin_width = mz_bin
        if self.verbose:
            print(f"[BIN] Fixed grid: {mz_min}–{mz_max} Da, step={mz_bin} Da "
                  f"→ {len(self.mz_axis)} bins")
        return self.mz_axis

    def build_adaptive_mz_axis(self, mz_min=None, mz_max=None,
                                target_resolving_power=None,
                                ppm_bin=None,
                                peak_picking_n_scans=50,
                                min_bin_width=0.01,
                                max_bin_width=1.0,
                                min_peaks_for_adaptation=200):
        """
        Adaptive m/z binning — bin width scales with m/z to match the
        instrument's resolving power (constant ppm width).

        Recommended for annotation workflows (LIPID MAPS, HMDB matching).

        Parameters
        ----------
        target_resolving_power : int or float
            Instrument R (e.g. 20 000 for MALDI-TOF, 50 000 for Orbitrap).
        ppm_bin : float
            Fixed ppm bin width (alternative to resolving_power).
        mz_min, mz_max : float — If None, inferred from data.
        peak_picking_n_scans : int — Spectra sampled for auto-estimation.
        min_bin_width, max_bin_width : float — Clamps on bin width (Da).

        Returns
        -------
        np.ndarray  — Non-uniform m/z axis (bin centres).
        dict        — Diagnostics.
        """
        if self.parser is None:
            raise RuntimeError("Load imzML first (load_imzml()).")

        if mz_min is None or mz_max is None:
            n_sample = min(20, len(self.coords))
            idx_sample = np.random.choice(len(self.coords), n_sample, replace=False)
            all_mz = []
            for i in idx_sample:
                mzs, _ = self.parser.getspectrum(i)
                if len(mzs):
                    all_mz.extend(mzs)
            if not all_mz:
                raise RuntimeError("No spectral data found.")
            mz_min = mz_min or float(np.percentile(all_mz, 1))
            mz_max = mz_max or float(np.percentile(all_mz, 99))
            if self.verbose:
                print(f"[BIN] Auto range: {mz_min:.1f}–{mz_max:.1f} Da")

        if target_resolving_power is not None:
            mode_label = f"resolving_power={target_resolving_power:.0f}"
            def _bw(m):
                return np.clip(m / target_resolving_power, min_bin_width, max_bin_width)

        elif ppm_bin is not None:
            mode_label = f"ppm={ppm_bin}"
            def _bw(m):
                return np.clip(m * ppm_bin / 1e6, min_bin_width, max_bin_width)

        else:
            mode_label = "auto (data-driven)"
            n_scans = min(peak_picking_n_scans, len(self.coords))
            idx_s = np.random.default_rng(42).choice(len(self.coords), n_scans, replace=False)
            spacings_by_mz = []
            for i in idx_s:
                mzs, intens = self.parser.getspectrum(i)
                if len(mzs) < 5:
                    continue
                mzs = np.sort(np.asarray(mzs, dtype=np.float64))
                thr = np.median(intens)
                mzs_peaks = mzs[np.asarray(intens) > thr]
                if len(mzs_peaks) < 2:
                    continue
                diff = np.diff(mzs_peaks)
                for j, d in enumerate(diff):
                    center = (mzs_peaks[j] + mzs_peaks[j+1]) / 2
                    if mz_min <= center <= mz_max:
                        spacings_by_mz.append((center, d))

            if len(spacings_by_mz) >= min_peaks_for_adaptation:
                spacings_by_mz.sort()
                mz_c = np.array([s[0] for s in spacings_by_mz])
                sp_c = np.array([s[1] for s in spacings_by_mz])
                ppm_est = np.median(sp_c / mz_c * 1e6)
                ppm_est = float(np.clip(ppm_est, 2.0, 500.0))
                if self.verbose:
                    print(f"[BIN] Auto estimated resolving power "
                          f"≈ {1e6/ppm_est:.0f}  ({ppm_est:.1f} ppm/bin)")
                def _bw(m):
                    return np.clip(m * ppm_est / 1e6, min_bin_width, max_bin_width)
            else:
                if self.verbose:
                    print(f"[BIN] Insufficient peaks for auto-adaptation. "
                          f"Falling back to fixed 0.1 Da.")
                return self.make_mz_axis(mz_min=mz_min, mz_max=mz_max, mz_bin=0.1), \
                       {'mode': 'fixed_fallback', 'n_bins': len(self.mz_axis)}

        mz_edges = [mz_min]
        while mz_edges[-1] < mz_max:
            w = _bw(mz_edges[-1])
            mz_edges.append(mz_edges[-1] + w)

        mz_edges = np.array(mz_edges, dtype=np.float64)
        mz_centers = (mz_edges[:-1] + mz_edges[1:]) / 2
        mz_bw = np.diff(mz_edges)

        self.mz_axis = mz_centers.astype(np.float32)
        self._mz_bin_edges = mz_edges.astype(np.float64)
        self._mz_bin_widths = mz_bw.astype(np.float64)
        self._mz_binning_mode = 'adaptive'
        self._mz_bin_width = None

        diag = {
            'mode': mode_label,
            'n_bins': len(mz_centers),
            'min_bw': float(mz_bw.min()),
            'max_bw': float(mz_bw.max()),
            'ppm_equiv_min': float(mz_bw.min() / mz_min * 1e6),
            'ppm_equiv_max': float(mz_bw.max() / mz_max * 1e6),
        }

        if self.verbose:
            print(f"[BIN] Adaptive axis: mode={mode_label}, "
                  f"{len(mz_centers)} bins, "
                  f"bin width {mz_bw.min():.4f}–{mz_bw.max():.4f} Da "
                  f"({diag['ppm_equiv_min']:.0f}–{diag['ppm_equiv_max']:.0f} ppm)")
        return self.mz_axis, diag

    # ──────────────────────────────────────────────────────────────────────────
    #  DATA CUBE
    # ──────────────────────────────────────────────────────────────────────────

    def build_datacube(self, mz_min=None, mz_max=None, mz_bin=0.1,
                       binning_mode='fixed', mem_efficient=False,
                       memmap_path=None, tic_only=False, use_tqdm=True,
                       normalization=None, chunk_size=5000):
        """
        Build the (height × width × n_bins) data cube from the imzML file.

        Parameters
        ----------
        mz_min, mz_max : float   — m/z range.
        mz_bin : float           — Fixed bin width in Da (fixed mode only).
        binning_mode : str       — 'fixed' (default) or 'adaptive'.
        normalization : str|None — 'TIC', 'RMS', 'MEDIAN', or None.
        memmap_path : str|None   — Path for memory-mapped file (large datasets).
        tic_only : bool          — Only compute TIC (faster).
        chunk_size : int         — Spectra per memory chunk.
        """
        if self.parser is None:
            raise RuntimeError("ImzML not loaded. Execute load_imzml().")

        if self.mz_axis is None:
            if binning_mode == 'adaptive':
                if self.verbose:
                    print("[BIN] Building adaptive m/z axis …")
                self.build_adaptive_mz_axis(mz_min=mz_min, mz_max=mz_max)
            else:
                if mz_min is None or mz_max is None:
                    raise ValueError("Define mz_min/mz_max or call make_mz_axis().")
                self.make_mz_axis(mz_min=mz_min, mz_max=mz_max, mz_bin=mz_bin)

        xs = self.coords[:, 0].astype(int)
        ys = self.coords[:, 1].astype(int)
        xmin, ymin = self._coord_min
        xs_rel = xs - xmin
        ys_rel = ys - ymin
        width  = xs_rel.max() + 1
        height = ys_rel.max() + 1
        n_bins = len(self.mz_axis)

        if self.verbose:
            mode_tag = getattr(self, '_mz_binning_mode', 'fixed')
            print(f"[CUBE] {width}×{height} px, {n_bins} bins "
                  f"({mode_tag} binning)  |  norm={normalization or 'None'}")

        if tic_only:
            tic = np.zeros((height, width), dtype=np.float32)
            it = tqdm(range(len(self.coords)), desc="TIC") if use_tqdm else range(len(self.coords))
            for idx in it:
                _, intens = self.parser.getspectrum(idx)
                tic[ys_rel[idx], xs_rel[idx]] += float(np.sum(intens))
            self.TIC = tic
            return self.TIC

        if memmap_path is not None:
            os.makedirs(os.path.dirname(memmap_path) or ".", exist_ok=True)
            data_cube = np.memmap(memmap_path, mode='w+', dtype=np.float32,
                                  shape=(height, width, n_bins))
            if self.verbose:
                print(f"[CUBE] memmap: {memmap_path}")
        else:
            data_cube = np.zeros((height, width, n_bins), dtype=np.float32)

        adaptive = getattr(self, '_mz_binning_mode', 'fixed') == 'adaptive'

        if adaptive:
            edges = self._mz_bin_edges
            def _bin_spectrum(mzs, intens):
                mzs = np.asarray(mzs, dtype=np.float64)
                intens = np.asarray(intens, dtype=np.float32)
                idx = np.searchsorted(edges, mzs, side='right') - 1
                valid = (idx >= 0) & (idx < n_bins)
                return idx[valid], intens[valid]
        else:
            bin_left  = float(self.mz_axis[0])
            bin_width = float((self.mz_axis[1] - self.mz_axis[0]) if n_bins > 1 else 1.0)
            def _bin_spectrum(mzs, intens):
                mzs = np.asarray(mzs, dtype=np.float32)
                intens = np.asarray(intens, dtype=np.float32)
                idx = np.floor((mzs - bin_left) / bin_width).astype(int)
                valid = (idx >= 0) & (idx < n_bins)
                return idx[valid], intens[valid]

        n_spectra = len(self.coords)
        it = tqdm(range(n_spectra), desc="Build datacube") if use_tqdm \
             else range(n_spectra)

        for idx in it:
            mzs, intens = self.parser.getspectrum(idx)
            if len(mzs) == 0:
                continue

            intens = np.array(intens, dtype=np.float32, copy=True)

            if normalization is not None:
                nm = normalization.upper()
                if nm == "TIC":
                    s = np.sum(intens)
                    if s > 0:
                        intens /= s
                elif nm == "RMS":
                    r = np.sqrt(np.mean(intens ** 2))
                    if r > 0:
                        intens /= r
                elif nm == "MEDIAN":
                    med = np.median(intens[intens > 0]) if np.any(intens > 0) else 1.0
                    intens /= max(med, 1e-12)
                else:
                    raise ValueError(f"Unknown normalization: {normalization}. "
                                     "Choose TIC, RMS, or MEDIAN.")

            bin_idx, bin_int = _bin_spectrum(mzs, intens)
            if len(bin_idx) == 0:
                continue
            bx, by = xs_rel[idx], ys_rel[idx]
            np.add.at(data_cube[by, bx], bin_idx, bin_int)

        self.TIC = data_cube.sum(axis=2)
        self.data_cube = data_cube

        if self.verbose:
            print(f"[CUBE] Done. TIC range: "
                  f"{self.TIC.min():.2e}–{self.TIC.max():.2e}")
        return self.data_cube

    # ──────────────────────────────────────────────────────────────────────────
    #  TIC VISUALISATION
    # ──────────────────────────────────────────────────────────────────────────

    def show_TIC_plotly(self, sigma=0.5, zoom_factor=1.0, cmap='magma', figsize=(800, 800)):
        """Display TIC image as interactive Plotly heatmap."""
        if self.TIC is None:
            raise RuntimeError("TIC missing. Run build_datacube().")

        tic_s = gaussian_filter(self.TIC, sigma=sigma)
        tic_disp = scipy_zoom(tic_s, zoom_factor, order=3) if zoom_factor != 1.0 else tic_s

        fig = go.Figure()
        fig.add_trace(go.Heatmap(
            z=tic_disp,
            colorscale=cmap,
            colorbar=dict(title="TIC Intensity"),
            hovertemplate="x=%{x}, y=%{y}<br>Intensity=%{z:.2e}<extra></extra>"
        ))
        fig.update_layout(
            title=f"TIC (Total Ion Current) – σ={sigma}, zoom={zoom_factor}x",
            width=figsize[0], height=figsize[1],
            xaxis=dict(showgrid=False, zeroline=False),
            yaxis=dict(scaleanchor="x", autorange='reversed'),
            template="plotly_dark"
        )
        fig.show()
        return tic_disp

    # ──────────────────────────────────────────────────────────────────────────
    #  ADAPTIVE TIC THRESHOLDING
    # ──────────────────────────────────────────────────────────────────────────

    def auto_TIC_threshold(self, method='auto', percentile=75,
                           show_diagnostics=True, min_cell_fraction=0.05,
                           max_cell_fraction=0.60):
        """
        Adaptive TIC thresholding — automatically selects the best strategy
        based on the actual intensity distribution.

        Methods: 'auto' (recommended), 'gmm', 'triangle', 'otsu', 'li',
                 'yen', 'isodata', 'percentile', 'local_gmm'.

        Returns
        -------
        float          — Best threshold value.
        dict           — Diagnostics (only when method='auto').
        """
        if self.TIC is None:
            raise RuntimeError("TIC missing. Run build_datacube().")

        from skimage.filters import (threshold_otsu, threshold_li,
                                     threshold_yen, threshold_isodata,
                                     threshold_triangle)
        from sklearn.mixture import GaussianMixture

        tic = self.TIC.copy()
        tic_flat = tic[tic > 0].flatten()
        if len(tic_flat) == 0:
            raise RuntimeError("TIC is all-zero — check datacube build.")

        log_tic = np.log1p(tic_flat)
        total_px = tic.size

        def _score_threshold(t, log_space=False):
            mask = tic > t
            frac = mask.sum() / total_px
            if frac < min_cell_fraction or frac > max_cell_fraction:
                return 0.0, frac
            fg = tic[mask].flatten()
            bg = tic[~mask].flatten()
            if len(fg) == 0 or len(bg) == 0:
                return 0.0, frac
            w_fg = len(fg) / total_px
            w_bg = len(bg) / total_px
            mu_fg, mu_bg = fg.mean(), bg.mean()
            mu_t = tic.mean()
            bcv = w_fg * (mu_fg - mu_t)**2 + w_bg * (mu_bg - mu_t)**2
            tcv = tic.var() + 1e-9
            sep = bcv / tcv
            optimal_frac = 0.15
            frac_score = np.exp(-((frac - optimal_frac)**2) / (2 * 0.12**2))
            score = 0.6 * sep + 0.4 * frac_score
            return float(score), float(frac)

        candidates = {}

        for name, fn in [('otsu', threshold_otsu), ('triangle', threshold_triangle),
                         ('li', threshold_li), ('yen', threshold_yen),
                         ('isodata', threshold_isodata)]:
            try:
                t = fn(tic)
                s, f = _score_threshold(t)
                candidates[name] = {'threshold': t, 'score': s, 'cell_frac': f}
            except Exception:
                pass

        t = np.percentile(tic_flat, percentile)
        s, f = _score_threshold(t)
        candidates['percentile'] = {'threshold': t, 'score': s, 'cell_frac': f}

        try:
            for n_comp in [2, 3]:
                log_vals = log_tic.reshape(-1, 1)
                gmm = GaussianMixture(n_components=n_comp, covariance_type='full',
                                      random_state=42, n_init=5, max_iter=300)
                gmm.fit(log_vals)
                means = gmm.means_.flatten()
                sorted_idx = np.argsort(means)
                mu_bg_log = means[sorted_idx[0]]
                mu_fg_log = means[sorted_idx[min(1, n_comp-1)]]
                log_grid = np.linspace(mu_bg_log, mu_fg_log, 500)
                proba_grid = gmm.predict_proba(log_grid.reshape(-1, 1))
                bg_comp = sorted_idx[0]
                fg_comp = sorted_idx[min(1, n_comp-1)]
                crossing = log_grid[np.argmin(np.abs(
                    proba_grid[:, bg_comp] - proba_grid[:, fg_comp]))]
                t_log = np.expm1(crossing)
                s, f = _score_threshold(t_log)
                key = f'gmm{n_comp}'
                candidates[key] = {'threshold': t_log, 'score': s,
                                   'cell_frac': f, 'n_components': n_comp}
        except Exception as e:
            if self.verbose:
                print(f"[WARN] GMM thresholding failed: {e}")

        if method == 'auto':
            valid = {k: v for k, v in candidates.items() if v['score'] > 0}
            if not valid:
                valid = candidates
            best_method = max(valid, key=lambda k: valid[k]['score'])
            best_threshold = valid[best_method]['threshold']
        elif method in candidates:
            best_method = method
            best_threshold = candidates[method]['threshold']
        elif method == 'gmm':
            key = 'gmm2' if 'gmm2' in candidates else 'gmm3'
            best_method = key
            best_threshold = candidates[key]['threshold']
        else:
            raise ValueError(f"Unknown method '{method}'.")

        if self.verbose:
            print(f"\n{'='*55}")
            print(f"  ADAPTIVE TIC THRESHOLDING  (method={method})")
            print(f"{'='*55}")
            for k, v in sorted(candidates.items(), key=lambda x: -x[1]['score']):
                marker = " ◀ SELECTED" if k == best_method else ""
                print(f"  {k:12s}  thresh={v['threshold']:9.2f}  "
                      f"score={v['score']:.3f}  "
                      f"cell_frac={v['cell_frac']:.3f}{marker}")
            print(f"{'='*55}\n")

        self._last_threshold_diagnostics = candidates

        if show_diagnostics:
            try:
                nbins = min(300, max(80, len(tic_flat) // 200))
                counts, edges = np.histogram(np.log1p(tic_flat), bins=nbins)
                bin_centers = (edges[:-1] + edges[1:]) / 2
                fig = go.Figure()
                fig.add_trace(go.Bar(
                    x=np.expm1(bin_centers), y=counts,
                    marker_color='rgba(100,149,237,0.6)',
                    name='TIC distribution',
                ))
                colors_map = {
                    'otsu': '#ff6b6b', 'triangle': '#ffd93d',
                    'li': '#6bcb77', 'yen': '#4d96ff',
                    'isodata': '#ff922b', 'percentile': '#cc5de8',
                    'gmm2': '#f06595', 'gmm3': '#f06595',
                }
                for k, v in candidates.items():
                    dash = 'solid' if k == best_method else 'dash'
                    width = 3 if k == best_method else 1.5
                    fig.add_vline(
                        x=v['threshold'],
                        line=dict(color=colors_map.get(k, '#aaa'), width=width, dash=dash),
                        annotation_text=f"{k} ({v['cell_frac']:.0%})",
                        annotation_font_size=10,
                        annotation_font_color=colors_map.get(k, '#aaa')
                    )
                fig.update_layout(
                    title=f"TIC Distribution — Threshold Candidates<br>"
                          f"<sup>Selected: <b>{best_method}</b>  "
                          f"threshold={best_threshold:.1f}</sup>",
                    xaxis_title="TIC intensity", yaxis_title="Pixel count",
                    template="plotly_dark", xaxis_type='log', height=450, bargap=0
                )
                fig.show()
            except Exception as e:
                if self.verbose:
                    print(f"[WARN] Diagnostics plot failed: {e}")

        if method == 'auto':
            return best_threshold, candidates
        return best_threshold

    # ──────────────────────────────────────────────────────────────────────────
    #  CELL SEGMENTATION
    # ──────────────────────────────────────────────────────────────────────────

    def watershed_or_dfs_from_mask(self, sigma=0, binary_mask=None, TIC_threshold=0,
                                    footprint=(3, 3), min_cell_size=2, max_cell_size=None,
                                    pixel_size_um=10, size_unit='px', remove_isolated=True,
                                    show_plotly=True, crop=None, verbose=True,
                                    colorscale='magma', color='red', width=0.8,
                                    method='watershed', intensity_tolerance=0.15):
        """
        Segment individual cells using Watershed or DFS algorithm.

        Parameters
        ----------
        sigma : float         — Gaussian smoothing before segmentation.
        binary_mask : array   — Pre-computed binary mask (or None to auto-threshold).
        TIC_threshold : float — Threshold if binary_mask is None.
        footprint : tuple     — Neighbourhood size for local maxima detection.
        min_cell_size : int   — Minimum cell size (in pixels or µm, see size_unit).
        max_cell_size : int   — Maximum cell size (None = no limit).
        pixel_size_um : float — Physical pixel size in µm.
        size_unit : str       — 'px' or 'um'.
        method : str          — 'watershed' (recommended) or 'dfs'.
        """
        self.pixel_size_um = pixel_size_um

        if self.TIC is None:
            raise RuntimeError("TIC missing. Run build_datacube() first.")

        tic = self.TIC.copy()
        if binary_mask is None:
            binary_mask = tic > TIC_threshold
            if verbose:
                print(f"[INFO] Mask created from TIC > {TIC_threshold}")

        if crop is not None:
            binary_mask = binary_mask[:crop[0], :crop[1]]
            tic = tic[:crop[0], :crop[1]]

        tic = gaussian_filter(tic, sigma=sigma)
        height_map = -tic
        height_map[~binary_mask] = 0

        if method.lower() == 'watershed':
            local_max_coords = peak_local_max(tic, labels=binary_mask, footprint=np.ones(footprint))
            local_max = np.zeros_like(tic, dtype=bool)
            if local_max_coords.size:
                local_max[tuple(local_max_coords.T)] = True
            markers = ndi_label(local_max)[0]
            labels_result = watershed(height_map, markers, mask=binary_mask)
            if verbose:
                print("[INFO] Watershed segmentation completed.")

        elif method.lower() == 'dfs':
            if verbose:
                print("[INFO] Using DFS-based segmentation.")

            rows, cols = binary_mask.shape
            labels_result = np.zeros_like(binary_mask, dtype=int)
            current_label = 0

            tic_norm = (tic - np.min(tic[binary_mask])) / (np.max(tic[binary_mask]) - np.min(tic[binary_mask]) + 1e-9)

            def dfs(r, c, label_id, seed_intensity):
                stack = [(r, c)]
                while stack:
                    x, y = stack.pop()
                    if (0 <= x < rows and 0 <= y < cols and
                        binary_mask[x, y] and labels_result[x, y] == 0):
                        if abs(tic_norm[x, y] - seed_intensity) <= intensity_tolerance:
                            labels_result[x, y] = label_id
                            for dx in [-1, 0, 1]:
                                for dy in [-1, 0, 1]:
                                    if dx != 0 or dy != 0:
                                        stack.append((x + dx, y + dy))

            local_max_coords = peak_local_max(tic, labels=binary_mask, footprint=np.ones(footprint))
            if verbose:
                print(f"[INFO] DFS seeding with {len(local_max_coords)} local maxima.")

            for coord in local_max_coords:
                i, j = coord
                if binary_mask[i, j] and labels_result[i, j] == 0:
                    current_label += 1
                    dfs(i, j, current_label, tic_norm[i, j])

            if verbose:
                print(f"[INFO] DFS segmentation completed: {current_label} initial regions found.")

            unassigned = np.argwhere((binary_mask == 1) & (labels_result == 0))
            if len(unassigned) > 0:
                from scipy.spatial import cKDTree
                assigned_coords = np.argwhere(labels_result > 0)
                assigned_labels = labels_result[labels_result > 0]
                tree = cKDTree(assigned_coords)
                for x, y in unassigned:
                    _, idx = tree.query((x, y))
                    labels_result[x, y] = assigned_labels[idx]
        else:
            raise ValueError(f"Unknown method '{method}'.")

        pixel_area_um2 = pixel_size_um ** 2
        min_size_px = min_cell_size
        max_size_px = max_cell_size

        if size_unit == 'um':
            if min_cell_size is not None:
                r_min = min_cell_size / 2
                min_size_px = np.pi * r_min ** 2 / pixel_area_um2
            if max_cell_size is not None:
                r_max = max_cell_size / 2
                max_size_px = np.pi * r_max ** 2 / pixel_area_um2

        def filter_cells_by_size(labels):
            new_labels = labels.copy()
            for lab in np.unique(labels):
                if lab == 0:
                    continue
                mask_lab = labels == lab
                cell_size_px = np.sum(mask_lab)
                if cell_size_px < min_size_px:
                    if remove_isolated:
                        new_labels[mask_lab] = 0
                elif max_size_px is not None and cell_size_px > max_size_px:
                    new_labels[mask_lab] = 0
            return new_labels

        labels_filtered = filter_cells_by_size(labels_result)
        self.cell_labels = labels_filtered

        if verbose:
            n_cells = len(np.unique(labels_filtered)) - 1
            print(f"{method.capitalize()} completed: {n_cells} cells detected.")

        if show_plotly:
            fig = go.Figure()
            fig.add_trace(go.Heatmap(
                z=tic, colorscale=colorscale, showscale=True,
                colorbar=dict(title="TIC Intensity"),
                hovertemplate="x=%{x}, y=%{y}<br>Intensity=%{z:.2e}<extra></extra>"
            ))
            for lab in np.unique(labels_filtered):
                if lab == 0:
                    continue
                mask_lab = labels_filtered == lab
                contours = measure.find_contours(mask_lab, 0.5)
                for contour in contours:
                    fig.add_trace(go.Scatter(
                        x=contour[:, 1], y=contour[:, 0], mode='lines',
                        line=dict(color=color, width=width), opacity=0.8,
                        hoverinfo='text', text=[f'Cell ID: {lab}'] * len(contour),
                        showlegend=False
                    ))
            fig.update_layout(
                title=f"Cell contours ({method.capitalize()})", width=800, height=800,
                yaxis=dict(scaleanchor="x", autorange='reversed'),
                xaxis=dict(showgrid=False, zeroline=False),
                template="plotly_dark", paper_bgcolor="black", plot_bgcolor="black",
                margin=dict(l=0, r=0, t=40, b=0)
            )
            fig.show()

        return labels_filtered

    # ──────────────────────────────────────────────────────────────────────────
    #  CELL MANAGEMENT
    # ──────────────────────────────────────────────────────────────────────────

    def remove_cells(self, remove_ids, show_plotly=True, verbose=True,
                     colorscale='jet', color='red', width=0.8, removed_color='gray'):
        """Remove specific cell IDs (artifacts) from the segmentation."""
        if self.cell_labels is None:
            raise RuntimeError("Cell segmentation missing.")

        cleaned_labels = self.cell_labels.copy()
        existing_ids = np.unique(cleaned_labels)
        valid_remove_ids = [rid for rid in remove_ids if rid in existing_ids]

        if len(valid_remove_ids) == 0:
            if verbose:
                print("[WARN] No valid cell IDs found to remove.")
            return cleaned_labels

        labels_original = cleaned_labels.copy()
        remove_mask = np.isin(cleaned_labels, valid_remove_ids)
        cleaned_labels[remove_mask] = 0
        self.cell_labels = cleaned_labels

        if verbose:
            remaining = np.unique(cleaned_labels)
            n_remaining = len(remaining[remaining > 0])
            print(f"[INFO] Removed {len(valid_remove_ids)} artefacts: {valid_remove_ids}")
            print(f"[INFO] {n_remaining} cells remaining.")

        if show_plotly:
            if not hasattr(self, "TIC") or self.TIC is None:
                raise RuntimeError("TIC missing.")

            tic = self.TIC.copy()
            fig = go.Figure()
            fig.add_trace(go.Heatmap(
                z=tic, colorscale=colorscale, showscale=True,
                colorbar=dict(title="TIC Intensity"),
            ))
            for rid in valid_remove_ids:
                mask_removed = labels_original == rid
                contours = measure.find_contours(mask_removed.astype(float), 0.5)
                for contour in contours:
                    fig.add_trace(go.Scatter(
                        x=contour[:, 1], y=contour[:, 0], mode='lines',
                        line=dict(color=removed_color, width=width), opacity=0.6,
                        hoverinfo='text', text=[f'Removed cell ID {rid}'] * len(contour),
                        showlegend=False
                    ))
            for lab in np.unique(cleaned_labels):
                if lab == 0:
                    continue
                mask_lab = cleaned_labels == lab
                contours = measure.find_contours(mask_lab.astype(float), 0.5)
                for contour in contours:
                    fig.add_trace(go.Scatter(
                        x=contour[:, 1], y=contour[:, 0], mode='lines',
                        line=dict(color=color, width=width), opacity=0.9,
                        hoverinfo='text', text=[f'Cell ID {lab}'] * len(contour),
                        showlegend=False
                    ))
            fig.update_layout(
                title="Segmentation after removing cells", width=800, height=800,
                yaxis=dict(scaleanchor="x", autorange='reversed'),
                xaxis=dict(showgrid=False, zeroline=False),
                template="plotly_dark", paper_bgcolor="black", plot_bgcolor="black",
                margin=dict(l=0, r=0, t=40, b=0)
            )
            fig.show()

        return cleaned_labels

    def remove_clone(self, clone_id, show_plotly=True, verbose=True,
                     colorscale='jet', color='cyan', removed_color='gray', width=0.8):
        """Remove all cells belonging to a given MSp/clone (by Class column)."""
        if self.last_cell_spectra is None or 'Class' not in self.last_cell_spectra.columns:
            raise RuntimeError("Run clustering first — 'Class' column missing.")
        if self.cell_labels is None:
            raise RuntimeError("Cell labels missing. Run segmentation first.")

        mask_clone = self.last_cell_spectra['Class'].astype(int) == int(clone_id)
        ids_to_remove = self.last_cell_spectra.loc[mask_clone, 'cell_id'].astype(int).tolist()

        if len(ids_to_remove) == 0:
            if verbose:
                print(f"[WARN] No cells found for clone {clone_id}.")
            return self.last_cell_spectra

        labels_original = self.cell_labels.copy()
        remove_mask = np.isin(self.cell_labels, ids_to_remove)
        self.cell_labels[remove_mask] = 0
        self.last_cell_spectra = self.last_cell_spectra[~mask_clone].reset_index(drop=True)

        if verbose:
            n_remaining = len(self.last_cell_spectra)
            remaining_clones = self.last_cell_spectra['Class'].nunique() if n_remaining > 0 else 0
            print(f"[INFO] Removed clone {clone_id}: {len(ids_to_remove)} cells removed.")
            print(f"[INFO] {n_remaining} cells remaining in {remaining_clones} clone(s).")

        return self.last_cell_spectra

    def combine_clones(self, clone_ids, new_label=None,
                       show_plotly=True, verbose=True,
                       color='cyan', contour_width=1.5):
        """
        Merge two or more molecular subpopulations (MSps / clones) into one.

        Useful when two MSps are biologically equivalent or when a small
        sub-cluster should be absorbed into a neighbour.

        Parameters
        ----------
        clone_ids : list of int — IDs of the MSps to merge (≥ 2 entries).
        new_label : int or None — Label for the merged MSp (default: min id).
        """
        if self.last_cell_spectra is None or 'Class' not in self.last_cell_spectra.columns:
            raise RuntimeError("Run clustering first — 'Class' column missing.")
        if len(clone_ids) < 2:
            raise ValueError("combine_clones() requires at least 2 clone IDs.")

        existing = set(self.last_cell_spectra['Class'].unique().astype(int))
        invalid = [c for c in clone_ids if int(c) not in existing]
        if invalid:
            raise ValueError(f"Clone IDs not found: {invalid}. Available: {sorted(existing)}")

        clone_ids = [int(c) for c in clone_ids]
        if new_label is None:
            new_label = min(clone_ids)
        new_label = int(new_label)

        df = self.last_cell_spectra.copy()
        mask_merge = df['Class'].astype(int).isin(clone_ids)
        df.loc[mask_merge, 'Class'] = new_label

        old_labels = sorted(df['Class'].unique().astype(int))
        label_map = {old: new for new, old in enumerate(old_labels)}
        df['Class'] = df['Class'].astype(int).map(label_map)
        new_label_mapped = label_map[new_label]

        prob_cols = [c for c in df.columns if c.startswith('prob_cluster_')]
        if prob_cols:
            merged_prob_cols = [f'prob_cluster_{c}' for c in clone_ids
                                if f'prob_cluster_{c}' in df.columns]
            if merged_prob_cols:
                df.loc[mask_merge, f'prob_cluster_{new_label_mapped}'] = \
                    df.loc[mask_merge, merged_prob_cols].sum(axis=1).clip(0, 1)
            keep_prob = [f'prob_cluster_{i}' for i in range(len(old_labels))]
            for pc_col in prob_cols:
                if pc_col not in keep_prob:
                    df.drop(columns=[pc_col], errors='ignore', inplace=True)

        self.last_cell_spectra = df

        for attr in ('clone_colors', 'clone_quality'):
            if hasattr(self, attr):
                delattr(self, attr)

        if verbose:
            n_merged = int(mask_merge.sum())
            n_new = int((df['Class'] == new_label_mapped).sum())
            remaining = sorted(df['Class'].unique().astype(int))
            print(f"[COMBINE] Merged MSps {clone_ids} → MSp {new_label_mapped} "
                  f"({n_merged} cells combined, new size={n_new}).")
            print(f"[COMBINE] Remaining MSps: {remaining}")

        return self.last_cell_spectra

    # ──────────────────────────────────────────────────────────────────────────
    #  SPECTRAL EXTRACTION
    # ──────────────────────────────────────────────────────────────────────────

    def extract_cell_spectra(self, agg='mean'):
        """
        Extract one representative spectrum per segmented cell.

        Parameters
        ----------
        agg : str — Aggregation: 'mean' (default), 'sum', or 'median'.

        Returns
        -------
        pd.DataFrame — Rows = cells, columns = metadata + m/z bins.
        """
        if self.cell_labels is None:
            raise RuntimeError("Cell segmentation missing.")
        if self.data_cube is None:
            raise RuntimeError("Data cube missing.")
        if self.mz_axis is None:
            raise RuntimeError("mz_axis missing.")

        labels = self.cell_labels
        meta = []
        spectra_matrix = []

        for lab in np.unique(labels):
            if lab == 0:
                continue
            mask = labels == lab
            pix_specs = self.data_cube[mask]
            if pix_specs.size == 0:
                continue
            spec = getattr(np, agg)(pix_specs, axis=0)
            props = measure.regionprops((labels == lab).astype(int))[0]
            centroid_row, centroid_col = props.centroid

            cell_meta = {
                'cell_id': int(lab),
                'area_px': int(props.area),
                'area_um': float(props.area * (self.pixel_size_um ** 2)),
                'centroid_row': float(centroid_row),
                'centroid_col': float(centroid_col),
            }

            if self.last_cell_spectra is not None:
                if 'Class' in self.last_cell_spectra.columns:
                    cell_meta['Class'] = self.last_cell_spectra.loc[
                        self.last_cell_spectra['cell_id'] == lab, 'Class'].values[0]
                if 'Score' in self.last_cell_spectra.columns:
                    cell_meta['Score'] = self.last_cell_spectra.loc[
                        self.last_cell_spectra['cell_id'] == lab, 'Score'].values[0]

            meta.append(cell_meta)
            spectra_matrix.append(spec)

        if len(spectra_matrix) == 0:
            df = pd.DataFrame(meta)
        else:
            df = pd.concat([
                pd.DataFrame(meta),
                pd.DataFrame(np.vstack(spectra_matrix),
                           columns=[f"{mz:.4f}" for mz in self.mz_axis])
            ], axis=1)

        self.last_cell_spectra = df
        if self.verbose:
            print(f"{len(df)} cell spectra extracted.")
        return df

    # ──────────────────────────────────────────────────────────────────────────
    #  QUALITY CONTROL
    # ──────────────────────────────────────────────────────────────────────────

    def compute_spectral_quality_metrics(self):
        """
        Compute QC metrics for each cell spectrum:
        SNR, sparsity, dynamic range, total intensity, spectral entropy, n_peaks.
        """
        if self.last_cell_spectra is None:
            raise RuntimeError("No cell spectra. Run extract_cell_spectra() first.")

        mz_cols = [c for c in self.last_cell_spectra.columns
                   if c not in ['cell_id', 'area_px', 'area_um', 'centroid_row',
                               'centroid_col', 'Class', 'Score', 'umap_0', 'umap_1']]

        X = self.last_cell_spectra[mz_cols].values

        qc_metrics = {}

        noise_level = np.percentile(X, 10, axis=1, keepdims=True)
        signal_level = np.mean(X, axis=1, keepdims=True)
        qc_metrics['SNR'] = (signal_level / (noise_level + 1e-9)).flatten()
        qc_metrics['sparsity'] = np.sum(X > 0, axis=1) / X.shape[1]
        qc_metrics['dynamic_range'] = np.log10(np.max(X, axis=1) / (np.percentile(X, 5, axis=1) + 1e-9))
        qc_metrics['total_intensity'] = np.sum(X, axis=1)

        X_norm = X / (X.sum(axis=1, keepdims=True) + 1e-9)
        qc_metrics['spectral_entropy'] = -np.sum(X_norm * np.log2(X_norm + 1e-9), axis=1)

        median_vals = np.median(X, axis=1, keepdims=True)
        qc_metrics['n_peaks'] = np.sum(X > median_vals, axis=1)

        for key, val in qc_metrics.items():
            self.last_cell_spectra[f'qc_{key}'] = val

        if self.verbose:
            print(f"Quality metrics computed for {len(self.last_cell_spectra)} cells.")
            print("\nQuality Metrics Summary:")
            for key in qc_metrics.keys():
                vals = qc_metrics[key]
                print(f"  {key}: mean={np.mean(vals):.3f}, std={np.std(vals):.3f}, "
                      f"min={np.min(vals):.3f}, max={np.max(vals):.3f}")

        self.qc_results = qc_metrics
        return qc_metrics

    def filter_low_quality_cells(self, snr_threshold=2.0, intensity_percentile=10,
                                  sparsity_min=0.05, entropy_min=None,
                                  show_plot=True, auto_remove=False):
        """
        Flag (and optionally remove) low-quality cells based on QC metrics.

        Parameters
        ----------
        snr_threshold : float       — Minimum SNR.
        intensity_percentile : float — Remove cells below this intensity percentile.
        sparsity_min : float        — Minimum fraction of non-zero bins.
        entropy_min : float|None    — Optional minimum spectral entropy.
        auto_remove : bool          — If True, removes flagged cells automatically.
        """
        if self.qc_results is None:
            print("[INFO] Computing quality metrics first...")
            self.compute_spectral_quality_metrics()

        df = self.last_cell_spectra.copy()

        flags = {}
        flags['low_snr'] = df['qc_SNR'] < snr_threshold

        intensity_thresh = np.percentile(df['qc_total_intensity'], intensity_percentile)
        flags['low_intensity'] = df['qc_total_intensity'] < intensity_thresh
        flags['low_sparsity'] = df['qc_sparsity'] < sparsity_min

        if entropy_min is not None:
            flags['low_entropy'] = df['qc_spectral_entropy'] < entropy_min

        combined_flag = np.zeros(len(df), dtype=bool)
        for flag_vals in flags.values():
            combined_flag |= flag_vals

        flagged_ids = df.loc[combined_flag, 'cell_id'].values

        if self.verbose:
            print(f"\n[QC FILTER RESULTS]")
            print(f"Total cells: {len(df)}")
            for flag_name, flag_vals in flags.items():
                n_flagged = np.sum(flag_vals)
                pct = 100 * n_flagged / len(df)
                print(f"  {flag_name}: {n_flagged} cells ({pct:.1f}%)")
            print(f"  TOTAL FLAGGED: {len(flagged_ids)} cells ({100*len(flagged_ids)/len(df):.1f}%)")

        if show_plot:
            import plotly.subplots as sp
            fig = sp.make_subplots(
                rows=2, cols=3,
                subplot_titles=('SNR Distribution', 'Total Intensity', 'Sparsity',
                               'Dynamic Range', 'Spectral Entropy', 'Number of Peaks')
            )
            metrics = ['qc_SNR', 'qc_total_intensity', 'qc_sparsity',
                      'qc_dynamic_range', 'qc_spectral_entropy', 'qc_n_peaks']
            for idx, metric in enumerate(metrics):
                row = idx // 3 + 1
                col = idx % 3 + 1
                fig.add_trace(go.Histogram(
                    x=df.loc[~combined_flag, metric], name='Good',
                    marker_color='green', opacity=0.6, showlegend=(idx == 0)
                ), row=row, col=col)
                fig.add_trace(go.Histogram(
                    x=df.loc[combined_flag, metric], name='Flagged',
                    marker_color='red', opacity=0.6, showlegend=(idx == 0)
                ), row=row, col=col)
                fig.update_xaxes(title_text=metric.replace('qc_', ''), row=row, col=col)
                fig.update_yaxes(title_text='Count', row=row, col=col)
            fig.update_layout(
                height=700, title_text="Quality Control Metrics - Good vs Flagged Cells",
                template="plotly_dark", barmode='overlay'
            )
            fig.show()

        if auto_remove and len(flagged_ids) > 0:
            print(f"\n[INFO] Auto-removing {len(flagged_ids)} low-quality cells...")
            self.remove_cells(flagged_ids.tolist(), show_plotly=False, verbose=False)

        return flagged_ids.tolist()

    # ──────────────────────────────────────────────────────────────────────────
    #  SOFT CLUSTERING (GMM)
    # ──────────────────────────────────────────────────────────────────────────

    def compute_soft_clustering(self, n_clusters=None, use_umap_space=True,
                                 n_components=2, n_neighbors=15, min_dist=0.1,
                                 metric='euclidean', use_scaler=True, random_state=0,
                                 covariance_type='tied', reg_covar=1e-6,
                                 max_clusters=10, min_clusters=2):
        """
        Soft clustering via Gaussian Mixture Models (GMM).

        Returns probabilistic membership (soft assignment) for each cell.

        Parameters
        ----------
        n_clusters : int|None — Number of clusters (None = auto via BIC/silhouette).
        covariance_type : str — 'tied' (recommended), 'full', 'diag', or 'spherical'.
        reg_covar : float     — Covariance regularisation (increase if ill-defined).
        max_clusters : int    — Upper bound for auto cluster selection.

        Returns
        -------
        soft_memberships : np.ndarray  shape (n_cells, n_clusters)
        hard_labels : np.ndarray       shape (n_cells,)
        confidence_scores : np.ndarray shape (n_cells,)
        """
        from sklearn.mixture import GaussianMixture
        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA
        from sklearn.metrics import silhouette_score

        if self.last_cell_spectra is None:
            raise RuntimeError("No cell spectra. Run extract_cell_spectra() first.")

        df = self.last_cell_spectra.copy()
        mz_cols = [c for c in df.columns if c not in [
            'cell_id', 'area_px', 'area_um', 'centroid_row', 'centroid_col',
            'Class', 'Score', 'umap_0', 'umap_1'
        ] and not c.startswith('qc_')]

        if len(mz_cols) == 0:
            raise RuntimeError("No m/z columns found.")

        X = df[mz_cols].values.astype(np.float64)
        if use_scaler:
            X = StandardScaler().fit_transform(X)

        try:
            import umap.umap_ as umap
            reducer = umap.UMAP(
                n_components=n_components, n_neighbors=n_neighbors,
                min_dist=min_dist, metric=metric, random_state=random_state
            )
            emb = reducer.fit_transform(X)
            if self.verbose:
                print("[INFO] UMAP reduction completed.")
        except Exception as e:
            if self.verbose:
                print(f"[WARN] UMAP failed: {e}")
            pca = PCA(n_components=min(10, X.shape[1]), random_state=random_state)
            emb = pca.fit_transform(X)[:, :n_components]
            if self.verbose:
                print("[INFO] Fallback to PCA reduction.")

        emb = emb.astype(np.float64)
        emb = (emb - emb.mean(axis=0)) / (emb.std(axis=0) + 1e-9)

        if n_clusters is None:
            bic_scores = []
            silhouette_scores = []
            candidates = range(min_clusters, min(max_clusters + 1, len(df) // 3))

            for k in candidates:
                try:
                    gmm_test = GaussianMixture(
                        n_components=k, covariance_type=covariance_type,
                        reg_covar=reg_covar, random_state=random_state,
                        max_iter=200, n_init=10
                    )
                    gmm_test.fit(emb)
                    bic_scores.append(gmm_test.bic(emb))
                    labels_test = gmm_test.predict(emb)
                    if len(np.unique(labels_test)) > 1:
                        silhouette_scores.append(silhouette_score(emb, labels_test))
                    else:
                        silhouette_scores.append(-1)
                except Exception as e:
                    if self.verbose:
                        print(f"[WARN] GMM failed for k={k}: {e}")
                    bic_scores.append(np.inf)
                    silhouette_scores.append(-1)

            valid_bic = [(i, s) for i, s in enumerate(bic_scores) if s != np.inf]
            if not valid_bic:
                raise RuntimeError("GMM failed for all k values.")

            best_idx_bic = min(valid_bic, key=lambda x: x[1])[0]
            best_k_bic = list(candidates)[best_idx_bic]

            valid_sil = [(i, s) for i, s in enumerate(silhouette_scores) if s > 0]
            if valid_sil:
                best_idx_sil = max(valid_sil, key=lambda x: x[1])[0]
                best_k_sil = list(candidates)[best_idx_sil]
                max_sil = silhouette_scores[best_idx_sil]
            else:
                best_k_sil = best_k_bic
                max_sil = -1

            n_clusters = best_k_sil if max_sil > 0.25 else best_k_bic

            if self.verbose:
                print(f"[INFO] Optimal clusters: {n_clusters}")

        try:
            gmm = GaussianMixture(
                n_components=n_clusters, covariance_type=covariance_type,
                reg_covar=reg_covar, random_state=random_state,
                max_iter=200, n_init=10
            )
            gmm.fit(emb)
        except Exception:
            if self.verbose:
                print("[INFO] Retrying with 'tied' covariance and higher regularization...")
            gmm = GaussianMixture(
                n_components=n_clusters, covariance_type='tied',
                reg_covar=1e-4, random_state=random_state,
                max_iter=200, n_init=10
            )
            gmm.fit(emb)

        soft_memberships = gmm.predict_proba(emb)
        hard_labels = gmm.predict(emb)
        confidence_scores = soft_memberships.max(axis=1)

        self.soft_memberships = soft_memberships

        df['umap_0'] = emb[:, 0]
        df['umap_1'] = emb[:, 1] if emb.shape[1] > 1 else 0.0
        df['Class'] = hard_labels.astype(int)
        df['soft_confidence'] = confidence_scores

        for k in range(n_clusters):
            df[f'prob_cluster_{k}'] = soft_memberships[:, k]

        self.last_cell_spectra = df

        if self.verbose:
            print(f"[INFO] Soft clustering completed: {n_clusters} clusters")
            print(f"  Mean confidence: {np.mean(confidence_scores):.3f}")
            print(f"  Cells with confidence < 0.5: {np.sum(confidence_scores < 0.5)}")

        return soft_memberships, hard_labels, confidence_scores

    def compute_fuzzy_cmeans(self, max_clusters=10, min_clusters=2, m=2.0,
                              max_iter=300, tol=1e-6, use_umap_space=True,
                              n_neighbors=15, min_dist=0.1, n_components=2,
                              random_state=0):
        """
        Soft clustering using Fuzzy C-Means (FCM).

        Falls back to a pure-numpy implementation if skfuzzy is not installed.

        Parameters
        ----------
        m : float — Fuzziness exponent (m=2 standard; m>3 very soft).
        """
        from sklearn.preprocessing import StandardScaler
        from sklearn.metrics import silhouette_score

        if self.last_cell_spectra is None:
            raise RuntimeError("No cell spectra. Run extract_cell_spectra() first.")

        df = self.last_cell_spectra.copy()
        meta_cols = {'cell_id','area_px','area_um','centroid_row','centroid_col',
                     'Class','Score','umap_0','umap_1','soft_confidence'}
        mz_cols = [c for c in df.columns
                   if c not in meta_cols
                   and not c.startswith('qc_') and not c.startswith('prob_cluster_')]

        X = df[mz_cols].values.astype(np.float64)
        X = StandardScaler().fit_transform(X)

        emb = X
        try:
            import umap.umap_ as umap_lib
            reducer = umap_lib.UMAP(n_components=n_components, n_neighbors=n_neighbors,
                                     min_dist=min_dist, random_state=random_state,
                                     metric='euclidean')
            emb = reducer.fit_transform(X)
            if self.verbose:
                print(f"[FCM] UMAP: {X.shape} → {emb.shape}")
        except Exception as e:
            if self.verbose:
                print(f"[FCM] UMAP unavailable ({e}), using scaled features.")

        def _numpy_fcm(data, c, m, max_iter, tol):
            n = data.shape[0]
            rng = np.random.default_rng(random_state)
            U = rng.random((c, n))
            U = U / U.sum(axis=0)
            for _ in range(max_iter):
                um = U ** m
                centers = (um @ data) / um.sum(axis=1, keepdims=True)
                dists = np.array([[np.linalg.norm(data[j] - centers[i]) ** 2
                                   for j in range(n)] for i in range(c)])
                dists = np.maximum(dists, 1e-10)
                exp = 2.0 / (m - 1)
                new_U = np.zeros_like(U)
                for i in range(c):
                    new_U[i] = 1.0 / np.sum((dists[i] / dists) ** exp, axis=0)
                delta = np.max(np.abs(new_U - U))
                U = new_U
                if delta < tol:
                    break
            return U

        candidates = list(range(min_clusters, max_clusters + 1))
        best_k, best_score = min_clusters, -1.0
        for c in candidates:
            try:
                try:
                    import skfuzzy
                    _, u, *_ = skfuzzy.cluster.cmeans(
                        emb.T, c=c, m=m, error=tol, maxiter=max_iter, seed=random_state)
                except ImportError:
                    u = _numpy_fcm(emb, c, m, max_iter, tol)
                hard = np.argmax(u, axis=0)
                if len(np.unique(hard)) < 2:
                    continue
                sc = silhouette_score(emb, hard, sample_size=min(1000, len(hard)))
                if sc > best_score:
                    best_score, best_k = sc, c
            except Exception:
                pass

        if self.verbose:
            print(f"[FCM] Optimal k={best_k}  (silhouette={best_score:.3f})")

        try:
            import skfuzzy
            _, u_final, *_ = skfuzzy.cluster.cmeans(
                emb.T, c=best_k, m=m, error=tol, maxiter=max_iter, seed=random_state)
        except ImportError:
            u_final = _numpy_fcm(emb, best_k, m, max_iter, tol)

        u_final = u_final.T
        hard_labels = np.argmax(u_final, axis=1)
        confidence  = u_final.max(axis=1)

        df['umap_0'] = emb[:, 0]
        df['umap_1'] = emb[:, 1] if emb.shape[1] > 1 else 0.0
        df['Class']  = hard_labels.astype(int)
        df['soft_confidence'] = confidence
        for k in range(best_k):
            df[f'prob_cluster_{k}'] = u_final[:, k]

        self.soft_memberships = u_final
        self.last_cell_spectra = df

        if self.verbose:
            print(f"[FCM] Done: {best_k} MSps, mean confidence={np.mean(confidence):.3f}")
        return u_final, hard_labels, confidence

    # ──────────────────────────────────────────────────────────────────────────
    #  CLONE QUALITY ASSESSMENT
    # ──────────────────────────────────────────────────────────────────────────

    def assess_clone_quality(self, min_confidence=0.6, min_cells_per_clone=5,
                            intensity_threshold_percentile=20, show_report=True):
        """
        Assess clone quality based on:
        - Soft clustering confidence scores
        - Within-clone spectral homogeneity (CV)
        - Between-clone separability (UMAP distance)
        - Clone size
        - Spectral intensity

        Returns dict with status: 'GOOD', 'QUESTIONABLE', or 'ARTIFACT'.
        """
        if self.soft_memberships is None:
            raise RuntimeError("Run compute_soft_clustering() first.")
        if 'Class' not in self.last_cell_spectra.columns:
            raise RuntimeError("No clone assignments found.")

        df = self.last_cell_spectra.copy()
        clones = sorted(df['Class'].unique())

        clone_quality = {}

        mz_cols = [c for c in df.columns if c not in [
            'cell_id', 'area_px', 'area_um', 'centroid_row', 'centroid_col',
            'Class', 'Score', 'umap_0', 'umap_1', 'soft_confidence'
        ] and not c.startswith('qc_') and not c.startswith('prob_cluster_')]

        X = df[mz_cols].values

        for cl in clones:
            mask = df['Class'] == cl
            n_cells = np.sum(mask)

            if n_cells < min_cells_per_clone:
                clone_quality[cl] = {
                    'status': 'LOW_SIZE', 'n_cells': n_cells,
                    'reason': f'Only {n_cells} cells (min={min_cells_per_clone})'
                }
                continue

            confidences = df.loc[mask, 'soft_confidence'].values
            mean_conf = np.mean(confidences)
            low_conf_ratio = np.sum(confidences < min_confidence) / n_cells

            clone_spectra = X[mask]
            cv_per_mz = np.std(clone_spectra, axis=0) / (np.mean(clone_spectra, axis=0) + 1e-9)
            mean_cv = np.mean(cv_per_mz[np.isfinite(cv_per_mz)])

            mean_intensity = np.mean(df.loc[mask, 'qc_total_intensity'])
            intensity_threshold = np.percentile(df['qc_total_intensity'], intensity_threshold_percentile)

            if 'umap_0' in df.columns:
                clone_center = df.loc[mask, ['umap_0', 'umap_1']].mean().values
                other_mask = df['Class'] != cl
                other_centers = df.loc[other_mask, ['umap_0', 'umap_1']].values
                if len(other_centers) > 0:
                    dists = np.linalg.norm(other_centers - clone_center, axis=1)
                    min_dist_to_other = np.min(dists)
                else:
                    min_dist_to_other = np.inf
            else:
                min_dist_to_other = None

            flags = []
            if mean_conf < 0.5:
                flags.append('LOW_CONFIDENCE')
            if low_conf_ratio > 0.4:
                flags.append('HIGH_UNCERTAINTY')
            if mean_cv > 1.5:
                flags.append('HIGH_VARIABILITY')
            if mean_intensity < intensity_threshold:
                flags.append('LOW_INTENSITY')
            if min_dist_to_other is not None and min_dist_to_other < 1.0:
                flags.append('POOR_SEPARATION')

            if len(flags) == 0:
                status = 'GOOD'
            elif len(flags) >= 3:
                status = 'ARTIFACT'
            else:
                status = 'QUESTIONABLE'

            clone_quality[cl] = {
                'status': status, 'n_cells': n_cells,
                'mean_confidence': mean_conf, 'low_conf_ratio': low_conf_ratio,
                'mean_cv': mean_cv, 'mean_intensity': mean_intensity,
                'min_dist_to_other': min_dist_to_other, 'flags': flags
            }

        self.clone_quality = clone_quality

        if show_report:
            print("\n" + "="*70)
            print("CLONE QUALITY ASSESSMENT REPORT")
            print("="*70)
            for cl in clones:
                info = clone_quality[cl]
                print(f"\n Clone {cl}: {info['status']}")
                print(f"   Cells: {info['n_cells']}")
                if 'mean_confidence' in info:
                    print(f"   Mean confidence: {info['mean_confidence']:.3f}")
                    print(f"   Low confidence ratio: {info['low_conf_ratio']:.3f}")
                    print(f"   Mean CV (variability): {info['mean_cv']:.3f}")
                    print(f"   Mean intensity: {info['mean_intensity']:.2e}")
                    if info['min_dist_to_other'] is not None:
                        print(f"   Min distance to other clones: {info['min_dist_to_other']:.3f}")
                    if info['flags']:
                        print(f"   ⚠️  Flags: {', '.join(info['flags'])}")
                else:
                    print(f"   ⚠️  {info['reason']}")
            print("\n" + "="*70)
            n_good = sum(1 for v in clone_quality.values() if v['status'] == 'GOOD')
            n_quest = sum(1 for v in clone_quality.values() if v['status'] == 'QUESTIONABLE')
            n_bad = sum(1 for v in clone_quality.values() if v['status'] in ['ARTIFACT', 'LOW_SIZE'])
            print(f"SUMMARY: {n_good} GOOD | {n_quest} QUESTIONABLE | {n_bad} ARTIFACT/LOW_SIZE")
            print("="*70 + "\n")

        return clone_quality

    # ──────────────────────────────────────────────────────────────────────────
    #  VISUALISATIONS
    # ──────────────────────────────────────────────────────────────────────────

    def show_umap_with_quality_overlay(self, figsize=(900, 700)):
        """UMAP scatter with marker opacity = soft confidence score."""
        if self.last_cell_spectra is None or 'umap_0' not in self.last_cell_spectra.columns:
            raise RuntimeError("UMAP not computed.")

        df = self.last_cell_spectra.copy()
        clones = sorted(df['Class'].unique())
        palette = pc.qualitative.Plotly
        clone_colors = {int(cl): palette[i % len(palette)] for i, cl in enumerate(clones)}
        self.clone_colors = clone_colors

        fig = go.Figure()
        for cl in clones:
            mask = df['Class'] == cl
            subset = df[mask]
            marker_sizes = 6 + 8 * subset['soft_confidence'].values
            fig.add_trace(go.Scatter(
                x=subset['umap_0'], y=subset['umap_1'],
                mode='markers',
                name=f'Clone {int(cl)}',
                marker=dict(
                    size=marker_sizes, color=clone_colors[int(cl)],
                    opacity=subset['soft_confidence'].values,
                    line=dict(width=0.5, color='white')
                ),
                text=[f"Cell {int(cid)}<br>Confidence: {conf:.2f}"
                      for cid, conf in zip(subset['cell_id'], subset['soft_confidence'])],
                hoverinfo='text'
            ))

        fig.update_layout(
            title="UMAP with Soft Clustering Confidence (opacity = confidence)",
            width=figsize[0], height=figsize[1],
            template="plotly_dark", xaxis_title="UMAP 1", yaxis_title="UMAP 2",
            legend_title="Clones"
        )
        fig.show()
        return fig

    def overlay_clusters_on_image_plotly(self, crop=None, TIC_cmap='gray',
                                         figsize=(900, 900), show_colorbar=False,
                                         contour_width=1.0, point_markers=True,
                                         marker_size=6, annotation_text=False,
                                         opacity_by_confidence=True):
        """Spatial TIC map with clone-coloured cell contours. Opacity = confidence."""
        if self.TIC is None or self.cell_labels is None:
            raise RuntimeError("TIC or cell_labels missing.")
        if self.last_cell_spectra is None or 'Class' not in self.last_cell_spectra.columns:
            raise RuntimeError("Clustering not performed.")

        tic = self.TIC.copy()
        labels = self.cell_labels.copy()
        if crop is not None:
            tic = tic[:crop[0], :crop[1]]
            labels = labels[:crop[0], :crop[1]]

        id_to_clone = dict(zip(
            self.last_cell_spectra['cell_id'].astype(int),
            self.last_cell_spectra['Class'].astype(int)
        ))
        id_to_confidence = {}
        if opacity_by_confidence and 'soft_confidence' in self.last_cell_spectra.columns:
            id_to_confidence = dict(zip(
                self.last_cell_spectra['cell_id'].astype(int),
                self.last_cell_spectra['soft_confidence']
            ))

        clones = sorted(self.last_cell_spectra['Class'].unique())
        if hasattr(self, "clone_colors"):
            clone_colors = self.clone_colors
        else:
            palette = pc.qualitative.Plotly
            clone_colors = {int(cl): palette[i % len(palette)] for i, cl in enumerate(clones)}

        fig = go.Figure()
        fig.add_trace(go.Heatmap(
            z=tic, colorscale=TIC_cmap, showscale=show_colorbar,
            hovertemplate="x=%{x}, y=%{y}<br>Intensity=%{z:.2e}<extra></extra>"
        ))

        for cid in np.unique(labels):
            if cid == 0 or cid not in id_to_clone:
                continue
            cl = id_to_clone[cid]
            color = clone_colors[int(cl)]
            mask = (labels == cid)
            opacity = 0.3 + 0.6 * id_to_confidence.get(cid, 1.0) if opacity_by_confidence else 0.9
            contours = measure.find_contours(mask.astype(float), 0.5)
            for contour in contours:
                conf_text = f" (conf={id_to_confidence.get(cid, 0):.2f})" if id_to_confidence else ""
                fig.add_trace(go.Scatter(
                    x=contour[:, 1], y=contour[:, 0], mode='lines',
                    line=dict(color=color, width=contour_width), opacity=opacity,
                    hoverinfo='text',
                    text=[f'Cell {int(cid)} (Clone {int(cl)}){conf_text}'] * len(contour),
                    showlegend=False
                ))

        for cl in clones:
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode='lines',
                line=dict(color=clone_colors[int(cl)], width=contour_width),
                name=f'Clone {int(cl)}'
            ))

        if point_markers:
            rows = self.last_cell_spectra['centroid_row'].values
            cols = self.last_cell_spectra['centroid_col'].values
            labels_cl = self.last_cell_spectra['Class'].values
            marker_colors = [clone_colors[int(c)] for c in labels_cl]
            fig.add_trace(go.Scatter(
                x=cols, y=rows, mode='markers',
                marker=dict(size=marker_size, color=marker_colors,
                           line=dict(width=1, color='black')),
                hoverinfo='text', showlegend=False
            ))

        fig.update_layout(
            title="TIC + Clone Contours (opacity = confidence)",
            width=figsize[0], height=figsize[1],
            yaxis=dict(scaleanchor="x", autorange='reversed'),
            xaxis=dict(showgrid=False, zeroline=False, visible=False),
            yaxis_visible=False, template="simple_white",
            margin=dict(l=0, r=0, t=40, b=0), legend_title="Clones"
        )
        fig.show()
        return fig

    def plot_cell_contours_by_clone(self, crop=None, contour_width=2.0,
                                     figsize=(900, 900),
                                     contour_color_override=None, show_tic=False):
        """Cell contours coloured by clone on white (or TIC) background."""
        if self.cell_labels is None:
            raise RuntimeError("cell_labels missing.")
        if self.last_cell_spectra is None or 'Class' not in self.last_cell_spectra.columns:
            raise RuntimeError("'Class' missing in last_cell_spectra.")

        labels = self.cell_labels.copy()
        if crop is not None:
            labels = labels[:crop[0], :crop[1]]

        id_to_clone = dict(zip(
            self.last_cell_spectra['cell_id'].astype(int),
            self.last_cell_spectra['Class'].astype(int)
        ))
        clones = sorted(self.last_cell_spectra['Class'].unique())
        palette = pc.qualitative.Plotly
        clone_colors = {int(cl): palette[i % len(palette)] for i, cl in enumerate(clones)}

        fig = go.Figure()
        if show_tic and self.TIC is not None:
            fig.add_trace(go.Heatmap(
                z=self.TIC, colorscale='gray', showscale=False,
                hovertemplate="x=%{x}, y=%{y}<br>TIC=%{z:.2e}<extra></extra>"
            ))
            bg_color = 'black'
        else:
            fig.add_shape(
                type="rect", x0=0, y0=0,
                x1=labels.shape[1], y1=labels.shape[0],
                line=dict(width=0), fillcolor="white", layer="below"
            )
            bg_color = 'white'

        for cid in np.unique(labels):
            if cid == 0 or cid not in id_to_clone:
                continue
            cl = id_to_clone[cid]
            color = contour_color_override if contour_color_override else clone_colors[int(cl)]
            mask = (labels == cid)
            for contour in measure.find_contours(mask.astype(float), 0.5):
                fig.add_trace(go.Scatter(
                    x=contour[:, 1], y=contour[:, 0], mode='lines',
                    line=dict(color=color, width=contour_width), opacity=1.0,
                    hoverinfo='text',
                    text=[f'Cell {int(cid)} (Clone {int(cl)})'] * len(contour),
                    name=f'Clone {int(cl)}', showlegend=False
                ))

        for cl in clones:
            leg_color = contour_color_override if contour_color_override else clone_colors[int(cl)]
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode='lines',
                line=dict(color=leg_color, width=contour_width), name=f'Clone {int(cl)}'
            ))

        fig.update_layout(
            title="Cell Contours Colored by Clone" + (" — TIC background" if show_tic else ""),
            width=figsize[0], height=figsize[1],
            yaxis=dict(scaleanchor="x", autorange='reversed'),
            xaxis=dict(showgrid=False, zeroline=False, visible=False),
            yaxis_visible=False, plot_bgcolor=bg_color, paper_bgcolor=bg_color,
            margin=dict(l=0, r=0, t=40, b=0), legend_title="Clones"
        )
        fig.show()
        return fig

    # ──────────────────────────────────────────────────────────────────────────
    #  SPECTRAL PLOTS
    # ──────────────────────────────────────────────────────────────────────────

    def global_cell_spectrum(self, agg='mean'):
        """Return aggregated spectrum across all cells."""
        if self.last_cell_spectra is None:
            raise RuntimeError("No cell spectra extracted.")
        meta_cols = {'cell_id','area_px','area_um','centroid_row','centroid_col',
                     'Class','Score','umap_0','umap_1','soft_confidence'}
        mz_cols = [c for c in self.last_cell_spectra.columns
                   if c not in meta_cols
                   and not c.startswith('qc_') and not c.startswith('prob_cluster_')]
        X = self.last_cell_spectra[mz_cols].values.astype(float)
        if agg == 'mean':   return np.mean(X, axis=0)
        if agg == 'median': return np.median(X, axis=0)
        if agg == 'max':    return np.max(X, axis=0)
        if agg == 'min':    return np.min(X, axis=0)
        return np.mean(X, axis=0)

    def plot_global_spectrum(self, agg='mean', show_std=True, normalize=False):
        """Plot the global mean/median spectrum across all cells with ±SD band."""
        if self.last_cell_spectra is None:
            raise RuntimeError("No cell spectra extracted.")
        meta_cols = {'cell_id','area_px','area_um','centroid_row','centroid_col',
                     'Class','Score','umap_0','umap_1','soft_confidence'}
        mz_cols = [c for c in self.last_cell_spectra.columns
                   if c not in meta_cols
                   and not c.startswith('qc_') and not c.startswith('prob_cluster_')]
        X = self.last_cell_spectra[mz_cols].values.astype(float)

        if agg == 'mean':   spectrum = np.mean(X, axis=0)
        elif agg == 'median': spectrum = np.median(X, axis=0)
        elif agg == 'max':  spectrum = np.max(X, axis=0)
        elif agg == 'min':  spectrum = np.min(X, axis=0)
        else:               spectrum = np.mean(X, axis=0)

        std_spec = np.std(X, axis=0) if agg in ('mean', 'median') else None

        if normalize and spectrum.max() > 0:
            if std_spec is not None:
                std_spec = std_spec / spectrum.max()
            spectrum = spectrum / spectrum.max()

        color = '#89b4fa'
        fig = go.Figure()

        if show_std and std_spec is not None:
            fig.add_trace(go.Scatter(
                x=np.concatenate([self.mz_axis, self.mz_axis[::-1]]),
                y=np.concatenate([spectrum + std_spec, (spectrum - std_spec)[::-1]]),
                fill='toself', fillcolor='rgba(137,180,250,0.15)',
                line=dict(width=0), showlegend=False, hoverinfo='skip',
                name='±1 SD'
            ))

        fig.add_trace(go.Scatter(
            x=self.mz_axis, y=spectrum, mode='lines',
            line=dict(width=1.5, color=color),
            name=f'{agg.capitalize()} spectrum  (n={len(X)})',
            hovertemplate='m/z=%{x:.4f}<br>Intensity=%{y:.2e}<extra></extra>'
        ))

        y_label = f"{'Normalized i' if normalize else 'I'}ntensity ({agg})"
        fig.update_layout(
            title=f"Global Cell Spectrum — {agg}  (n={len(X)} cells)",
            xaxis_title="m/z", yaxis_title=y_label,
            template="plotly_dark", height=450, hovermode="x unified",
            legend=dict(orientation='h', y=-0.15)
        )
        fig.show()
        return spectrum

    def plot_clone_spectra_with_uncertainty(self, agg='mean', show_std=True,
                                             normalize=False, mode='overlay'):
        """
        Plot MSp spectra with ±SD bands.

        Parameters
        ----------
        mode : str — 'overlay' (all on one figure) or 'individual' (subplots).
        """
        if self.last_cell_spectra is None:
            raise RuntimeError("No cell spectra extracted.")
        if 'Class' not in self.last_cell_spectra.columns:
            raise RuntimeError("Clustering not performed.")

        meta_cols = {'cell_id','area_px','area_um','centroid_row','centroid_col',
                     'Class','Score','umap_0','umap_1','soft_confidence'}
        mz_cols = [c for c in self.last_cell_spectra.columns
                   if c not in meta_cols
                   and not c.startswith('qc_') and not c.startswith('prob_cluster_')]

        clones = sorted(self.last_cell_spectra['Class'].unique())
        if hasattr(self, 'clone_colors'):
            clone_colors = self.clone_colors
        else:
            palette = pc.qualitative.Plotly
            clone_colors = {int(cl): palette[i % len(palette)] for i, cl in enumerate(clones)}

        def _compute(cl):
            data = self.last_cell_spectra.loc[
                self.last_cell_spectra['Class'] == cl, mz_cols].values.astype(float)
            if agg == 'mean':   s = np.mean(data, axis=0)
            elif agg == 'median': s = np.median(data, axis=0)
            elif agg == 'max':  s = np.max(data, axis=0)
            elif agg == 'min':  s = np.min(data, axis=0)
            else:               s = np.mean(data, axis=0)
            sd = np.std(data, axis=0) if agg in ('mean','median') else None
            if normalize and s.max() > 0:
                if sd is not None: sd = sd / s.max()
                s = s / s.max()
            return s, sd, len(data)

        y_label = f"{'Norm. i' if normalize else 'I'}ntensity ({agg})"

        if mode == 'overlay':
            fig = go.Figure()
            for cl in clones:
                s, sd, n = _compute(cl)
                color = clone_colors[int(cl)]
                hex_c = color.lstrip('#')
                r, g_, b = int(hex_c[0:2],16), int(hex_c[2:4],16), int(hex_c[4:6],16)
                if show_std and sd is not None:
                    fig.add_trace(go.Scatter(
                        x=np.concatenate([self.mz_axis, self.mz_axis[::-1]]),
                        y=np.concatenate([s+sd, (s-sd)[::-1]]),
                        fill='toself', fillcolor=f'rgba({r},{g_},{b},0.15)',
                        line=dict(width=0), showlegend=False, hoverinfo='skip'
                    ))
                fig.add_trace(go.Scatter(
                    x=self.mz_axis, y=s, mode='lines',
                    line=dict(width=2, color=color),
                    name=f"MSp {int(cl)}  (n={n})",
                    hovertemplate='m/z=%{x:.4f}<br>I=%{y:.2e}<extra></extra>'
                ))
            fig.update_layout(
                title=f"MSp Spectra — {agg} overlay{'  (normalised)' if normalize else ''}",
                xaxis_title="m/z", yaxis_title=y_label,
                template="plotly_dark", height=500,
                hovermode="x unified", legend_title="MSp"
            )
            fig.show()
            return fig

        from plotly.subplots import make_subplots
        n_clones = len(clones)
        fig = make_subplots(
            rows=n_clones, cols=1, shared_xaxes=True,
            vertical_spacing=0.04,
            subplot_titles=[f"MSp {int(cl)}" for cl in clones]
        )
        for i, cl in enumerate(clones, 1):
            s, sd, n = _compute(cl)
            color = clone_colors[int(cl)]
            hex_c = color.lstrip('#')
            r, g_, b = int(hex_c[0:2],16), int(hex_c[2:4],16), int(hex_c[4:6],16)
            if show_std and sd is not None:
                fig.add_trace(go.Scatter(
                    x=np.concatenate([self.mz_axis, self.mz_axis[::-1]]),
                    y=np.concatenate([s+sd, (s-sd)[::-1]]),
                    fill='toself', fillcolor=f'rgba({r},{g_},{b},0.15)',
                    line=dict(width=0), showlegend=False, hoverinfo='skip'
                ), row=i, col=1)
            fig.add_trace(go.Scatter(
                x=self.mz_axis, y=s, mode='lines',
                line=dict(width=1.8, color=color),
                name=f"MSp {int(cl)}  (n={n})",
                hovertemplate='m/z=%{x:.4f}<br>I=%{y:.2e}<extra></extra>'
            ), row=i, col=1)
            fig.update_yaxes(title_text=y_label, row=i, col=1)

        fig.update_xaxes(title_text="m/z", row=n_clones, col=1)
        fig.update_layout(
            title=f"MSp Spectra — {agg} individual{'  (normalised)' if normalize else ''}",
            template="plotly_dark", height=max(300, 220 * n_clones),
            hovermode="x unified", showlegend=True, legend_title="MSp"
        )
        fig.show()
        return fig

    # ──────────────────────────────────────────────────────────────────────────
    #  ION SPATIAL DISTRIBUTION
    # ──────────────────────────────────────────────────────────────────────────

    def show_mz_distribution_in_cells_plotly(self, mz_value, tolerance=0.1,
                                              cmap='magma', figsize=(700, 700),
                                              smoothing_sigma=1.0, clip_min=0.0,
                                              clip_max=None, show_contours=True,
                                              contour_color='cyan', contour_width=1.5,
                                              contour_opacity=0.9):
        """
        Display spatial distribution of a selected ion (m/z) within cells.

        Parameters
        ----------
        mz_value : float     — Central m/z value.
        tolerance : float    — ± integration window (Da).
        smoothing_sigma : float — Gaussian smoothing (0 = off).
        """
        if self.data_cube is None or self.cell_labels is None:
            raise RuntimeError("Data cube or cell labels missing.")

        mz_idx = np.where(
            (self.mz_axis >= mz_value - tolerance) & (self.mz_axis <= mz_value + tolerance)
        )[0]
        if len(mz_idx) == 0:
            raise ValueError(f"No bin found for m/z={mz_value:.2f} ± {tolerance}")

        pic_cube = self.data_cube[:, :, mz_idx].sum(axis=2)

        if self.last_cell_spectra is not None:
            valid_cell_ids = set(self.last_cell_spectra['cell_id'].astype(int))
            cell_mask = np.isin(self.cell_labels, list(valid_cell_ids))
        else:
            cell_mask = self.cell_labels > 0

        pic_masked = pic_cube * cell_mask

        if smoothing_sigma > 0:
            pic_masked = gaussian_filter(pic_masked, sigma=smoothing_sigma)

        if clip_max is not None:
            pic_masked = np.clip(pic_masked, clip_min, clip_max)

        fig = go.Figure()
        fig.add_trace(go.Heatmap(
            z=pic_masked, colorscale=cmap,
            colorbar=dict(title=f"Intensity m/z={mz_value:.2f}"),
            hovertemplate="x=%{x}, y=%{y}<br>Intensity=%{z:.2e}<extra></extra>"
        ))

        if show_contours:
            for lab in np.unique(self.cell_labels):
                if lab == 0:
                    continue
                if self.last_cell_spectra is not None and lab not in valid_cell_ids:
                    continue
                mask = self.cell_labels == lab
                for contour in measure.find_contours(mask, 0.5):
                    fig.add_trace(go.Scatter(
                        x=contour[:, 1], y=contour[:, 0], mode='lines',
                        line=dict(color=contour_color, width=contour_width),
                        opacity=contour_opacity, hoverinfo='text',
                        text=[f'Cell ID: {lab}'] * len(contour), showlegend=False
                    ))

        fig.update_layout(
            title=f"Ion Distribution m/z={mz_value:.2f} (σ={smoothing_sigma})",
            width=figsize[0], height=figsize[1],
            yaxis=dict(scaleanchor="x", autorange='reversed'),
            xaxis=dict(showgrid=False, zeroline=False),
            template="plotly_dark", paper_bgcolor="black", plot_bgcolor="black",
            margin=dict(l=0, r=0, t=40, b=0)
        )
        fig.show()
        return pic_masked

    def plot_multi_ion_overlay(self, mz_list, tolerance=0.1, colors=None,
                                smoothing_sigma=0.8, opacity=0.75,
                                show_contours=True, contour_color='white',
                                contour_width=0.8, figsize=(900, 900),
                                normalize_each=True):
        """
        Overlay multiple ions on the spatial map, each with a distinct colour.

        Parameters
        ----------
        mz_list : list of float — m/z values to display (max 6 recommended).
        colors : list of str    — HTML colours (auto-assigned if None).
        normalize_each : bool   — Normalise each ion map to [0,1].
        """
        if self.data_cube is None or self.cell_labels is None:
            raise RuntimeError("Data cube or cell labels missing.")

        default_colors = ['#ff4d6d','#00f5d4','#fee440','#7b2fff',
                          '#ff7c43','#00b4d8','#80b918','#f72585']
        if colors is None:
            colors = [default_colors[i % len(default_colors)] for i in range(len(mz_list))]
        elif len(colors) < len(mz_list):
            colors = list(colors) + [default_colors[i % len(default_colors)]
                                     for i in range(len(mz_list) - len(colors))]

        valid_cell_ids = set(self.last_cell_spectra['cell_id'].astype(int)) \
            if self.last_cell_spectra is not None else None

        ion_maps = {}
        fig = go.Figure()

        fig.add_shape(type="rect", x0=0, y0=0,
                      x1=self.cell_labels.shape[1], y1=self.cell_labels.shape[0],
                      line=dict(width=0), fillcolor="black", layer="below")

        for i, mz_val in enumerate(mz_list):
            idx = np.where(
                (self.mz_axis >= mz_val - tolerance) &
                (self.mz_axis <= mz_val + tolerance)
            )[0]
            if len(idx) == 0:
                if self.verbose:
                    print(f"[WARN] m/z={mz_val:.4f}: no bin found — skipped.")
                continue

            pic = self.data_cube[:, :, idx].sum(axis=2).astype(float)

            if valid_cell_ids is not None:
                cell_mask = np.isin(self.cell_labels, list(valid_cell_ids))
            else:
                cell_mask = self.cell_labels > 0
            pic *= cell_mask

            if smoothing_sigma > 0:
                pic = gaussian_filter(pic, sigma=smoothing_sigma)

            vmax = pic.max()
            if normalize_each and vmax > 0:
                pic = pic / vmax

            ion_maps[mz_val] = pic

            hex_c = colors[i].lstrip('#')
            r, g, b = int(hex_c[0:2], 16), int(hex_c[2:4], 16), int(hex_c[4:6], 16)
            colorscale_i = [
                [0.0, f'rgba({r},{g},{b},0)'],
                [1.0, f'rgba({r},{g},{b},1)']
            ]
            fig.add_trace(go.Heatmap(
                z=pic, colorscale=colorscale_i, showscale=False,
                opacity=opacity, zmin=0, zmax=1 if normalize_each else vmax,
                hovertemplate=f"m/z={mz_val:.4f}<br>x=%{{x}}, y=%{{y}}<br>I=%{{z:.3f}}<extra></extra>",
                name=f"m/z {mz_val:.4f}"
            ))

        if show_contours and self.cell_labels is not None:
            for lab in np.unique(self.cell_labels):
                if lab == 0:
                    continue
                if valid_cell_ids is not None and lab not in valid_cell_ids:
                    continue
                for contour in measure.find_contours((self.cell_labels == lab).astype(float), 0.5):
                    fig.add_trace(go.Scatter(
                        x=contour[:, 1], y=contour[:, 0], mode='lines',
                        line=dict(color=contour_color, width=contour_width),
                        opacity=0.6, hoverinfo='skip', showlegend=False
                    ))

        for i, mz_val in enumerate(mz_list):
            if mz_val not in ion_maps:
                continue
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode='markers',
                marker=dict(size=12, color=colors[i], symbol='square'),
                name=f"m/z {mz_val:.4f}"
            ))

        fig.update_layout(
            title=f"Multi-Ion Spatial Overlay ({len(ion_maps)} ions)",
            width=figsize[0], height=figsize[1],
            yaxis=dict(scaleanchor="x", autorange='reversed'),
            xaxis=dict(showgrid=False, zeroline=False, visible=False),
            yaxis_visible=False, paper_bgcolor='black', plot_bgcolor='black',
            margin=dict(l=0, r=0, t=45, b=0),
            legend=dict(bgcolor='rgba(20,20,30,0.85)', bordercolor='#444',
                        borderwidth=1, font=dict(color='white', size=11)),
            legend_title=dict(text="Ions", font=dict(color='white'))
        )
        fig.show()
        return ion_maps

    # ──────────────────────────────────────────────────────────────────────────
    #  DISTRIBUTION PLOTS
    # ──────────────────────────────────────────────────────────────────────────

    def plot_distribution_by_clone(self, mz_value=None, tolerance=0.1,
                                    metric=None, plot_type='violin',
                                    points='scatter', normalize=False,
                                    figsize=(900, 550)):
        """
        Violin / box / combined plot of an ion intensity or QC metric by MSp.

        Parameters
        ----------
        mz_value : float|None — Plot ion intensity ±tolerance.
        metric : str|None     — Column name or 'TIC' for TIC at centroid.
        plot_type : str       — 'violin', 'box', or 'both'.
        points : str          — 'scatter', 'outliers', or 'none'.
        """
        if self.last_cell_spectra is None or 'Class' not in self.last_cell_spectra.columns:
            raise RuntimeError("Run clustering first.")

        df = self.last_cell_spectra.copy()
        clones = sorted(df['Class'].astype(int).unique())
        palette = pc.qualitative.Plotly
        if hasattr(self, 'clone_colors'):
            clone_colors = self.clone_colors
        else:
            clone_colors = {int(cl): palette[i % len(palette)] for i, cl in enumerate(clones)}

        if mz_value is not None and mz_value > 0:
            meta_cols = {'cell_id','area_px','area_um','centroid_row','centroid_col',
                         'Class','Score','umap_0','umap_1','soft_confidence'}
            mz_cols_f = [c for c in df.columns
                         if c not in meta_cols
                         and not c.startswith('qc_') and not c.startswith('prob_cluster_')]
            mz_floats = np.array([float(c) for c in mz_cols_f])
            idx = np.where(np.abs(mz_floats - mz_value) <= tolerance)[0]
            if len(idx) == 0:
                raise ValueError(f"No bin found for m/z={mz_value:.4f} ±{tolerance}")
            values = df[np.array(mz_cols_f)[idx]].sum(axis=1).values
            ylabel = f"Intensity  m/z={mz_value:.4f} ±{tolerance}"
            title  = f"Ion distribution — m/z={mz_value:.4f}"
        elif metric == 'TIC' and self.TIC is not None:
            rows = np.clip(df['centroid_row'].values.astype(int), 0, self.TIC.shape[0]-1)
            cols = np.clip(df['centroid_col'].values.astype(int), 0, self.TIC.shape[1]-1)
            values = self.TIC[rows, cols]
            ylabel = "TIC at centroid"
            title  = "TIC distribution by MSp"
        elif metric and metric in df.columns:
            values = df[metric].values.astype(float)
            ylabel = metric
            title  = f"{metric} distribution by MSp"
        else:
            raise ValueError("Provide mz_value > 0 or a valid metric name.")

        if normalize:
            vmin, vmax = values.min(), values.max()
            if vmax > vmin:
                values = (values - vmin) / (vmax - vmin)
            ylabel = f"{ylabel} (norm.)"

        df['_val'] = values
        pt_map = {'scatter': 'all', 'outliers': 'outliers', 'none': False}
        pts_arg = pt_map.get(points, 'all')

        fig = go.Figure()
        for cl in clones:
            mask = df['Class'].astype(int) == cl
            y = df.loc[mask, '_val'].values
            color = clone_colors[int(cl)]
            hex_c = color.lstrip('#')
            r, g_, b_ = int(hex_c[0:2],16), int(hex_c[2:4],16), int(hex_c[4:6],16)
            fill_rgba = f'rgba({r},{g_},{b_},0.35)'
            line_rgba = f'rgba({r},{g_},{b_},1.0)'
            hover = [f"Cell {cid}<br>MSp {int(cl)}<br>val={v:.4f}"
                     for cid, v in zip(df.loc[mask, 'cell_id'].astype(int).values, y)]

            if plot_type in ('violin', 'both'):
                fig.add_trace(go.Violin(
                    y=y, name=f"MSp {int(cl)}",
                    box_visible=(plot_type == 'both'),
                    meanline_visible=True, points=pts_arg,
                    fillcolor=fill_rgba, line_color=line_rgba, opacity=0.85,
                    text=hover, hoverinfo='text',
                    legendgroup=f"cl{int(cl)}", showlegend=(plot_type != 'both')
                ))
            if plot_type == 'box':
                box_pts = 'all' if points == 'scatter' else pts_arg
                fig.add_trace(go.Box(
                    y=y, name=f"MSp {int(cl)}",
                    boxpoints=box_pts,
                    marker=dict(color=line_rgba, size=4, opacity=0.6),
                    line=dict(color=line_rgba), fillcolor=fill_rgba,
                    text=hover, hoverinfo='text',
                    legendgroup=f"cl{int(cl)}"
                ))

        fig.update_layout(
            title=title, yaxis_title=ylabel, xaxis_title="MSp",
            template="plotly_dark", width=figsize[0], height=figsize[1],
            violinmode='group', boxmode='group',
            legend_title="MSp", margin=dict(t=50, b=50)
        )
        fig.show()
        return fig

    # ──────────────────────────────────────────────────────────────────────────
    #  DIFFERENTIAL ANALYSIS
    # ──────────────────────────────────────────────────────────────────────────

    def run_differential_analysis(self, top_n=30, pval_threshold=0.05,
                                   fc_threshold=1.5, method='kruskal',
                                   fdr_correction=True,
                                   show_volcano=True, show_heatmap=True,
                                   figsize=(1000, 700)):
        """
        Differential analysis of m/z features across MSps.

        Parameters
        ----------
        top_n : int            — Top differential ions to return.
        pval_threshold : float — p-value threshold (raw or adjusted, see fdr_correction).
        fc_threshold : float   — Minimum fold-change.
        method : str           — 'kruskal' (non-parametric) or 'anova'.
        fdr_correction : bool  — If True (default), apply Benjamini-Hochberg FDR correction
                                 and filter on adjusted p-value.
                                 If False, use raw p-values (no multiple-testing correction).
        show_volcano : bool    — Show volcano plot.
        show_heatmap : bool    — Show Z-score heatmap.

        Returns
        -------
        top_markers : pd.DataFrame
        result_df   : pd.DataFrame (full results, always contains both pval and pval_adj)
        """
        from scipy import stats as scipy_stats

        if self.last_cell_spectra is None or 'Class' not in self.last_cell_spectra.columns:
            raise RuntimeError("Run clustering first.")

        df = self.last_cell_spectra
        clones = sorted(df['Class'].astype(int).unique())
        if len(clones) < 2:
            raise RuntimeError("Need at least 2 clones for differential analysis.")

        meta_cols = {'cell_id','area_px','area_um','centroid_row','centroid_col',
                     'Class','Score','umap_0','umap_1','soft_confidence'}
        mz_cols = [c for c in df.columns
                   if c not in meta_cols
                   and not c.startswith('qc_')
                   and not c.startswith('prob_cluster_')]

        X = df[mz_cols].values.astype(float)
        groups = [X[df['Class'].astype(int) == cl] for cl in clones]

        stat_vals, p_vals = [], []
        for j in range(len(mz_cols)):
            cols_j = [g[:, j] for g in groups if len(g) > 0]
            if len(cols_j) < 2:
                stat_vals.append(0.0); p_vals.append(1.0); continue
            try:
                if method == 'kruskal':
                    s, p = scipy_stats.kruskal(*cols_j)
                else:
                    s, p = scipy_stats.f_oneway(*cols_j)
            except Exception:
                s, p = 0.0, 1.0
            stat_vals.append(float(s)); p_vals.append(float(p))

        stat_arr = np.array(stat_vals)
        p_arr    = np.array(p_vals)

        # ── BH FDR correction (optional) ────────────────────────────────────
        if fdr_correction:
            n_tests = len(p_arr)
            order   = np.argsort(p_arr)
            rank    = np.empty_like(order); rank[order] = np.arange(1, n_tests + 1)
            p_adj   = np.minimum(1.0, p_arr * n_tests / rank)
            for i in range(n_tests - 2, -1, -1):
                p_adj[order[i]] = min(p_adj[order[i]], p_adj[order[i + 1]])
            p_used      = p_adj          # threshold applied on adjusted p
            p_used_col  = 'pval_adj'
            pval_label  = 'adj. p-value (BH FDR)'
            thresh_label = f"FDR={pval_threshold}"
        else:
            p_adj       = p_arr.copy()   # pval_adj column = raw p (no correction)
            p_used      = p_arr
            p_used_col  = 'pval'
            pval_label  = 'raw p-value'
            thresh_label = f"p={pval_threshold}"

        global_means = np.mean(X, axis=0) + 1e-12
        clone_means  = np.vstack([np.mean(groups[i], axis=0) for i in range(len(clones))])
        fc           = clone_means.max(axis=0) / global_means
        best_clone   = np.array([clones[i] for i in clone_means.argmax(axis=0)])

        result_df = pd.DataFrame({
            'mz':          [float(c) for c in mz_cols],
            'stat':        stat_arr,
            'pval':        p_arr,
            'pval_adj':    p_adj,
            'fold_change': fc,
            'best_clone':  best_clone,
            'neg_log10_p': -np.log10(np.clip(p_used, 1e-300, 1.0))
        })

        sig = result_df[(p_used < pval_threshold) &
                        (result_df['fold_change'] >= fc_threshold if fc_threshold > 0 else True)]
        sig = sig.sort_values('neg_log10_p', ascending=False)
        top_markers = sig.head(top_n).reset_index(drop=True)

        if self.verbose:
            corr_tag = "BH FDR" if fdr_correction else "raw p"
            print(f"[DA] {len(sig)} significant ions  "
                  f"({corr_tag}<{pval_threshold}, FC≥{fc_threshold})")
            for _, row in top_markers.head(5).iterrows():
                print(f"     m/z={row['mz']:.4f}  "
                      f"{'adj-p' if fdr_correction else 'p'}={row[p_used_col]:.2e}  "
                      f"FC={row['fold_change']:.2f}  best→Clone {int(row['best_clone'])}")

        palette = pc.qualitative.Plotly
        clone_colors = {int(cl): palette[i % len(palette)] for i, cl in enumerate(clones)}

        if show_volcano:
            is_sig = (p_used < pval_threshold) & \
                     (result_df['fold_change'] >= fc_threshold)
            fig_v = go.Figure()
            mask_ns = ~is_sig
            fig_v.add_trace(go.Scatter(
                x=np.log2(result_df.loc[mask_ns, 'fold_change'] + 1e-6),
                y=result_df.loc[mask_ns, 'neg_log10_p'],
                mode='markers', marker=dict(size=4, color='#585b70', opacity=0.5),
                name='n.s.',
                text=[f"m/z={m:.4f}" for m in result_df.loc[mask_ns, 'mz']],
                hoverinfo='text'
            ))
            for cl in clones:
                mask_cl = is_sig & (result_df['best_clone'].astype(int) == cl)
                if mask_cl.sum() == 0:
                    continue
                fig_v.add_trace(go.Scatter(
                    x=np.log2(result_df.loc[mask_cl, 'fold_change'] + 1e-6),
                    y=result_df.loc[mask_cl, 'neg_log10_p'],
                    mode='markers',
                    marker=dict(size=7, color=clone_colors[cl],
                                line=dict(width=0.5, color='white')),
                    name=f'Clone {cl}',
                    text=[f"m/z={m:.4f}<br>FC={fc_:.2f}<br>"
                          f"{'adj-p' if fdr_correction else 'p'}={p:.2e}"
                          for m, fc_, p in zip(result_df.loc[mask_cl, 'mz'],
                                              result_df.loc[mask_cl, 'fold_change'],
                                              result_df.loc[mask_cl, p_used_col])],
                    hoverinfo='text'
                ))
            fig_v.add_hline(y=-np.log10(pval_threshold),
                            line=dict(color='#f38ba8', dash='dash', width=1),
                            annotation_text=thresh_label)
            fig_v.add_vline(x=np.log2(fc_threshold + 1e-6),
                            line=dict(color='#fab387', dash='dash', width=1),
                            annotation_text=f"FC={fc_threshold}")
            for _, row in top_markers.head(10).iterrows():
                fig_v.add_annotation(
                    x=np.log2(row['fold_change'] + 1e-6),
                    y=row['neg_log10_p'],
                    text=f"  {row['mz']:.2f}", showarrow=False,
                    font=dict(size=9, color='white'), xanchor='left'
                )
            fig_v.update_layout(
                title=f"Volcano — {method.upper()}  |  {pval_label}  threshold={pval_threshold}",
                xaxis_title="log₂(Fold Change)",
                yaxis_title=f"-log₁₀({pval_label})",
                template="plotly_dark", width=figsize[0], height=figsize[1],
                legend_title="Best MSp"
            )
            fig_v.show()

        if show_heatmap and len(top_markers) >= 2:
            top_mz_cols = [f"{float(m):.4f}" for m in top_markers['mz']]
            df_sorted = df.sort_values('Class')
            Z = df_sorted[top_mz_cols].values.astype(float)
            means = Z.mean(axis=0); stds = Z.std(axis=0) + 1e-12
            Z_norm = np.clip((Z - means) / stds, -3, 3)

            corr_str = "BH FDR" if fdr_correction else "raw p"
            fig_h = go.Figure(go.Heatmap(
                z=Z_norm.T,
                x=list(range(len(df_sorted))),
                y=[f"{float(m):.2f}" for m in top_markers['mz']],
                colorscale='RdBu_r', zmid=0,
                colorbar=dict(title="Z-score"),
                hovertemplate="Cell %{x}<br>m/z=%{y}<br>Z=%{z:.2f}<extra></extra>"
            ))
            fig_h.update_layout(
                title=(f"Top-{len(top_markers)} marker ions — Z-score heatmap "
                       f"(sorted by MSp  |  {method.upper()}, {corr_str}<{pval_threshold}, FC≥{fc_threshold})"),
                xaxis_title="Cells (sorted by MSp)", yaxis_title="m/z",
                template="plotly_dark",
                width=figsize[0], height=max(400, 20 * len(top_markers) + 150),
                yaxis=dict(autorange='reversed')
            )
            fig_h.show()

        return top_markers, result_df

    # ──────────────────────────────────────────────────────────────────────────
    #  CLUSTER MANAGEMENT
    # ──────────────────────────────────────────────────────────────────────────

    def remove_clusters(self, clusters_to_remove, cluster_column='Class',
                        inplace=True, verbose=True):
        """
        Remove one or more cluster IDs from last_cell_spectra.

        Parameters
        ----------
        clusters_to_remove : int or list[int]
        inplace : bool — If True (default), update self.last_cell_spectra.
        """
        if self.last_cell_spectra is None:
            raise RuntimeError("No cell data. Run extract_cell_spectra() first.")
        if cluster_column not in self.last_cell_spectra.columns:
            raise ValueError(f"Column '{cluster_column}' not found.")

        if isinstance(clusters_to_remove, (int, np.integer)):
            clusters = [int(clusters_to_remove)]
        else:
            clusters = [int(c) for c in clusters_to_remove]

        df = self.last_cell_spectra.copy()
        mask = ~df[cluster_column].isin(clusters)
        filtered = df[mask].reset_index(drop=True)
        removed_count = len(df) - len(filtered)
        if verbose:
            print(f"Removed {removed_count} cells belonging to clusters {clusters}. "
                  f"Remaining cells: {len(filtered)}")
        if inplace:
            self.last_cell_spectra = filtered
        return filtered

    # ──────────────────────────────────────────────────────────────────────────
    #  EXPORT
    # ──────────────────────────────────────────────────────────────────────────

    def export_cells_to_csv(self, out_csv_path, include_spectra=True,
                            annotate=False, annotation_mode='pos',
                            annotation_ppm=10.0, annotation_top_n=3,
                            ions_only=False, round_intensities=4):
        """
        Export cell data to CSV.

        Parameters
        ----------
        out_csv_path : str         — Output file path.
        include_spectra : bool     — Include per-cell ion intensity columns.
        annotate : bool            — Annotate top differential ions (LIPID MAPS/HMDB).
        annotation_mode : str      — 'pos', 'neg', or 'both'.
        annotation_ppm : float     — Mass accuracy for annotation.
        ions_only : bool           — Export only ion intensity matrix.
        round_intensities : int    — Decimal places for intensity values.
        """
        if self.last_cell_spectra is None:
            raise RuntimeError("No cell data to export.")

        df = self.last_cell_spectra.copy()

        meta_cols_set = {'cell_id', 'area_px', 'area_um', 'centroid_row',
                         'centroid_col', 'Class', 'Score', 'umap_0', 'umap_1',
                         'soft_confidence'}
        qc_cols = [c for c in df.columns if c.startswith('qc_')]
        prob_cols = [c for c in df.columns if c.startswith('prob_cluster_')]
        mz_cols = [c for c in df.columns
                   if c not in meta_cols_set
                   and not c.startswith('qc_')
                   and not c.startswith('prob_cluster_')]

        if round_intensities is not None and mz_cols:
            df[mz_cols] = df[mz_cols].round(round_intensities)

        if ions_only:
            out_df = df[mz_cols].copy()
        elif not include_spectra:
            keep = [c for c in df.columns if c not in mz_cols]
            out_df = df[keep].copy()
        else:
            out_df = df.copy()

        os.makedirs(os.path.dirname(out_csv_path) or ".", exist_ok=True)
        out_df.to_csv(out_csv_path, index=False)

        if self.verbose:
            n_cells = len(out_df)
            n_cols = len(out_df.columns)
            size_kb = os.path.getsize(out_csv_path) // 1024 if os.path.exists(out_csv_path) else 0
            print(f"[EXPORT] CSV: {out_csv_path}  "
                  f"({n_cells} cells, {n_cols} columns, ~{size_kb} KB)")
        return out_df

    def export_cells_to_imzML(self, out_imzML_path):
        """
        Export segmented cell pixels to imzML format.

        Parameters
        ----------
        out_imzML_path : str — Output imzML path.
        """
        if self.last_cell_spectra is None:
            raise RuntimeError("No cell data to export.")
        if self.data_cube is None:
            raise RuntimeError("Data cube missing.")
        if self.cell_labels is None:
            raise RuntimeError("Cell labels missing.")

        os.makedirs(os.path.dirname(out_imzML_path) or ".", exist_ok=True)
        height, width = self.data_cube.shape[:2]

        with ImzMLWriter(out_imzML_path) as writer:
            for y in range(height):
                for x in range(width):
                    cell_id = self.cell_labels[y, x]
                    if cell_id > 0:
                        mzs = self.mz_axis
                        intens = self.data_cube[y, x, :].astype(np.float32)
                        abs_x = x + self._coord_min[0]
                        abs_y = y + self._coord_min[1]
                        writer.addSpectrum(mzs, intens, coords=(abs_x, abs_y))

        if self.verbose:
            n_pixels = np.sum(self.cell_labels > 0)
            print(f"imzML export written: {out_imzML_path} ({n_pixels} pixels exported)")


# ══════════════════════════════════════════════════════════════════════════════
#  LIPID / METABOLITE ANNOTATION
# ══════════════════════════════════════════════════════════════════════════════

# fmt: off
_ADDUCT_TABLE = {
    "[M+H]+"      : (+1,  1, +1.007276),
    "[M+Na]+"     : (+1,  1, +22.989218),
    "[M+K]+"      : (+1,  1, +38.963158),
    "[M+NH4]+"    : (+1,  1, +18.034374),
    "[M+Li]+"     : (+1,  1, +7.016003),
    "[M+2H]2+"    : (+2,  1, +1.007276),
    "[M+H-H2O]+"  : (+1,  1, +1.007276 - 18.010565),
    "[M+H-NH3]+"  : (+1,  1, +1.007276 - 17.026549),
    "[M+ACN+H]+"  : (+1,  1, +1.007276 + 41.052764),
    "[M+2Na-H]+"  : (+1,  1, +44.971160),
    "[M-H]-"      : (-1,  1, -1.007276),
    "[M+Cl]-"     : (-1,  1, +34.969402),
    "[M+FA-H]-"   : (-1,  1, +44.997654),
    "[M+Ac-H]-"   : (-1,  1, +59.013305),
    "[M-H-H2O]-"  : (-1,  1, -1.007276 - 18.010565),
    "[M+Br]-"     : (-1,  1, +78.918885),
    "[M-2H]2-"    : (-2,  1, -1.007276),
}
# fmt: on

_LIPID_DB = [
    # ── Glycerophosphocholines (PC) ──────────────────────────────────────────
    ("PC(32:0)",    "C40H80NO8P",  733.5622, "GP","PC","HMDB0000564","LPC","both"),
    ("PC(34:1)",    "C42H82NO8P",  759.5778, "GP","PC","HMDB0007859","LPC","both"),
    ("PC(34:2)",    "C42H80NO8P",  757.5622, "GP","PC","HMDB0008107","LPC","both"),
    ("PC(36:2)",    "C44H84NO8P",  785.5935, "GP","PC","HMDB0008973","LPC","both"),
    ("PC(36:4)",    "C44H80NO8P",  781.5622, "GP","PC","HMDB0011380","LPC","both"),
    ("PC(38:4)",    "C46H84NO8P",  809.5935, "GP","PC","HMDB0011065","LPC","both"),
    ("PC(38:6)",    "C46H80NO8P",  805.5622, "GP","PC","HMDB0013847","LPC","both"),
    ("PC(36:1)",    "C44H86NO8P",  787.6091, "GP","PC","HMDB0010383","LPC","both"),
    ("PC(34:0)",    "C42H84NO8P",  761.5935, "GP","PC","HMDB0015659","LPC","both"),
    ("PC(38:5)",    "C46H82NO8P",  807.5778, "GP","PC","HMDB0005372","LPC","both"),
    ("PC(40:6)",    "C48H84NO8P",  833.5935, "GP","PC","HMDB0010983","LPC","both"),
    ("PC(40:4)",    "C48H88NO8P",  837.6248, "GP","PC","HMDB0010380","LPC","both"),
    ("PC(16:0/0:0)","C24H50NO7P",  495.3353, "GP","LPC","HMDB0010382","LPC","both"),
    ("PC(18:0/0:0)","C26H54NO7P",  523.3666, "GP","LPC","HMDB0010386","LPC","both"),
    ("PC(18:1/0:0)","C26H52NO7P",  521.3510, "GP","LPC","HMDB0002815","LPC","both"),
    # ── Glycerophosphoethanolamines (PE) ─────────────────────────────────────
    ("PE(36:4)",    "C41H74NO8P",  739.5153, "GP","PE","HMDB0009398","LPE","both"),
    ("PE(38:4)",    "C43H78NO8P",  767.5466, "GP","PE","HMDB0009404","LPE","both"),
    ("PE(38:6)",    "C43H74NO8P",  763.5153, "GP","PE","HMDB0009405","LPE","both"),
    ("PE(34:1)",    "C39H76NO8P",  717.5309, "GP","PE","HMDB0009384","LPE","both"),
    ("PE(36:2)",    "C41H78NO8P",  743.5466, "GP","PE","HMDB0009397","LPE","both"),
    ("PE(34:0)",    "C39H78NO8P",  719.5466, "GP","PE","HMDB0009383","LPE","both"),
    ("PE(40:6)",    "C45H78NO8P",  791.5466, "GP","PE","HMDB0009410","LPE","both"),
    ("PE(16:0/0:0)","C21H44NO7P",  453.2883, "GP","LPE","HMDB0011503","LPE","both"),
    ("PE(18:0/0:0)","C23H48NO7P",  481.3196, "GP","LPE","HMDB0011483","LPE","both"),
    # ── Sphingomyelins (SM) ───────────────────────────────────────────────────
    ("SM(d18:1/16:0)","C39H79N2O6P", 702.5668,"SP","SM","HMDB0001348","SM","both"),
    ("SM(d18:1/18:0)","C41H83N2O6P", 730.5981,"SP","SM","HMDB0001348","SM","both"),
    ("SM(d18:1/24:1)","C47H93N2O6P", 812.6763,"SP","SM","HMDB0001348","SM","both"),
    ("SM(d18:1/24:0)","C47H95N2O6P", 814.6920,"SP","SM","HMDB0001348","SM","both"),
    # ── Ceramides (Cer) ───────────────────────────────────────────────────────
    ("Cer(d18:1/16:0)","C34H67NO3",  537.5121,"SP","Cer","HMDB0004949","Cer","both"),
    ("Cer(d18:1/18:0)","C36H71NO3",  565.5434,"SP","Cer","HMDB0006106","Cer","both"),
    ("Cer(d18:1/24:1)","C42H81NO3",  647.6217,"SP","Cer","HMDB0006107","Cer","both"),
    ("Cer(d18:1/24:0)","C42H83NO3",  649.6373,"SP","Cer","HMDB0004949","Cer","both"),
    # ── Triacylglycerols (TG) ─────────────────────────────────────────────────
    ("TG(16:0/16:0/16:0)","C51H98O6",  806.7341,"GL","TG","HMDB0005356","TG","pos"),
    ("TG(16:0/18:1/18:1)","C55H102O6", 858.7654,"GL","TG","HMDB0005357","TG","pos"),
    ("TG(18:0/18:1/18:1)","C57H106O6", 886.7967,"GL","TG","HMDB0005375","TG","pos"),
    # ── Phosphatidylserines (PS) ──────────────────────────────────────────────
    ("PS(34:1)",    "C40H76NO10P", 761.5231,"GP","PS","HMDB0008185","LPS","neg"),
    ("PS(36:1)",    "C42H80NO10P", 789.5544,"GP","PS","HMDB0008239","LPS","neg"),
    ("PS(38:4)",    "C44H78NO10P", 811.5231,"GP","PS","HMDB0008388","LPS","neg"),
    # ── Phosphatidylinositols (PI) ────────────────────────────────────────────
    ("PI(34:1)",    "C43H79O13P",  834.5295,"GP","PI","HMDB0009783","LPI","neg"),
    ("PI(36:2)",    "C45H81O13P",  860.5452,"GP","PI","HMDB0009783","LPI","neg"),
    ("PI(38:4)",    "C47H83O13P",  886.5608,"GP","PI","HMDB0009783","LPI","neg"),
    # ── Fatty acids (FA) ──────────────────────────────────────────────────────
    ("FA(16:0)",    "C16H32O2",    256.2402,"FA","FA","HMDB0000220","FA","neg"),
    ("FA(18:0)",    "C18H36O2",    284.2715,"FA","FA","HMDB0000827","FA","neg"),
    ("FA(18:1)",    "C18H34O2",    282.2559,"FA","FA","HMDB0000207","FA","neg"),
    ("FA(20:4)",    "C20H32O2",    304.2402,"FA","FA","HMDB0001043","FA","neg"),
    ("FA(22:6)",    "C22H32O2",    328.2402,"FA","FA","HMDB0002183","FA","neg"),
    # ── Cholesterol & sterols ─────────────────────────────────────────────────
    ("Cholesterol", "C27H46O",     386.3549,"ST","ST","HMDB0000067","ST","both"),
    ("Cholesteryl ester(18:1)","C45H78O2",670.5966,"ST","SE","HMDB0010359","SE","pos"),
    # ── Common metabolites ────────────────────────────────────────────────────
    ("ATP",         "C10H16N5O13P3",  506.9957,"NuA","NuA","HMDB0000538","","both"),
    ("ADP",         "C10H15N5O10P2",  427.0126,"NuA","NuA","HMDB0001341","","both"),
    ("NAD+",        "C21H27N7O14P2",  663.1096,"Vit","CoE","HMDB0000902","","both"),
    ("Glutathione", "C10H17N3O6S",    307.0838,"AA","Pep","HMDB0000125","","neg"),
    ("L-Glutamine", "C5H10N2O3",      146.0691,"AA","AA","HMDB0000641","","both"),
    ("L-Tryptophan","C11H12N2O2",     204.0899,"AA","AA","HMDB0000929","","both"),
    ("Glucose",     "C6H12O6",        180.0634,"Sug","Sug","HMDB0000122","","neg"),
    ("Cholesterol", "C27H46O",        386.3549,"ST","ST","HMDB0000067","ST","both"),
]

_PATHWAY_MAP = {
    "PC":  ["Glycerophospholipid metabolism","Membrane composition"],
    "PE":  ["Glycerophospholipid metabolism","Membrane composition","Autophagy"],
    "PS":  ["Glycerophospholipid metabolism","Apoptosis signalling"],
    "PI":  ["Glycerophospholipid metabolism","PI3K/Akt signalling"],
    "SM":  ["Sphingolipid metabolism","Membrane rafts"],
    "Cer": ["Sphingolipid metabolism","Apoptosis signalling","Ceramide signalling"],
    "TG":  ["Triacylglycerol metabolism","Lipid storage"],
    "FA":  ["Fatty acid metabolism","Beta-oxidation"],
    "LPC": ["Lysophospholipid signalling","Inflammation"],
    "ST":  ["Sterol metabolism","Membrane fluidity"],
    "SE":  ["Cholesterol esterification","Lipid storage"],
    "NuA": ["Purine metabolism","Energy metabolism"],
    "CoE": ["TCA cycle","Redox metabolism"],
}


class LipidAnnotator:
    """
    Fast exact-mass annotation of MS peaks against the embedded lipid /
    metabolite database (LIPID MAPS / HMDB).

    Supports positive and negative ionisation modes with all common adducts.

    Usage
    -----
    >>> ann = LipidAnnotator(mode='pos', ppm_tolerance=5.0)
    >>> results = ann.annotate_peak_list(mz_list)
    >>> ann.annotate_clones(pipeline)
    """

    def __init__(self, mode: str = 'pos', ppm_tolerance: float = 5.0,
                 min_dbe: float = -0.5, require_formula: bool = True):
        if mode not in ('pos', 'neg', 'both'):
            raise ValueError("mode must be 'pos', 'neg', or 'both'.")
        self.mode = mode
        self.ppm_tol = ppm_tolerance
        self.min_dbe = min_dbe
        self.require_formula = require_formula

        self._adducts = {
            name: delta
            for name, (charge, mult, delta) in _ADDUCT_TABLE.items()
            if (mode == 'both') or
               (mode == 'pos' and charge > 0) or
               (mode == 'neg' and charge < 0)
        }
        self._db = self._build_index()

    @staticmethod
    def _parse_formula(formula: str):
        import re
        pat = re.compile(r'([A-Z][a-z]?)(\d*)')
        counts = {}
        for elem, cnt in pat.findall(formula):
            counts[elem] = counts.get(elem, 0) + (int(cnt) if cnt else 1)
        return counts

    @staticmethod
    def _dbe(formula_counts: dict) -> float:
        C = formula_counts.get('C', 0)
        H = formula_counts.get('H', 0)
        N = formula_counts.get('N', 0)
        return 1 + C - H / 2.0 + N / 2.0

    def _build_index(self):
        records = []
        for row in _LIPID_DB:
            name, formula, neutral_mass, cls, subcls, hmdb, lmid, lipid_mode = row
            if lipid_mode != 'both':
                if self.mode != 'both' and lipid_mode != self.mode:
                    continue
            for adduct_name, adduct_delta in self._adducts.items():
                obs_mz = neutral_mass + adduct_delta
                if obs_mz <= 0:
                    continue
                formula_counts = self._parse_formula(formula) if formula else None
                dbe_val = self._dbe(formula_counts) if formula_counts else None
                records.append({
                    'name': name, 'formula': formula,
                    'neutral_mass': neutral_mass, 'obs_mz': obs_mz,
                    'adduct': adduct_name, 'class': cls, 'subclass': subcls,
                    'hmdb': hmdb, 'lipidmaps': lmid, 'dbe': dbe_val,
                })
        return records

    def annotate_peak(self, mz: float, min_hits: int = 1):
        """Annotate a single m/z value. Returns list of hits sorted by ppm error."""
        hits = []
        for rec in self._db:
            ppm = abs(mz - rec['obs_mz']) / rec['obs_mz'] * 1e6
            if ppm <= self.ppm_tol:
                if rec['dbe'] is not None and rec['dbe'] < self.min_dbe:
                    continue
                hit = rec.copy()
                hit['mz_query'] = mz
                hit['ppm_error'] = round(ppm, 4)
                hits.append(hit)
        hits.sort(key=lambda x: x['ppm_error'])
        return hits

    def annotate_peak_list(self, mz_array, top_n: int = 3,
                           min_intensity_array=None,
                           intensity_threshold_percentile: float = 0.0):
        """
        Annotate an array of m/z values.

        Returns pd.DataFrame with columns: mz_query, name, formula,
        neutral_mass, obs_mz, adduct, ppm_error, class, subclass,
        dbe, hmdb, lipidmaps, pathways, confidence_score.
        """
        mz_array = np.asarray(mz_array, dtype=np.float64)
        if min_intensity_array is not None and intensity_threshold_percentile > 0:
            intens = np.asarray(min_intensity_array, dtype=np.float64)
            thresh = np.percentile(intens[intens > 0], intensity_threshold_percentile)
            mz_array = mz_array[intens > thresh]

        all_hits = []
        for mz in mz_array:
            hits = self.annotate_peak(mz)[:top_n]
            all_hits.extend(hits)

        if not all_hits:
            return pd.DataFrame()

        df = pd.DataFrame(all_hits)
        df['pathways'] = df['subclass'].map(
            lambda sc: '; '.join(_PATHWAY_MAP.get(sc, [])))
        df['confidence_score'] = 1.0 / (1.0 + df['ppm_error'])
        df = df.sort_values(['mz_query', 'ppm_error']).reset_index(drop=True)
        return df

    def annotate_clones(self, pipeline, top_markers_df=None, n_top_ions: int = 30,
                        ppm_filter: float = None, show_table: bool = True,
                        show_enrichment: bool = True, figsize: tuple = (1100, 700)):
        """
        Annotate differential ions per MSp/clone.

        Workflow: extract top ions per clone → annotate → compute class
        enrichment → optional enrichment heatmap.
        """
        if pipeline.last_cell_spectra is None:
            raise RuntimeError("Run clustering first.")
        if 'Class' not in pipeline.last_cell_spectra.columns:
            raise RuntimeError("No 'Class' column — run compute_soft_clustering().")

        df = pipeline.last_cell_spectra.copy()
        meta_cols = {'cell_id','area_px','area_um','centroid_row','centroid_col',
                     'Class','Score','umap_0','umap_1','soft_confidence'}
        mz_cols = [c for c in df.columns
                   if c not in meta_cols
                   and not c.startswith('qc_') and not c.startswith('prob_cluster_')]
        mz_vals = np.array([float(c) for c in mz_cols])

        clones = sorted(df['Class'].unique())
        import plotly.colors as pc_colors
        palette = pc_colors.qualitative.Plotly
        clone_colors = {int(cl): palette[i % len(palette)] for i, cl in enumerate(clones)}

        if top_markers_df is not None:
            query_mzs = np.array(top_markers_df['mz'].values, dtype=float)
        else:
            query_mzs_list = []
            for cl in clones:
                mask = df['Class'] == cl
                X_cl = df.loc[mask, mz_cols].values.astype(float)
                X_ot = df.loc[~mask, mz_cols].values.astype(float)
                mean_cl = X_cl.mean(axis=0) + 1e-9
                mean_ot = X_ot.mean(axis=0) + 1e-9
                fc = np.log2(mean_cl / mean_ot)
                top_idx = np.argsort(np.abs(fc))[::-1][:n_top_ions]
                query_mzs_list.extend(mz_vals[top_idx].tolist())
            query_mzs = np.unique(query_mzs_list)

        ann_df = self.annotate_peak_list(query_mzs)
        if ann_df.empty:
            print(f"[ANN] No annotations found within {self.ppm_tol} ppm.")
            return ann_df, {}

        if ppm_filter is not None:
            ann_df = ann_df[ann_df['ppm_error'] <= ppm_filter].copy()

        clone_match = []
        for mz_q in ann_df['mz_query'].values:
            col = f"{mz_q:.4f}"
            if col not in df.columns:
                col = mz_cols[np.argmin(np.abs(mz_vals - mz_q))]
            mean_per_clone = {int(cl): df.loc[df['Class']==cl, col].mean() for cl in clones}
            best_clone = max(mean_per_clone, key=mean_per_clone.get)
            clone_match.append(best_clone)
        ann_df['best_clone'] = clone_match

        enrichment = {}
        for cl in clones:
            sub = ann_df[ann_df['best_clone'] == cl]
            counts = sub['subclass'].value_counts().to_dict()
            enrichment[int(cl)] = counts

        if show_table:
            print(f"\n{'='*80}")
            print(f"  LIPID / METABOLITE ANNOTATION RESULTS  "
                  f"({len(ann_df)} hits, mode={self.mode}, ppm≤{self.ppm_tol})")
            print(f"{'='*80}")
            display_cols = ['mz_query','name','adduct','ppm_error',
                            'class','subclass','confidence_score','best_clone','pathways']
            display_cols = [c for c in display_cols if c in ann_df.columns]
            with pd.option_context('display.max_rows', 60,
                                   'display.max_colwidth', 40,
                                   'display.float_format', '{:.4f}'.format):
                print(ann_df[display_cols].to_string(index=False))
            print(f"{'='*80}\n")

        if show_enrichment and len(clones) >= 2:
            all_classes = sorted({cls for cl_dict in enrichment.values() for cls in cl_dict})
            heatmap_z = np.array([[enrichment.get(int(cl), {}).get(cls, 0)
                                   for cls in all_classes] for cl in clones], dtype=float)
            col_max = heatmap_z.max(axis=0) + 1e-9
            heatmap_zn = heatmap_z / col_max

            fig = go.Figure(go.Heatmap(
                z=heatmap_zn, x=all_classes,
                y=[f"MSp {int(cl)}" for cl in clones],
                colorscale='Viridis', zmin=0, zmax=1,
                colorbar=dict(title="Rel. enrichment"),
                text=heatmap_z.astype(int).astype(str),
                texttemplate="%{text}",
                hovertemplate="MSp %{y}<br>Class: %{x}<br>Count: %{text}<extra></extra>"
            ))
            fig.update_layout(
                title="Lipid Class Enrichment per Molecular Subpopulation",
                xaxis_title="Lipid / Metabolite class",
                yaxis_title="MSp (clone)",
                template="plotly_dark",
                width=figsize[0], height=figsize[1],
                xaxis=dict(tickangle=-40)
            )
            fig.show()

        return ann_df, enrichment

    def search(self, query: str, top: int = 20):
        """Free-text search in the lipid database by name or class."""
        q = query.lower()
        hits = [r for r in self._db
                if q in r['name'].lower()
                or q in r['class'].lower()
                or q in r['subclass'].lower()]
        df_out = pd.DataFrame(hits[:top])
        if df_out.empty:
            print(f"[SEARCH] No results for '{query}'.")
        return df_out


# ── Convenience method attached to SCMSIPipeline ──────────────────────────────
def _annotate_ions(self, mode='pos', ppm_tolerance=5.0,
                   top_markers_df=None, n_top_ions=30,
                   show_table=True, show_enrichment=True):
    """
    Shortcut: annotate marker ions against LIPID MAPS / HMDB database.

    Parameters
    ----------
    mode : 'pos', 'neg', or 'both'
    ppm_tolerance : float — Mass accuracy in ppm.
    """
    annotator = LipidAnnotator(mode=mode, ppm_tolerance=ppm_tolerance)
    return annotator.annotate_clones(
        self, top_markers_df=top_markers_df, n_top_ions=n_top_ions,
        show_table=show_table, show_enrichment=show_enrichment
    )

SCMSIPipeline.annotate_ions = _annotate_ions

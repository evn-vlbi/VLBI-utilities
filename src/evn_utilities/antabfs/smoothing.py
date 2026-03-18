"""Robust smoothing for Tsys data.

Provides sigma-clipping and median-based smoothing that is resilient
to large outliers commonly found in raw Tsys measurements.  Enabled
by default in non-interactive (batch) mode.

Strategy:
    1. Sigma-clip using the MAD (median absolute deviation) as a robust
       scale estimator — this flags extreme outliers without being affected
       by them.
    2. Replace flagged values with a local median interpolation.
    3. Apply a Savitzky-Golay (or rolling-median) smoothing pass to produce
       a clean Tsys curve that preserves genuine trends.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from functools import partial

import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import savgol_filter

_MAD_SCALE = 1.4826


def mad(data: np.ndarray) -> float:
    """Median Absolute Deviation (MAD), scaled to match Gaussian sigma."""
    med = np.median(data)
    return float(_MAD_SCALE * np.median(np.abs(data - med)))


def sigma_clip(
    data: np.ndarray,
    sigma: float = 3.0,
    max_iter: int = 5,
) -> np.ndarray:
    """Return a boolean mask where ``True`` marks *good* (non-outlier) values.

    Uses iterative MAD-based clipping.  Convergence is tested by comparing
    good-value counts (cheaper than ``np.array_equal``).  ``np.median`` is
    used instead of ``np.nanmedian`` on the already-masked subset, avoiding
    the NaN-scanning overhead.
    """
    mask = np.isfinite(data) & (data > 0)
    n_good = int(mask.sum())
    for _ in range(max_iter):
        if n_good < 3:
            break
        sub = data[mask]                        # one copy, reused below
        med = float(np.median(sub))
        scale = float(_MAD_SCALE * np.median(np.abs(sub - med)))
        if scale == 0.0:
            break
        new_mask = mask & (np.abs(data - med) < sigma * scale)
        n_new = int(new_mask.sum())
        if n_new == n_good:                     # converged — count comparison
            break
        mask, n_good = new_mask, n_new
    return mask


def smooth_tsys_channel(values: np.ndarray, times: np.ndarray | None = None,
                        block: np.ndarray | None = None, sigma: float = 3.0,
                        window: int = 0, method: str = "savgol") -> np.ndarray:
    """Smooth a single BBC-channel Tsys time series.

    Parameters
    ----------
    values : 1-D array
        Raw Tsys measurements.
    times : 1-D array or None
        Timestamps (unused; kept for API compatibility).
    block : 1-D array or None
        Scan-block indices.  Per-block sigma-clipping when provided.
    sigma : float
        Clipping threshold in MAD-scaled sigma units.
    window : int
        Smoothing window size (0 = auto: max(5, n//20)).
    method : str
        ``"savgol"`` for Savitzky-Golay, ``"median"`` for rolling median.

    Returns
    -------
    ndarray
        Smoothed Tsys values with outliers replaced.
    """
    y = np.array(values, dtype=float)
    n = len(y)
    if n < 3:
        return y

    # ------------------------------------------------------------------
    # Step 1: per-block sigma-clip with global-scale fallback for small
    # blocks (< 3 points cannot form a reliable local scale estimate).
    # ------------------------------------------------------------------
    good = np.ones(n, dtype=bool)
    if block is not None:
        blk = np.asarray(block)
        processed = np.zeros(n, dtype=bool)
        for b in np.unique(blk):
            idx = np.where(blk == b)[0]
            if len(idx) < 3:
                continue
            good[idx] = sigma_clip(y[idx], sigma=sigma)
            processed[idx] = True

        unprocessed = ~processed
        if unprocessed.any() and processed.any():
            proc_good = good & processed
            if proc_good.sum() >= 3:
                sub = y[proc_good]
                gmed = float(np.median(sub))
                gscale = float(_MAD_SCALE * np.median(np.abs(sub - gmed)))
                if gscale > 0.0:
                    good[unprocessed] = np.abs(y[unprocessed] - gmed) < sigma * gscale
    else:
        good = sigma_clip(y, sigma=sigma)

    # ------------------------------------------------------------------
    # Step 2: replace bad values with the local median of good neighbours.
    # Pre-build a sorted array of good positions and use searchsorted for
    # O(log n) window lookup instead of slicing + boolean masking.
    # ------------------------------------------------------------------
    y_clean = y.copy()
    bad_idx = np.where(~good)[0]
    if len(bad_idx) and good.any():
        global_median = float(np.median(y[good]))
        good_pos = np.where(good)[0]            # sorted good-value positions
        for bi in bad_idx:
            lo = int(np.searchsorted(good_pos, bi - 20, "left"))
            hi = int(np.searchsorted(good_pos, bi + 21, "right"))
            nearby = good_pos[lo:hi]
            y_clean[bi] = float(np.median(y[nearby])) if len(nearby) else global_median

    # ------------------------------------------------------------------
    # Step 3: smoothing pass.
    # ``window |= 1`` sets the LSB in one bitwise op, guaranteeing odd.
    # ``np.clip`` is called with ``out=smoothed`` for an in-place update.
    # ------------------------------------------------------------------
    if window == 0:
        window = max(5, n // 20)
    window |= 1                                 # ensure odd
    window = min(window, n)

    if method == "savgol" and n >= window >= 5:
        poly_order = min(3, window - 1)
        try:
            smoothed = savgol_filter(y_clean, window, poly_order)
        except Exception:
            smoothed = median_filter(y_clean, size=window)
    else:
        smoothed = median_filter(y_clean, size=window)

    np.clip(smoothed, 0.1, None, out=smoothed)  # in-place
    return smoothed


def smooth_tsys_matrix(tsys_matrix: np.ndarray, times: np.ndarray | None = None,
                       block: np.ndarray | None = None, sigma: float = 3.0,
                       window: int = 0, method: str = "savgol") -> np.ndarray:
    """Smooth all BBC channels in parallel using a thread pool.

    Each channel row is independent, so they are distributed across
    worker threads.  NumPy / SciPy release the GIL for most operations,
    allowing true parallel execution even in CPython.

    Parameters
    ----------
    tsys_matrix : ndarray, shape (n_channels, n_times)
        One row per BBC channel.
    times, block, sigma, window, method
        Passed to :func:`smooth_tsys_channel`.

    Returns
    -------
    ndarray
        Smoothed matrix with the same shape.
    """
    n_ch = tsys_matrix.shape[0]
    _fn = partial(smooth_tsys_channel, times=times, block=block, sigma=sigma, window=window, method=method)
    max_workers = min(n_ch, os.cpu_count() or 1)
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        rows = list(pool.map(_fn, tsys_matrix))
    return np.vstack(rows)

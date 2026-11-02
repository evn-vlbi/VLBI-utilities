"""Robust smoothing for Tsys data.

Provides sigma-clipping and median-based smoothing that is resilient
to large outliers commonly found in raw Tsys measurements.  Enabled
by default in non-interactive mode.

Smoothing is performed **per scan**: each scan block is processed on its
own and no information ever crosses a scan boundary.  This way a noisy or
mispointed scan cannot distort its neighbours, and genuine Tsys level
changes between scans are preserved as hard steps instead of being blurred.

Strategy, applied independently to every scan:
    1. Sigma-clip using the MAD (median absolute deviation) as a robust
       scale estimator — this flags extreme outliers without being affected
       by them.
    2. Replace flagged values with the scan's robust median so the smoother
       is not pulled by the outliers.
    3. Fit a straight line over the scan (degree-1 Savitzky-Golay) — or a
       rolling median — to produce a clean within-scan trend.  A line cannot
       bend, so no spurious parabolas appear on short, noisy scans.  Scans too
       short to smooth collapse to their robust median (a flat level).
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from functools import partial

import numpy as np
from scipy.ndimage import median_filter

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


def _smooth_scan(y: np.ndarray, sigma: float, window: int, method: str) -> np.ndarray:
    """Smooth one scan's Tsys series, robust to outliers, with no cross-scan info.

    Outliers are sigma-clipped within this scan only, then a robust straight-line
    fit (or rolling median) over the remaining good points produces a smooth
    within-scan trend.  Scans with too few good points to fit collapse to their
    robust median (a flat level).

    Parameters
    ----------
    y : 1-D array
        Raw Tsys values for a single scan.
    sigma, window, method
        As documented on :func:`smooth_tsys_channel`.

    Returns
    -------
    ndarray
        Smoothed values for the scan, same length as *y*.
    """
    m = len(y)
    if m == 0:
        return y.copy()

    good = sigma_clip(y, sigma=sigma) if m >= 3 else (np.isfinite(y) & (y > 0))
    if not good.any():
        return y.copy()

    scan_med = float(np.median(y[good]))
    if m < 5 or int(good.sum()) < 5:
        return np.full(m, scan_med)             # too few points → flat level

    x = np.arange(m, dtype=float)

    if method != "median":
        # Default ("savgol"): a least-squares straight line — i.e. a degree-1
        # Savitzky-Golay over the whole scan — fitted to the good (non-outlier)
        # points only and evaluated at every sample.  A line cannot bend, so it
        # never invents the spurious parabolas a higher-order fit produces on
        # short, noisy scans (here a cubic differs from a line by only ~0.06 K
        # yet bends in ~80% of scans), and outliers are excluded from the fit.
        xoff = float(x[good].mean())
        slope, intercept = np.polyfit(x[good] - xoff, y[good], 1)
        return slope * (x - xoff) + intercept

    # method == "median": rolling median over the scan; outliers are first
    # replaced by the scan median so they cannot leak into the window result.
    y_clean = y.copy()
    y_clean[~good] = scan_med
    w = window if window > 0 else m
    w = min(w, m)
    if w % 2 == 0:
        w -= 1
    return median_filter(y_clean, size=max(w, 3))


def smooth_tsys_channel(values: np.ndarray, times: np.ndarray | None = None,
                        block: np.ndarray | None = None, sigma: float = 3.0,
                        window: int = 0, method: str = "savgol") -> np.ndarray:
    """Smooth a single BBC-channel Tsys time series, scan by scan.

    Every scan (block) is smoothed independently by :func:`_smooth_scan`:
    outlier rejection and the smoothing pass never cross a scan boundary, so a
    noisy scan cannot distort its neighbours and genuine level changes between
    scans are kept as hard steps.

    Parameters
    ----------
    values : 1-D array
        Raw Tsys measurements.
    times : 1-D array or None
        Timestamps (unused; kept for API compatibility).
    block : 1-D array or None
        Scan-block index per sample.  When None the whole series is treated as
        a single scan.
    sigma : float
        Clipping threshold in MAD-scaled sigma units.
    window : int
        Smoothing window size (0 = auto: the full length of each scan).
    method : str
        ``"savgol"`` for Savitzky-Golay, ``"median"`` for rolling median.

    Returns
    -------
    ndarray
        Smoothed Tsys values with outliers replaced, clipped to >= 0.1 K.
    """
    y = np.array(values, dtype=float)
    n = len(y)
    if n < 3:
        return y

    if block is None:
        out = _smooth_scan(y, sigma, window, method)
    else:
        out = y.copy()
        blk = np.asarray(block)
        for b in np.unique(blk):
            idx = np.where(blk == b)[0]
            out[idx] = _smooth_scan(y[idx], sigma, window, method)

    np.clip(out, 0.1, None, out=out)
    return out


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

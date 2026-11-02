"""Robust smoothing for Tsys data.

Provides sigma-clipping and robust smoothing that is resilient to large
outliers commonly found in raw Tsys measurements.  Enabled by default in
non-interactive mode.  Two spans are available:

``span="scan"`` — every scan block is processed on its own and no
information crosses a scan boundary.  Strategy, per scan:
    1. Sigma-clip using the MAD (median absolute deviation) as a robust
       scale estimator.
    2. Replace flagged values with the scan's robust median so the smoother
       is not pulled by the outliers.
    3. Fit a straight line over the scan (degree-1 Savitzky-Golay) — or a
       rolling median — to produce a clean within-scan trend.  Scans too
       short to smooth collapse to their robust median (a flat level).

``span="source"`` — Tsys must evolve smoothly in time for one source, but
never mixes sources.  Scans are grouped by source (from the log's
``:source=`` lines) into *segments*; a new segment starts when the same
source returns after more than ``max_gap`` seconds.  Strategy, per segment:
    1. Robust median (sigma-clipped) and mean time of every scan.
    2. Scan medians that deviate by more than ``sigma`` robust sigmas from a
       robust (Theil-Sen) line through their 4 nearest same-source scans are
       dropped (bad scans); their samples get the curve value from the
       neighbouring good scans.
    3. Smoothing spline in time through the remaining scan medians, weighted
       by the number of good samples (``scipy.interpolate.make_smoothing_spline``).
       The smoothing level is the smoothest curve whose rms deviation from the
       medians stays within their scan-to-scan scatter (estimated in step 2).
       Segments with 2-4 scans are joined by straight segments between scan
       medians; a lone scan is smoothed exactly as with ``span="scan"``.
    4. The curve is evaluated at every sample of the segment (clamped to the
       first/last scan-median time, so it never extrapolates).
    5. Guard: if the curve misses any accepted scan median by more than 20%
       of its level, Tsys is not smooth for this source (e.g. it jumps between
       very different levels from scan to scan) and the segment is smoothed
       per scan instead (with ``method`` / ``window``).  Such fallbacks are
       counted and reported by :func:`smooth_tsys_matrix`.
    Otherwise ``method`` and ``window`` are not used in this span.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from functools import partial

import numpy as np
from scipy.interpolate import make_smoothing_spline
from scipy.ndimage import median_filter

_MAD_SCALE = 1.4826
SMOOTH_SPANS = ("scan", "source")
_MIN_SPLINE_POINTS = 5  # make_smoothing_spline needs at least 5 points; fewer -> straight line / flat level
_MIN_SCATTER_FRAC = 0.005  # scan-median scatter floor, as a fraction of the Tsys level (0.25 K at 50 K)
_LOG_LAM_MIN, _LOG_LAM_MAX, _BISECT_STEPS = -8.0, 8.0, 16  # smoothing-parameter search (x in hours)
_MAX_MISS_FRAC = 0.20  # curve missing an accepted scan median by more than this -> segment smoothed per scan


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

    Uses iterative MAD-based clipping.  Non-finite and non-positive values are
    always bad.  Convergence is tested by comparing good-value counts.
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


# ---------------------------------------------------------------------- #
# span = "scan"
# ---------------------------------------------------------------------- #

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


# ---------------------------------------------------------------------- #
# span = "source"
# ---------------------------------------------------------------------- #

def source_segments(times: np.ndarray, block: np.ndarray, scan_sources: dict[int, str],
                    max_gap: float) -> list[list[np.ndarray]]:
    """Group sample indices into same-source segments of scans.

    Parameters
    ----------
    times : Unix timestamps per sample.
    block : scan number per sample.
    scan_sources : ``{scan number: source name}``; scans missing from it form their own one-scan segment.
    max_gap : seconds; a same-source scan starting more than this after the previous one starts a new segment.

    Returns
    -------
    list of segments; each segment is a list of per-scan index arrays, in time order.
    """
    t = np.asarray(times, dtype=float)
    blk = np.asarray(block)
    by_source: dict[str, list[tuple[float, float, np.ndarray]]] = {}
    for b in np.unique(blk):
        idx = np.where(blk == b)[0]
        source = scan_sources.get(int(b), f"<no source: scan {int(b)}>")
        by_source.setdefault(source, []).append((float(t[idx].min()), float(t[idx].max()), idx))

    segments: list[list[np.ndarray]] = []
    for scans in by_source.values():
        scans.sort(key=lambda s: s[0])
        segment = [scans[0][2]]
        last_end = scans[0][1]
        for start, end, idx in scans[1:]:
            if start - last_end > max_gap:
                segments.append(segment)
                segment = []
            segment.append(idx)
            last_end = end
        segments.append(segment)
    return segments


def _smoothest_spline(x: np.ndarray, y: np.ndarray, w: np.ndarray, target_rms: float):
    """Smoothest spline whose weighted rms residual to (x, y) does not exceed *target_rms* (Reinsch criterion).

    The residual grows monotonically with the smoothing parameter lam, so lam is found by bisection in log10
    space over [_LOG_LAM_MIN, _LOG_LAM_MAX].  (GCV is unreliable with the few points a segment usually has: it
    tends to interpolate the noisy scan medians.)
    """
    def fit(log_lam: float):
        spline = make_smoothing_spline(x, y, w=w, lam=10.0 ** log_lam)
        return spline, float(np.sqrt(np.sum(w * (y - spline(x)) ** 2) / np.sum(w)))

    lo, hi = _LOG_LAM_MIN, _LOG_LAM_MAX
    spline_hi, rms_hi = fit(hi)
    if rms_hi <= target_rms:
        return spline_hi                                   # even a near-straight line is good enough
    spline_lo, rms_lo = fit(lo)
    if rms_lo > target_rms:
        return spline_lo                                   # nothing reaches the target: follow the data closely
    for _ in range(_BISECT_STEPS):
        mid = 0.5 * (lo + hi)
        spline_mid, rms_mid = fit(mid)
        if rms_mid <= target_rms:
            lo, spline_lo = mid, spline_mid
        else:
            hi = mid
    return spline_lo


def _fit_curve(x: np.ndarray, y: np.ndarray, w: np.ndarray, target_rms: float):
    """Return a callable ``f(x_new)`` fitted to (x, y, weights w).

    >= _MIN_SPLINE_POINTS points: smoothest spline within *target_rms* of the points; 2-4: straight segments
    joining the points (too few to estimate a smoothing level, and a single line through 3-4 scans can miss every
    scan by several K); 1: flat level.  Outside [x[0], x[-1]] the value is held constant (no extrapolation).
    *x* must be strictly increasing.
    """
    if len(x) >= _MIN_SPLINE_POINTS:
        spline = _smoothest_spline(x, y, w, target_rms)
        return lambda xn: spline(np.clip(xn, x[0], x[-1]))
    if len(x) >= 2:
        return lambda xn: np.interp(xn, x, y)
    return lambda xn: np.full(np.shape(xn), float(y[0]))


def _neighbour_prediction(x: np.ndarray, med: np.ndarray, i: int) -> float:
    """Predict med[i] from its 4 nearest neighbouring scans (excluding i) with a Theil-Sen line.

    The 5-scan window around i is shifted to stay inside the segment, so edge scans are predicted by
    extrapolating the trend of the next scans instead of being biased by it.  Theil-Sen (median of pairwise
    slopes) is robust to one bad neighbour.  Requires len(x) >= 5.
    """
    start = min(max(i - 2, 0), len(x) - 5)
    nb = [j for j in range(start, start + 5) if j != i]
    xs, ys = x[nb], med[nb]
    slopes = [(ys[b] - ys[a]) / (xs[b] - xs[a]) for a in range(4) for b in range(a + 1, 4)]
    slope = float(np.median(slopes))
    intercept = float(np.median(ys - slope * xs))
    return slope * float(x[i]) + intercept


def _good_scan_medians(x: np.ndarray, med: np.ndarray, sigma: float) -> tuple[np.ndarray, float]:
    """Return (mask of scan medians consistent with their neighbours, scatter of one scan median in K).

    The reference for scan i is predicted from its neighbours only (:func:`_neighbour_prediction`), *excluding
    scan i itself*: a reference that includes the point reproduces smooth data exactly, which makes most
    residuals 0, the MAD 0, and the test meaningless.  A median is rejected if it deviates from its reference by
    more than *sigma* robust sigmas of all such deviations.  The scale is floored at _MIN_SCATTER_FRAC of the
    typical level and converted to the scatter of a single median (the prediction adds ~30% of the variance).
    With fewer than _MIN_SPLINE_POINTS scans nothing is rejected and the scatter is 0.
    *x* must be strictly increasing.
    """
    n = len(med)
    if n < _MIN_SPLINE_POINTS:
        return np.ones(n, dtype=bool), 0.0
    ref = np.array([_neighbour_prediction(x, med, i) for i in range(n)])
    resid = med - ref
    scale = max(mad(resid), _MIN_SCATTER_FRAC * float(np.median(np.abs(med))))
    if scale == 0.0:
        return np.ones(n, dtype=bool), 0.0
    return np.abs(resid - np.median(resid)) <= sigma * scale, scale / np.sqrt(1.3)


def _smooth_segment(y: np.ndarray, t: np.ndarray, segment: list[np.ndarray], sigma: float) -> np.ndarray | None:
    """Smooth one same-source segment: spline through robust scan medians, evaluated at every sample.

    Returns the smoothed values for ``np.concatenate(segment)`` (same order), or None if no scan has good data
    or if the data are not smooth for this source: the curve misses an accepted scan median by more than
    _MAX_MISS_FRAC of its level (e.g. Tsys jumping between very different levels from scan to scan).  The caller
    then smooths each scan of the segment on its own.
    """
    med_t: list[float] = []
    med_y: list[float] = []
    weight: list[float] = []
    for idx in segment:
        ys = y[idx]
        good = sigma_clip(ys, sigma=sigma) if len(ys) >= 3 else (np.isfinite(ys) & (ys > 0))
        if good.any():
            med_t.append(float(np.mean(t[idx][good])))
            med_y.append(float(np.median(ys[good])))
            weight.append(float(good.sum()))
    if not med_t:
        return None

    order = np.argsort(med_t)
    x0 = med_t[order[0]]
    x = (np.asarray(med_t)[order] - x0) / 3600.0          # hours since first scan median: well-conditioned
    yv = np.asarray(med_y)[order]
    wv = np.asarray(weight)[order]
    keep = np.concatenate(([True], np.diff(x) > 0))      # spline needs strictly increasing x
    x, yv, wv = x[keep], yv[keep], wv[keep] / wv[keep].mean()

    # Reject outlier scans (e.g. RFI or a mispointed scan) BEFORE fitting: the reference comes from a robust line
    # through the neighbouring scans, which a single bad scan cannot pull (a fitted curve can, and then rejects
    # its good neighbours instead).
    good, scatter = _good_scan_medians(x, yv, sigma)
    curve = _fit_curve(x[good], yv[good], wv[good], target_rms=scatter)
    if np.any(np.abs(curve(x[good]) - yv[good]) > _MAX_MISS_FRAC * np.abs(yv[good])):
        return None

    idx_all = np.concatenate(segment)
    return curve((t[idx_all] - x0) / 3600.0)


# ---------------------------------------------------------------------- #
# Public API
# ---------------------------------------------------------------------- #

def smooth_tsys_channel(values: np.ndarray, times: np.ndarray | None = None,
                        block: np.ndarray | None = None, sigma: float = 3.0,
                        window: int = 0, method: str = "savgol", span: str = "scan",
                        segments: list[list[np.ndarray]] | None = None,
                        fallbacks: list[int] | None = None) -> np.ndarray:
    """Smooth a single BBC-channel Tsys time series.

    Parameters
    ----------
    values : 1-D array
        Raw Tsys measurements (NaN / <= 0 values are ignored by the fits).
    times : 1-D array or None
        Unix timestamps; required for ``span="source"``.
    block : 1-D array or None
        Scan-block index per sample.  When None the whole series is treated as
        a single scan (``span="scan"`` only).
    sigma : float
        Clipping threshold in MAD-scaled sigma units.
    window : int
        Smoothing window size (0 = the full length of each scan); ``span="scan"`` with ``"median"`` only.
    method : str
        ``"savgol"`` for a per-scan straight line, ``"median"`` for rolling median (``span="scan"`` only).
    span : str
        ``"scan"`` (each scan alone) or ``"source"`` (smooth across same-source scans, see module docstring).
    segments : list or None
        Output of :func:`source_segments`; required for ``span="source"``.
    fallbacks : list or None
        If given, the indices (into *segments*) of segments that were smoothed per scan instead (span="source",
        Tsys not smooth for that source) are appended to it.

    Returns
    -------
    ndarray
        Smoothed Tsys values, clipped to >= 0.1 K.  Samples that cannot be fitted keep their input value.
    """
    y = np.array(values, dtype=float)
    n = len(y)
    if n < 3:
        return y

    if span == "source":
        if times is None or segments is None:
            raise ValueError('span="source" needs times and segments')
        t = np.asarray(times, dtype=float)
        out = y.copy()
        for seg_i, segment in enumerate(segments):
            # A lone scan has nothing to share with: smooth it exactly as span="scan" does (keeps its slope).
            vals = _smooth_segment(y, t, segment, sigma) if len(segment) > 1 else None
            if vals is not None:
                out[np.concatenate(segment)] = vals
                continue
            if len(segment) >= _MIN_SPLINE_POINTS and fallbacks is not None:
                fallbacks.append(seg_i)
            for idx in segment:
                out[idx] = _smooth_scan(y[idx], sigma, window, method)
    elif block is None:
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
                       window: int = 0, method: str = "savgol", span: str = "scan",
                       scan_sources: dict[int, str] | None = None, max_gap: float = 3600.0,
                       labels: list[str] | None = None, report: bool = False) -> np.ndarray:
    """Smooth all BBC channels in parallel using a thread pool.

    Parameters
    ----------
    tsys_matrix : ndarray, shape (n_channels, n_times)
        One row per BBC channel.
    times, block, sigma, window, method, span
        Passed to :func:`smooth_tsys_channel`.
    scan_sources : dict or None
        ``{scan number: source name}`` for ``span="source"`` (missing scans are smoothed on their own).
    max_gap : float
        Seconds; same-source scans further apart are smoothed as separate segments (``span="source"``).
    labels : list or None
        Channel labels used in the report (default: ``ch1``, ``ch2``...).
    report : bool
        Print which sources/channels fell back to per-scan smoothing (``span="source"``).

    Returns
    -------
    ndarray
        Smoothed matrix with the same shape.
    """
    if span not in SMOOTH_SPANS:
        raise ValueError(f"span must be one of {SMOOTH_SPANS}, got {span!r}")
    segments = None
    if span == "source":
        if times is None or block is None:
            raise ValueError('span="source" needs times and block')
        segments = source_segments(times, block, scan_sources or {}, max_gap)
    n_ch = tsys_matrix.shape[0]
    fallbacks: list[list[int]] = [[] for _ in range(n_ch)]
    _fn = partial(smooth_tsys_channel, times=times, block=block, sigma=sigma, window=window, method=method,
                  span=span, segments=segments)
    max_workers = min(n_ch, os.cpu_count() or 1)
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        rows = list(pool.map(lambda i: _fn(tsys_matrix[i], fallbacks=fallbacks[i]), range(n_ch)))
    if report and segments is not None:
        _report_fallbacks(segments, fallbacks, np.asarray(block), scan_sources or {}, labels)
    return np.vstack(rows)


def _report_fallbacks(segments: list[list[np.ndarray]], fallbacks: list[list[int]], block: np.ndarray,
                      scan_sources: dict[int, str], labels: list[str] | None) -> None:
    """Print one line per source whose Tsys was not smooth (segment smoothed per scan), with the channels."""
    per_source: dict[str, set[str]] = {}
    for ch_i, seg_list in enumerate(fallbacks):
        label = labels[ch_i] if labels and ch_i < len(labels) else f"ch{ch_i + 1}"
        for seg_i in seg_list:
            source = scan_sources.get(int(block[segments[seg_i][0][0]]), "?")
            per_source.setdefault(source, set()).add(label)
    for source, chans in sorted(per_source.items()):
        print(f"  Smooth span source: Tsys of {source} is not smooth across its scans (curve misses a scan by "
              f"> {_MAX_MISS_FRAC:.0%}) in {', '.join(sorted(chans))} -> smoothed per scan there.")

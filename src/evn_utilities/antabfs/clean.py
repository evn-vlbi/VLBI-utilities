"""Detail-preserving Tsys cleaning (smoothing method ``"clean"``).

Goal: remove RFI outliers and recover the true Tsys level, while keeping every sample that looks like a real
measurement *unchanged* (second-scale structure is not averaged away).  Output = baseline + kept detail.

Assumptions about the data (checked on EVN/LBA FS logs, 1 s and 15 s cadence):
    * Clean Tsys has a small, nearly constant *relative* noise per channel (radiometer + cal noise, ~0.5-3%).
    * RFI appears as spikes of either sign (with continuous cal, RFI in one cal phase can push Tsys down too),
      dense bursts, and samples the FS logs as <= 0 or absurdly large.  Those invalid samples, and any sample the
      caller passes as NaN (e.g. manual flags), are *missing*, never data.
    * Slow trends are real: they must be followed, never flattened.

Algorithm, per channel:
    1. Local statistics in a sliding window of ``window_s`` seconds (default 60 s; converted to samples with each
       scan's own cadence, minimum 5):
         location = shorth (centre of the shortest interval holding half the samples): it locks onto the dense
                    clean cluster even with ~50% lopsided contamination, where a median is dragged towards junk;
         noise    = MAD of first differences / sqrt(2): insensitive to slow trends, inflated by spikes/bursts.
    2. Clean relative noise of the channel ``rel_clean`` = 10th percentile of noise/location over all windows,
       floored at ``REL_NOISE_FLOOR`` (0.5%; Tsys is never calibrated better than that, and stations that log
       very stable Tsys would otherwise get an absurdly tight clip).
    3. A window is *trusted* if its relative noise <= ``TRUST_FACTOR`` * rel_clean.  Per scan, the first rule that
       applies gives the baseline level ``base0`` and the noise ``sigma_k`` (in K) used for clipping:
         "window": the scan has trusted windows -> base0 / noise are the trusted windows' location / noise,
                   interpolated in time across untrusted stretches (contaminated stretches never set the level);
         "scan":   no trusted window, but the whole scan has relative noise <= ``TRUST_SCAN_FACTOR`` * rel_clean
                   -> the scan's own level and noise (a noisier but coherent scan keeps its detail);
         "source": (span "source" only) the level is interpolated in time from same-source scans that have a
                   level, if the scan's own valid samples support it at least half as well as their own shorth
                   level (else the scan keeps its own level, as "own");
         "own":    otherwise -> the scan's shorth level with the strict clean noise;
         "empty":  no valid sample and nothing to borrow -> NaN (the caller keeps its input values).
    4. Clip: samples with |y - base0| > ``sigma`` * noise are outliers.
    5. Final baseline (kinds "window"/"scan"): rolling median of the kept samples, again only from trusted windows
       and interpolated across the rest, then a short boxcar (``window/3``) to remove median-filter steps.
    6. Output: kept samples in trusted windows -> the measured value, unchanged.  Kept samples in untrusted windows
       -> baseline + w * (value - baseline), w = min(1, (TRUST_FACTOR * noise / local scatter)^2) with the local
       scatter = MAD of kept residuals over ``DETAIL_SAMPLES`` samples (shrinks burst-inflated scatter towards
       the baseline).  Outliers and missing samples -> baseline.
"""

from __future__ import annotations

import warnings

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.ndimage import uniform_filter1d

MAX_VALID_TSYS = 10000.0     # K; larger (or <= 0) logged values are not measurements (same limit as prefilter)
DEFAULT_WINDOW_S = 60.0      # s; local statistics window when the caller passes window <= 0
MIN_WINDOW_SAMPLES = 5       # shorth needs >= 5 samples
REL_NOISE_FLOOR = 0.005      # floor of the clean relative noise (0.5%)
CLEAN_PERCENTILE = 10.0      # percentile of the window relative noise taken as the channel's clean noise
TRUST_FACTOR = 2.5           # window trusted if its relative noise <= TRUST_FACTOR * clean noise
TRUST_SCAN_FACTOR = 5.0      # whole scan accepted as its own level if relative noise <= this * clean noise
SOURCE_SUPPORT_FRAC = 0.5    # a borrowed same-source level must explain >= this fraction of the own-level support
DETAIL_SAMPLES = 15          # samples in the local-scatter window used to shrink detail in untrusted stretches
_MAD_SCALE = 1.4826
_SHORTH_TO_SIGMA = 1.0 / 1.349  # width of the shortest half of a Gaussian sample = 1.349 sigma

KINDS = ("window", "scan", "source", "own", "empty")


# ---------------------------------------------------------------------- #
# Rolling robust statistics (NaN = missing)
# ---------------------------------------------------------------------- #

def _windows(y: np.ndarray, w: int) -> np.ndarray:
    """Centred sliding windows of odd length *w* over *y* (NaN-padded at the edges), shape (len(y), w)."""
    h = w // 2
    padded = np.concatenate([np.full(h, np.nan), y, np.full(h, np.nan)])
    return sliding_window_view(padded, 2 * h + 1)


def _nanmedian_rows(v: np.ndarray) -> np.ndarray:
    """Row-wise median ignoring NaN; NaN for all-NaN rows (without the numpy warning)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(v, axis=1)


def rolling_median(y: np.ndarray, w: int) -> np.ndarray:
    """Centred rolling median of odd window *w*, ignoring NaN (fewer samples at the edges)."""
    return _nanmedian_rows(_windows(y, w))


def _shorth_rows(v_sorted: np.ndarray) -> np.ndarray:
    """Shorth location per row of an ascending-sorted 2-D array with NaN last; NaN if < MIN_WINDOW_SAMPLES finite."""
    n_finite = np.sum(np.isfinite(v_sorted), axis=1)
    loc = np.full(len(v_sorted), np.nan)
    for n in np.unique(n_finite):               # few distinct counts (scan edges, missing samples)
        if n < MIN_WINDOW_SAMPLES:
            continue
        rows = np.where(n_finite == n)[0]
        h = n // 2 + 1
        vv = v_sorted[rows, :n]
        width = vv[:, h - 1:] - vv[:, :n - h + 1]
        j = np.argmin(width, axis=1)
        ar = np.arange(len(rows))
        loc[rows] = 0.5 * (vv[ar, j] + vv[ar, j + h - 1])
    return loc


def _diff_noise_rows(d_windows: np.ndarray) -> np.ndarray:
    """Per-row noise from first differences: MAD(d) / sqrt(2); NaN with fewer than 4 differences."""
    med = _nanmedian_rows(d_windows)
    noise = _MAD_SCALE * _nanmedian_rows(np.abs(d_windows - med[:, None])) / np.sqrt(2.0)
    noise[np.sum(np.isfinite(d_windows), axis=1) < MIN_WINDOW_SAMPLES - 1] = np.nan
    return noise


def _first_differences(y: np.ndarray) -> np.ndarray:
    """Difference of each finite sample to the previous finite one (NaN elsewhere and at the first finite one)."""
    d = np.full(len(y), np.nan)
    idx = np.where(np.isfinite(y))[0]
    if len(idx) > 1:
        d[idx[1:]] = np.diff(y[idx])
    return d


def rolling_location_noise(y: np.ndarray, w: int) -> tuple[np.ndarray, np.ndarray]:
    """Rolling (shorth location, first-difference noise) over odd window *w*; NaN = missing / not enough samples."""
    loc = _shorth_rows(np.sort(_windows(y, w), axis=1))
    noise = _diff_noise_rows(_windows(_first_differences(y), w))
    return loc, noise


def location_noise(y: np.ndarray) -> tuple[float, float]:
    """(shorth location, first-difference noise) of all finite samples of *y*; NaN if too few samples."""
    finite = y[np.isfinite(y)]
    if len(finite) < MIN_WINDOW_SAMPLES:
        return float("nan"), float("nan")
    loc = float(_shorth_rows(np.sort(finite)[None, :])[0])
    noise = float(_diff_noise_rows(np.diff(finite)[None, :])[0])
    return loc, noise


# ---------------------------------------------------------------------- #
# Per-scan steps
# ---------------------------------------------------------------------- #

def valid_mask(y: np.ndarray) -> np.ndarray:
    """True where *y* is a usable measurement: finite, > 0 and < MAX_VALID_TSYS."""
    return np.isfinite(y) & (y > 0) & (y < MAX_VALID_TSYS)


def window_samples(t: np.ndarray, window_s: float, n_max: int) -> int:
    """Odd number of samples spanning *window_s* seconds at the median cadence of *t*, in [5, n_max] if possible."""
    dt = float(np.median(np.diff(t))) if len(t) > 1 else 1.0
    n = max(MIN_WINDOW_SAMPLES, int(round(window_s / max(dt, 1e-3))))
    n = min(n, max(n_max, 1))
    return n if n % 2 else max(n - 1, 1)


def _interp_from(t: np.ndarray, values: np.ndarray, ok: np.ndarray, fallback: float) -> np.ndarray:
    """*values* at the *ok* samples, linearly interpolated in time elsewhere (held at the ends); fallback if none."""
    if ok.any():
        return np.interp(t, t[ok], values[ok])
    return np.full(len(t), fallback)


def _scan_state(y: np.ndarray, t: np.ndarray, window_s: float, rel_clean: float) -> dict:
    """Step 3 for one scan without cross-scan information.

    Returns a dict: kind (one of KINDS), base0 (K, per sample), rel (clip noise / level, per sample),
    trusted (bool per sample), valid (bool per sample), own_level (shorth of the valid samples, NaN if too few).
    """
    m = len(y)
    valid = valid_mask(y)
    yv = np.where(valid, y, np.nan)
    own_level, own_noise = location_noise(yv)
    state = dict(kind="own", base0=np.full(m, own_level), rel=np.full(m, rel_clean), trusted=np.zeros(m, bool),
                 valid=valid, own_level=own_level)
    if not valid.any():
        state["kind"] = "empty"
        return state
    if not np.isfinite(own_level):             # 1-4 valid samples: their median, strict noise
        state["own_level"] = float(np.median(y[valid]))
        state["base0"] = np.full(m, state["own_level"])
        return state

    w = window_samples(t, window_s, m)
    loc, noise = rolling_location_noise(yv, w)
    ratio = noise / loc
    trusted = np.isfinite(ratio) & (ratio <= TRUST_FACTOR * rel_clean)
    if trusted.any():
        rel = _interp_from(t, ratio, trusted, rel_clean)
        rel = uniform_filter1d(rel, size=max(3, w // 3), mode="nearest")
        state.update(kind="window", base0=_interp_from(t, loc, trusted, own_level), rel=np.maximum(rel, rel_clean),
                     trusted=trusted)
        return state
    own_rel = own_noise / own_level if np.isfinite(own_noise) else np.inf
    if own_rel <= TRUST_SCAN_FACTOR * rel_clean:
        state.update(kind="scan", rel=np.full(m, max(own_rel, rel_clean)), trusted=np.ones(m, bool))
    return state


def _support(y: np.ndarray, valid: np.ndarray, level: float, noise_k: float, sigma: float) -> int:
    """Number of valid samples within sigma * noise_k of *level*."""
    return int(np.sum(valid & (np.abs(y - level) <= sigma * noise_k)))


def _borrow_source_levels(y: np.ndarray, t: np.ndarray, block: np.ndarray, states: dict, segments: list,
                          rel_clean: float, sigma: float) -> None:
    """Step 3 "source": give scans of kind "own"/"empty" the time-interpolated level of same-source scans.

    Modifies *states* in place.  Two passes per segment:
        A. A scan with valid samples only takes a borrowed level if its samples support it at least
           SOURCE_SUPPORT_FRAC as well as its own shorth level (support = valid samples within the same band,
           sigma * clean noise at the borrowed level, around each level).  Otherwise it keeps its own level and
           strict noise (kind "own", ``kept_own``): its data say the source level does not apply.
        B. With the levels settled in A (a kept own level is a neighbour too), every remaining "own"/"empty" scan
           gets the level interpolated in time from the scans with a level (held at the segment ends).
    """
    for segment in segments:
        keys = [int(block[idx[0]]) for idx in segment]
        t_mid = np.array([float(np.mean(t[idx])) for idx in segment])
        for check_support in (True, False):                     # pass A, then pass B
            has_level = np.array([states[k]["kind"] in ("window", "scan") or states[k].get("kept_own", False)
                                  for k in keys])
            if not has_level.any():
                break
            levels = np.array([float(np.median(states[k]["base0"])) for k in keys])
            for i, k in enumerate(keys):
                st = states[k]
                if st["kind"] not in ("own", "empty") or st.get("kept_own", False):
                    continue
                if check_support and st["kind"] == "empty":
                    continue
                borrowed = float(np.interp(t_mid[i], t_mid[has_level], levels[has_level]))
                ys = y[segment[i]]
                if check_support:
                    band = rel_clean * borrowed
                    sup_borrowed = _support(ys, st["valid"], borrowed, band, sigma)
                    sup_own = _support(ys, st["valid"], st["own_level"], band, sigma)
                    st["kept_own"] = sup_borrowed < SOURCE_SUPPORT_FRAC * sup_own
                    continue
                st.update(kind="source", base0=np.full(len(ys), borrowed), rel=np.full(len(ys), rel_clean))


def _finish_scan(y: np.ndarray, t: np.ndarray, st: dict, sigma: float,
                 window_s: float) -> tuple[np.ndarray, np.ndarray]:
    """Steps 4-6 for one scan.  Returns (output values, kept-sample mask)."""
    m = len(y)
    if st["kind"] == "empty":
        return np.full(m, np.nan), np.zeros(m, bool)
    valid, trusted, base0 = st["valid"], st["trusted"], st["base0"]
    noise = st["rel"] * base0
    good = valid & (np.abs(y - base0) <= sigma * noise)

    base = base0
    if st["kind"] in ("window", "scan") and good.sum() >= MIN_WINDOW_SAMPLES:
        w = window_samples(t, window_s, m)
        yg = np.where(good, y, np.nan)
        windows = _windows(yg, w)
        med = _nanmedian_rows(windows)
        enough = np.sum(np.isfinite(windows), axis=1) >= max(MIN_WINDOW_SAMPLES, w // 4)
        base = _interp_from(t, med, np.isfinite(med) & enough & trusted, float(np.median(base0)))
        base = uniform_filter1d(base, size=max(3, w // 3), mode="nearest")

    resid = y - base
    local_scatter = _MAD_SCALE * rolling_median(np.where(good, np.abs(resid), np.nan), DETAIL_SAMPLES)
    with np.errstate(divide="ignore", invalid="ignore"):
        weight = np.clip((TRUST_FACTOR * noise / local_scatter) ** 2, 0.0, 1.0)
    weight = np.where(trusted, 1.0, np.where(np.isfinite(weight), weight, 0.0))
    return np.where(good, base + weight * resid, base), good


# ---------------------------------------------------------------------- #
# Public API
# ---------------------------------------------------------------------- #

def clean_relative_noise(y: np.ndarray, t: np.ndarray, block: np.ndarray, window_s: float) -> float:
    """Step 2: clean relative noise of one channel (CLEAN_PERCENTILE of window noise/location, >= REL_NOISE_FLOOR)."""
    ratios = []
    for b in np.unique(block):
        idx = np.where(block == b)[0]
        yv = np.where(valid_mask(y[idx]), y[idx], np.nan)
        if np.sum(np.isfinite(yv)) < MIN_WINDOW_SAMPLES:
            continue
        loc, noise = rolling_location_noise(yv, window_samples(t[idx], window_s, len(idx)))
        ratios.append(noise / loc)
    if not ratios:
        return REL_NOISE_FLOOR
    all_ratios = np.concatenate(ratios)
    all_ratios = all_ratios[np.isfinite(all_ratios)]
    if not len(all_ratios):
        return REL_NOISE_FLOOR
    return max(REL_NOISE_FLOOR, float(np.percentile(all_ratios, CLEAN_PERCENTILE)))


def clean_tsys_channel(values: np.ndarray, times: np.ndarray, block: np.ndarray, sigma: float = 3.0,
                       window_s: float = 0.0, segments: list[list[np.ndarray]] | None = None,
                       stats: dict | None = None) -> np.ndarray:
    """Clean one BBC-channel Tsys series: outliers -> baseline, real measurements kept (see module docstring).

    Parameters
    ----------
    values : 1-D array
        Tsys in K.  NaN, <= 0 and >= MAX_VALID_TSYS samples are missing (pass manual flags / fills as NaN).
    times : 1-D array
        Unix timestamps, same length, non-decreasing within each scan.
    block : 1-D array
        Scan number per sample.
    sigma : float
        Clip threshold in units of the local clean noise.
    window_s : float
        Local statistics window in seconds (<= 0: DEFAULT_WINDOW_S).
    segments : list or None
        Same-source scan groups from ``smoothing.source_segments`` (enables the "source" rule); None = per scan only.
    stats : dict or None
        If given, filled with ``rel_clean`` (float), ``kinds`` ({kind: n scans}) and ``n_clipped`` (valid samples
        replaced by the baseline).

    Returns
    -------
    ndarray
        Cleaned values, same length; NaN where no estimate exists (scans of kind "empty").
    """
    y = np.asarray(values, dtype=float)
    t = np.asarray(times, dtype=float)
    blk = np.asarray(block)
    window_s = window_s if window_s > 0 else DEFAULT_WINDOW_S
    out = np.full(len(y), np.nan)
    if len(y) == 0:
        return out

    rel_clean = clean_relative_noise(y, t, blk, window_s)
    scans = {int(b): np.where(blk == b)[0] for b in np.unique(blk)}
    states = {k: _scan_state(y[idx], t[idx], window_s, rel_clean) for k, idx in scans.items()}
    if segments is not None:
        _borrow_source_levels(y, t, blk, states, segments, rel_clean, sigma)

    kinds = {kind: 0 for kind in KINDS}
    n_clipped = 0
    for k, idx in scans.items():
        out[idx], good = _finish_scan(y[idx], t[idx], states[k], sigma, window_s)
        kinds[states[k]["kind"]] += 1
        n_clipped += int(np.sum(states[k]["valid"] & ~good))
    if stats is not None:
        stats.update(rel_clean=rel_clean, kinds=kinds, n_clipped=n_clipped)
    return out

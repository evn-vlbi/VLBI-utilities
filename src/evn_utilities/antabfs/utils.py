"""Shared utility functions for antabfs processing.

Contains Tcal lookup, data pre-filtering, and outlier detection helpers
used across the antabfs pipeline.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy import stats as _scipy_stats

_MAD_SCALE = 1.4826


def get_rxg_dirs() -> list[str]:
    """Return RXG search directories from the user config, or ['.'] as fallback.

    Reads ``$XDG_CONFIG_HOME/evn/evn_rxg.config`` (defaulting to
    ``~/.config/evn/evn_rxg.config``).  Each non-blank, non-comment line
    is treated as a directory path.  If the file does not exist or contains
    no valid paths, returns ``['.']`` (current working directory).
    """
    xdg_config = os.environ.get("XDG_CONFIG_HOME", os.path.expanduser("~/.config"))
    config_file = os.path.join(xdg_config, "evn", "evn_rxg.config")

    if os.path.isfile(config_file):
        dirs = []
        with open(config_file) as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    dirs.append(line)
        if dirs:
            return dirs

    return ["."]


def get_tcal(lo_freq: float, pol: str, freq: float, station: str,
             rxg_dirs: list[str] | None = None, rxg_files: list[str] | None = None,
             debug: bool = False) -> float:
    """Interpolate Tcal from RXG files for a given frequency and polarization.

    Parameters
    ----------
    lo_freq : float
        Local oscillator frequency in MHz.
    pol : str
        Polarization string (first 3 chars matched, e.g. ``'rcp'``).
    freq : float
        Sky frequency of the channel in MHz.
    station : str
        Two-letter station code (case-insensitive).
    rxg_dirs : list[str] or None
        Directories to search for RXG files.  If ``None``, uses
        :func:`get_rxg_dirs` to determine the search path.
    rxg_files : list[str] or None
        Explicit list of RXG filenames.  Searched in each of *rxg_dirs*.
        If ``None``, all ``.rxg`` files in the search dirs are scanned.
    debug : bool
        Print diagnostic messages.

    Returns
    -------
    float
        Interpolated Tcal value (0.0 if not found).
    """
    if rxg_dirs is None:
        rxg_dirs = get_rxg_dirs()

    rxg_list: list[str] = []
    for cal_dir in rxg_dirs:
        cal_dir = cal_dir.rstrip("/") + "/"
        if rxg_files:
            for f in rxg_files:
                candidate = cal_dir + f
                if os.path.isfile(candidate):
                    rxg_list.append(candidate)
        else:
            try:
                rxg_list.extend(
                    cal_dir + f
                    for f in os.listdir(cal_dir)
                    if f.endswith(".rxg")
                )
            except FileNotFoundError:
                continue

    tcal = 0.0
    file_ok = False
    st_code = station[0].upper() + station[1].lower()

    for filename in rxg_list:
        if st_code not in filename and not rxg_files:
            continue

        with open(filename) as fh:
            lines = fh.read().splitlines()

        file_ok = False
        for i, line in enumerate(lines):
            if line[:5] == "range":
                parts = line.split()
                rmin, rmax = float(parts[1]), float(parts[2])
                if rmin <= lo_freq <= rmax:
                    file_ok = True
                    if debug:
                        print(f"Using {filename} RXG file")
                else:
                    break
            elif line[:5] == "fixed":
                fval = float(line.split()[1])
                rmin, rmax = fval - 10, fval + 10
                if rmin <= lo_freq <= rmax:
                    file_ok = True
                    if debug:
                        print(f"Using {filename} RXG file")
                else:
                    break

            if file_ok and i < len(lines) - 1:
                if lines[i][:3] == pol and lines[i + 1][:3] == pol:
                    f1 = float(lines[i].split()[1])
                    f2 = float(lines[i + 1].split()[1])
                    if f1 <= freq <= f2:
                        t1 = float(lines[i].split()[2])
                        t2 = float(lines[i + 1].split()[2])
                        tcal = t1 + (freq - f1) * (t2 - t1) / (f2 - f1)
                        break

        if tcal != 0:
            break

    if not file_ok:
        print(f"tcal = {tcal:g}")
        print("A suitable rxg_file was not found. "
              'Maybe tcal is inside LOG file ("caltemp" tag)')

    return tcal


def prefilter(tsys_matrix: np.ndarray, block: np.ndarray, max_limit: float = 10000.0) -> np.ndarray:
    """Replace negative or excessively large Tsys values with geometric means.

    For each BBC channel (row), bad values are replaced by the geometric mean
    of good values within the same scan block.  The per-block geometric mean
    is computed in one vectorised pass using ``np.bincount``.  Channels are
    processed in parallel via a thread pool.

    Parameters
    ----------
    tsys_matrix : ndarray, shape (n_channels, n_times)
        Transposed Tsys array (one row per BBC channel).
    block : ndarray, shape (n_times,)
        Scan block index for each time sample.
    max_limit : float
        Upper threshold for valid Tsys.

    Returns
    -------
    ndarray
        Filtered copy of *tsys_matrix*.
    """
    tsys = np.array(tsys_matrix, dtype=float)
    block = np.asarray(block)

    # Pre-compute block inverse index once; reused for every channel.
    _, inv = np.unique(block, return_inverse=True)
    n_blks = int(inv.max()) + 1

    def _filter_row(row: np.ndarray) -> None:
        """Replace bad values in *row* in-place (row is a view of tsys)."""
        valid = (row > 0) & (row < max_limit)
        if valid.all():
            return

        # Log-space values; NaN where invalid so nanmean ignores them.
        with np.errstate(divide="ignore", invalid="ignore"):
            log_vals = np.where(valid, np.log(row), np.nan)

        # Per-block sum of log-values and count of valid values via bincount.
        log_sum = np.bincount(inv, weights=np.where(valid, log_vals, 0.0), minlength=n_blks)
        cnt = np.bincount(inv, weights=valid.astype(np.float64), minlength=n_blks)
        with np.errstate(invalid="ignore"):
            blk_log_gm = np.where(cnt > 0, log_sum / cnt, np.nan)

        bad_pos = np.where(~valid)[0]
        bad_blk = inv[bad_pos]
        gm_log = blk_log_gm[bad_blk].copy()  # per-bad-point log geomean

        # Handle blocks where every value is bad (NaN in gm_log).
        nan_mask = np.isnan(gm_log)
        if nan_mask.any():
            global_log = float(np.nanmean(log_vals)) if valid.any() else 0.0
            for k in np.where(nan_mask)[0]:
                bi = int(bad_blk[k])
                gm_log[k] = next(
                    (blk_log_gm[nb] for nb in range(bi + 1, n_blks)
                     if not np.isnan(blk_log_gm[nb])),
                    global_log,
                )

        row[bad_pos] = np.exp(gm_log)

    # Rows of a C-contiguous 2-D array are non-overlapping memory regions;
    # modifying different rows from different threads is safe.
    n_ch = tsys.shape[0]
    max_workers = min(n_ch, os.cpu_count() or 1)
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        list(pool.map(_filter_row, tsys))   # tsys rows are views

    return tsys


def sm_fit(x: np.ndarray, y: np.ndarray, sigma: float = 3.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Robust Theil-Sen linear fit with MAD-based prediction bounds.

    Theil-Sen uses the median of all pairwise slopes and is therefore
    insensitive to even a large fraction of high outliers.  Bounds are
    derived from the MAD of the fit residuals rather than a fixed
    percentage, so they adapt to the actual scatter in the data.

    Returns
    -------
    fitted, lower, upper
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if len(x) > 1:
        xc = x - x.mean()
        res = _scipy_stats.theilslopes(y, xc)
        fitted = res.slope * xc + res.intercept
        residuals = y - fitted
        resid_med = float(np.median(residuals))
        scale = _MAD_SCALE * float(np.median(np.abs(residuals - resid_med)))
        if scale == 0.0:
            scale = max(1.0, 0.01 * float(np.median(np.abs(fitted))))
        return fitted, fitted - sigma * scale, fitted + sigma * scale
    if len(x) == 1:
        return np.array(y), np.array(y * 0.9), np.array(y * 1.1)
    return np.array([]), np.array([]), np.array([])


def compute_outliers(block: np.ndarray, x: np.ndarray, y: np.ndarray,
                     sigma: float = 3.0) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                                                   np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Identify outliers by fitting per-block robust linear models.

    Uses a two-pass strategy so that small scan blocks (1–3 points) with
    embedded outliers do not distort the fit or the bounds:

    Pass 1 — per-block median as a rough fit (vectorised via
    ``np.unique(return_inverse=True)``), global pooled MAD to exclude
    obvious outliers before fitting.

    Pass 2 — per-block Theil-Sen regression on the cleaned data, evaluated
    at *all* block positions.  Bounds use the global MAD of clean residuals.

    Returns
    -------
    fit, low, up, in_x, in_y, out_x, out_y
    """
    block = np.asarray(block)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # Inverse index maps every sample to its 0-based block position.
    _, inv = np.unique(block, return_inverse=True)
    n_blks = int(inv.max()) + 1

    # ------------------------------------------------------------------
    # Pass 1: per-block median → rough_fit (fully vectorised).
    # ------------------------------------------------------------------
    block_meds = np.array([np.median(y[inv == i]) for i in range(n_blks)])
    rough_fit = block_meds[inv]

    rough_resid = y - rough_fit
    r_med = float(np.median(rough_resid))
    global_scale = _MAD_SCALE * float(np.median(np.abs(rough_resid - r_med)))
    if global_scale == 0.0:
        global_scale = max(1.0, 0.01 * float(np.median(np.abs(rough_fit))))

    pre_inlier = np.abs(rough_resid) <= 4.0 * global_scale

    # ------------------------------------------------------------------
    # Pass 2: Theil-Sen per block on pre-filtered clean data.
    # Results written directly into pre-allocated fit_arr.
    # ------------------------------------------------------------------
    fit_arr = np.empty_like(y)

    for bi in range(n_blks):
        mask = inv == bi
        idx = np.where(mask)[0]
        clean = mask & pre_inlier
        n_clean = int(clean.sum())

        if n_clean >= 2:
            xb, yb = x[clean], y[clean]
            xc_mean = xb.mean()
            res = _scipy_stats.theilslopes(yb, xb - xc_mean)
            fit_arr[idx] = res.slope * (x[idx] - xc_mean) + res.intercept
        elif n_clean == 1:
            fit_arr[idx] = float(y[clean][0])
        else:
            fit_arr[idx] = float(rough_fit[idx[0]])

    # ------------------------------------------------------------------
    # Final bounds: global MAD of clean residuals from the pass-2 fit.
    # ------------------------------------------------------------------
    final_resid = y - fit_arr
    clean_resid = final_resid[pre_inlier]
    if len(clean_resid) >= 3:
        cr_med = float(np.median(clean_resid))
        final_scale = _MAD_SCALE * float(np.median(np.abs(clean_resid - cr_med)))
        if final_scale == 0.0:
            final_scale = global_scale
    else:
        final_scale = global_scale

    low = fit_arr - sigma * final_scale
    up  = fit_arr + sigma * final_scale
    in_mask = (y >= low) & (y <= up)

    return fit_arr, low, up, x[in_mask], y[in_mask], x[~in_mask], y[~in_mask]


def modify_data(x: np.ndarray, y: np.ndarray, block: np.ndarray,
                xmax: float, xmin: float, ymax: float, ymin: float) -> np.ndarray:
    """Replace points inside a selection rectangle with geometric-mean estimates."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    block = np.asarray(block)

    keep = (x < xmin) | (x > xmax) | (y < ymin) | (y > ymax) | (y < 0)
    selected = np.where(~keep)[0]
    if not len(selected):
        return y

    _, inv = np.unique(block, return_inverse=True)
    n_blks = int(inv.max()) + 1

    for i in selected:
        bi = int(inv[i])
        # Good candidates: kept values in the same block.
        cand_idx = np.where(keep & (inv == bi))[0]
        if not len(cand_idx):
            # Search forward through neighbouring blocks.
            for delta in range(1, min(21, n_blks - bi)):
                cand_idx = np.where(keep & (inv == bi + delta))[0]
                if len(cand_idx):
                    break
            if not len(cand_idx):
                cand_idx = np.where(keep)[0]
        if len(cand_idx):
            y[i] = np.exp(np.mean(np.log(np.abs(y[cand_idx]))))

    return y

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
        Tcal linearly interpolated from the first matching RXG file's table; outside the table the closest
        (end) value is used. 0.0 if no RXG file covers *lo_freq* or it has no entries for *pol*.
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

    file_ok = False
    st_code = station[0].upper() + station[1].lower()

    for filename in rxg_list:
        if st_code not in filename and not rxg_files:
            continue
        table = read_tcal_table(filename, lo_freq, pol)
        if table is None:
            continue
        file_ok = True
        if debug:
            print(f"Using {filename} RXG file")
        if not table:
            continue
        freqs, tcals = zip(*sorted(table))
        if debug and not freqs[0] <= freq <= freqs[-1]:
            print(f"  freq={freq:.2f} MHz ({pol}) outside the Tcal table [{freqs[0]:.2f}, {freqs[-1]:.2f}] MHz "
                  f"of {filename}: using the closest table value.")
        # np.interp clamps to the end values outside the table, i.e. the closest frequency.
        return float(np.interp(freq, freqs, tcals))

    # Missing RXG files are reported once (as a fatal error) by the CLI, not per channel here.
    if debug:
        reason = "no Tcal entries for this polarization" if file_ok else "no RXG file found"
        print(f"{reason}: LO={lo_freq:.2f} MHz, pol={pol}, freq={freq:.2f} MHz (tcal=0).")
    return 0.0


def read_tcal_table(filename: str, lo_freq: float, pol: str) -> list[tuple[float, float]] | None:
    """Read the ``(freq_mhz, tcal_k)`` Tcal-table entries of polarization *pol* from an RXG file.

    Returns ``None`` if the file's LO ``range`` (or ``fixed`` LO +/- 10 MHz) does not contain *lo_freq*,
    otherwise the (possibly empty) list of entries whose line starts with *pol* (e.g. ``'rcp'``).
    """
    with open(filename) as fh:
        lines = fh.read().splitlines()

    lo_ok: bool | None = None
    table: list[tuple[float, float]] = []
    for line in lines:
        parts = line.split()
        if lo_ok is None and line[:5] == "range":
            lo_ok = float(parts[1]) <= lo_freq <= float(parts[2])
        elif lo_ok is None and line[:5] == "fixed":
            lo_ok = abs(float(parts[1]) - lo_freq) <= 10
        elif line.startswith("end_tcal_table"):
            break
        elif lo_ok and line[:3] == pol and len(parts) >= 3:
            table.append((float(parts[1]), float(parts[2])))
        if lo_ok is False:
            return None
    return table if lo_ok else None


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

"""ANTAB file writer.

Writes the Tsys data section of an ANTAB file, including scan markers,
DPFU/POLY blocks, and inline Tsys-from-log comments.
"""

from __future__ import annotations

import datetime
import io
from functools import lru_cache

import numpy as np


@lru_cache(maxsize=None)
def _format_time(unix_ts: float) -> str:
    """Format a Unix timestamp as ``'DDD HH:MM.mm'``.

    Result is cached with ``lru_cache`` so repeated calls for the same
    timestamp (duplicate-check vs. write path) cost only a dict lookup.
    Rounds the fractional minute to 2 decimal places and carries to the
    next minute / hour / day on overflow.
    """
    dt = datetime.datetime.fromtimestamp(unix_ts, tz=datetime.timezone.utc)
    d = dt.timetuple().tm_yday
    h = dt.hour
    m = round(dt.minute + dt.second / 60.0 + dt.microsecond / 60e6, 2)

    if m >= 60.0:
        m -= 60.0
        h += 1
        if h >= 24:
            h -= 24
            d += 1

    return f"{d:03d} {h:02d}:{m:05.2f}"


def _flush_tsys_log(buf: io.StringIO, tsys_log: list[list[float]], idx: int,
                    time_tsys_log: float, time_limit: float, tsys_line: np.ndarray,
                    data_idx: int, seen: set[str]) -> tuple[int, float]:
    """Write Tsys-from-log comment entries up to *time_limit* (Unix seconds).

    Replaces the former ``_write_tsys_log_before`` / ``_write_tsys_log_until``
    pair.  Both had identical logic; the old ``_before`` variant converted the
    Unix timestamp through a formatted string back to a fractional day just to
    do a numeric comparison — this function compares Unix timestamps directly,
    eliminating that round-trip.
    """
    while idx < len(tsys_log) and time_tsys_log <= time_limit:
        ts_key = _format_time(time_tsys_log)
        if ts_key not in seen:
            seen.add(ts_key)
            entry = tsys_log[idx]
            vals = np.asarray(entry[1:], dtype=float)
            ref = np.asarray(tsys_line[data_idx], dtype=float)
            m = min(len(vals), len(ref))
            buf.write(
                f"\n! {ts_key}"
                + "".join(f" {v:.1f}" for v in vals[:m][np.isfinite(ref[:m])])
            )
        idx += 1
        if idx < len(tsys_log):
            time_tsys_log = tsys_log[idx][0]

    return idx, time_tsys_log


def write_antab(file_out: str, header: list[str], index_line: list[str], scan_line: list[str],
                tsys_line: np.ndarray, block: list[int], time: list[float],
                tsys_log: list[list[float]], setup_time: list[list],
                dpfu_lines: dict[str, list[str]], polyelev_line: dict[str, list[str]],
                station_name: str) -> None:
    """Append Tsys data to an ANTAB file.

    All output is accumulated in an in-memory ``io.StringIO`` buffer and
    written to disk in a single ``f.write()`` call, minimising system-call
    overhead.  Additional optimisations:

    * ``scan_dict`` precomputes the block → scan-comment mapping for O(1)
      lookup instead of an O(n_scans) search at each block boundary.
    * Tsys-from-log flushing uses direct Unix-timestamp comparison, removing
      the string-round-trip in the former ``_write_tsys_log_before``.
    * ``np.isfinite`` replaces ``try: int(j_val)`` for valid-value filtering.
    * ``lru_cache`` on ``_format_time`` amortises repeated timestamp formatting.

    Parameters
    ----------
    file_out : str
        Path to the output ANTAB file (opened in append mode).
    header : list[str]
        Per-setup header comment blocks.
    index_line : list[str]
        Per-setup INDEX= lines.
    scan_line : list[str]
        Scan comment lines with source/time info.
    tsys_line : ndarray
        2-D array of Tsys values, shape ``(n_times, n_channels)``.
    block : list[int]
        Scan block index per time sample.
    time : list[float]
        UNIX timestamps per sample.
    tsys_log : list[list[float]]
        Tsys read directly from the log file (first element is timestamp).
    setup_time : list[list]
        ``[[timestamp, setup_name], ...]`` boundaries.
    dpfu_lines : dict
        DPFU line strings keyed by setup.
    polyelev_line : dict
        POLY line strings keyed by setup.
    station_name : str
        Two-letter station code.
    """
    # Precompute scan dict: block_num → scan comment string  O(n_scans) once
    scan_dict: dict[int, str] = {}
    for scan in scan_line:
        try:
            scan_dict[int(scan.split("=")[1].split(" ")[0])] = scan
        except (IndexError, ValueError):
            pass

    buf = io.StringIO()
    tsys_log_idx = 0
    setup_idx = 0
    time_tsys_log = tsys_log[0][0] if tsys_log else 0.0
    seen: set[str] = set()

    for i, (blk_i, t_i) in enumerate(zip(block, time)):
        # --- Setup transitions ---
        if setup_idx < len(setup_time):
            setup = setup_time[setup_idx][1]
            if setup_time[setup_idx][0] <= t_i:
                # Cache the per-setup line lists; avoids repeated .get() calls
                setup_dpfu = dpfu_lines.get(setup, [])
                setup_poly = polyelev_line.get(setup, [])
                for dl, pl in zip(setup_dpfu, setup_poly):
                    if setup_idx > 0:
                        buf.write("\n/\n")
                    buf.write(dl + " " + pl)
                buf.write("/\n")
                buf.write(f"TSYS {station_name} FT = 1.0 TIMEOFF=0\n")
                buf.write(index_line[setup_idx][:-1] + "\n")
                buf.write("/\n")
                for h_line in header[setup_idx]:
                    buf.write(h_line)
                setup_idx += 1

        # --- Scan marker: flush tsys_log up to this data point, then write ---
        if i == 0 or blk_i != block[i - 1]:
            if blk_i in scan_dict:
                # Flush log entries that precede this scan (direct Unix compare)
                tsys_log_idx, time_tsys_log = _flush_tsys_log(buf, tsys_log, tsys_log_idx,
                    time_tsys_log, t_i, tsys_line, i, seen)
                buf.write(scan_dict[blk_i])

        # --- Tsys-from-log interleaved up to current data timestamp ---
        if tsys_log_idx < len(tsys_log):
            tsys_log_idx, time_tsys_log = _flush_tsys_log(buf, tsys_log, tsys_log_idx,
                time_tsys_log, t_i, tsys_line, i, seen)

        # --- Main Tsys line ---
        ts_key = _format_time(t_i)
        if ts_key in seen:
            continue
        seen.add(ts_key)

        row = tsys_line[i]
        buf.write(f"\n{ts_key}" + "".join(f" {v:.1f}" for v in row[np.isfinite(row)]))

    # --- Remaining Tsys-from-log entries after all data ---
    while tsys_log_idx < len(tsys_log):
        ts_key = _format_time(time_tsys_log)
        if ts_key not in seen:
            seen.add(ts_key)
            buf.write(f"\n! {ts_key}"
                      + "".join(f" {v:.1f}" for v in tsys_log[tsys_log_idx][1:]))
        tsys_log_idx += 1
        if tsys_log_idx >= len(tsys_log):
            break
        time_tsys_log = tsys_log[tsys_log_idx][0]

    buf.write("\n/\n")

    with open(file_out, "a") as f:
        f.write(buf.getvalue())

"""Command-line interface for antabfs.

Replaces the old ``sys.argv`` parsing with ``argparse`` and orchestrates
the full ANTAB-file generation pipeline.
"""

from __future__ import annotations

import argparse
import bisect
import os
import subprocess
import sys

import numpy as np

from . import __version__ as version
from .header import AntabHeader
from .logfile import LogFile
from .plotting import MultiChannelSelection
from .smoothing import smooth_tsys_matrix
from .utils import compute_outliers, prefilter
from .writer import write_antab


def build_parser() -> argparse.ArgumentParser:
    """Build and return the argument parser."""
    try:
        from rich_argparse import RichHelpFormatter as Formatter
    except ImportError:
        Formatter = argparse.HelpFormatter

    parser = argparse.ArgumentParser(prog="antabfs", formatter_class=Formatter,
        description="Generate ANTAB-format amplitude calibration files from Field System log files.")
    parser.add_argument("logfile", help="Path to the Field System log file.")
    parser.add_argument("-r", "--rxg-dir", default=None,
        help="Directory containing RXG calibration files. Overrides any config-file paths. "
             "If not set, directories are read from $XDG_CONFIG_HOME/evn/evn_rxg.config "
             "(or ~/.config/evn/evn_rxg.config), falling back to the current directory.")
    parser.add_argument("-f", "--rxg-files", nargs="+", default=None,
        help="Explicit list of RXG filenames to use.")
    parser.add_argument("-o", "--output", default=None,
        help="Output ANTAB filename (default: <experiment><station>.antabfs).")
    parser.add_argument("-i", "--interactive", action="store_true", default=False,
        help="Enable interactive (matplotlib) outlier editing.")
    parser.add_argument("--no-smooth", action="store_true", default=False,
        help="Disable automatic Tsys smoothing (smoothing is ON by default).")
    parser.add_argument("--sigma", type=float, default=3.0,
        help="Sigma-clipping threshold for smoothing (default: %(default)s).")
    parser.add_argument("--smooth-window", type=int, default=0,
        help="Smoothing window size (0 = auto, default: %(default)s).")
    parser.add_argument("--smooth-method", choices=["savgol", "median"], default="savgol",
        help="Smoothing method: savgol or median (default: %(default)s).")
    parser.add_argument("--vlbeer", action="store_true", default=False,
        help="Upload the output ANTAB file to vlbeer.ira.inaf.it (as user 'evn') "
             "into vlbi_arch/mmmYY/ where mmmYY is the observation month and year.")
    parser.add_argument("-d", "--debug", action="store_true", default=False,
        help="Print debug messages.")
    parser.add_argument("-V", "--version", action="version", version=f"%(prog)s {version}")
    return parser


def main(argv: list[str] | None = None) -> None:
    """Main entry point for antabfs."""
    parser = build_parser()
    args = parser.parse_args(argv)

    log_path = args.logfile
    if not os.path.isfile(log_path):
        print(f"Error: log file '{log_path}' not found.", file=sys.stderr)
        sys.exit(1)

    # --- Parse the log file ---
    print(f"Reading log file: {log_path}")
    log_f = LogFile(log_path, rxg_dir=args.rxg_dir, rxg_files=args.rxg_files, debug=args.debug)

    header, index_line, scan_line, tsys_data, block_data, time_data, tsys_log, setup_time = \
        log_f.get_log_data()

    if not tsys_data:
        print("No Tsys data found in the log file.", file=sys.stderr)
        sys.exit(1)

    # --- Build output filename ---
    antab_file = args.output or f"{log_f.experiment()}{log_f.station().lower()}.antabfs"

    # --- Build ANTAB header ---
    antab_header = AntabHeader(log_f)
    dpfu_lines, polyelev_lines = antab_header.write_antab_preamble(antab_file)

    # --- Prepare Tsys matrix ---
    tsys_t = prefilter(np.array(tsys_data, dtype=float).T, np.array(block_data))
    block_array = np.array(block_data)
    time_array = np.array(time_data)

    # --- Determine BBC labels from header ---
    bbc_labels: list[str] = [
        h_line.split("=")[1].split(":")[0].strip()
        for h in header for h_line in h.split("\n")
        if "Column" in h_line and len(h_line.split("=")) > 1
    ]

    # --- Process each setup/part ---
    part_boundaries = [(0, len(block_data))]
    if len(setup_time) > 1:
        part_boundaries = []
        for bp_i in range(len(setup_time)):
            t_start = setup_time[bp_i][0]
            t_end = setup_time[bp_i + 1][0] if bp_i + 1 < len(setup_time) else time_data[-1] + 1
            idx_start = _find_index(time_data, t_start)
            idx_end = _find_index(time_data, t_end)
            if idx_start is not None and idx_end is not None:
                part_boundaries.append((idx_start, idx_end))
            elif idx_start is not None:
                part_boundaries.append((idx_start, len(time_data)))

    for part_num, (i_start, i_end) in enumerate(part_boundaries):
        if i_start >= i_end:
            continue

        n_channels = tsys_t.shape[0]
        block_part = block_array[i_start:i_end]
        time_part = time_array[i_start:i_end]
        tsys_part = tsys_t[:, i_start:i_end]

        if args.interactive:
            # MultiChannelSelection auto-replaces outliers and writes back in-place.
            MultiChannelSelection(tsys_t[:, i_start:i_end], time_part, block_part,
                bbc_labels[:n_channels], sigma=args.sigma, part_num=part_num)
        elif not args.no_smooth:
            print(f"  Smoothing part {part_num + 1} ({n_channels} channels, {i_end - i_start} samples)...")
            tsys_t[:, i_start:i_end] = smooth_tsys_matrix(tsys_part, times=time_part,
                block=block_part, sigma=args.sigma, window=args.smooth_window, method=args.smooth_method)
        else:
            # No smoothing: detect outliers and replace with local median.
            for ch_i in range(n_channels):
                y = tsys_part[ch_i].copy()
                n = len(y)
                _, low, up, _, _, _, _ = compute_outliers(block_part, time_part, y, sigma=args.sigma)
                out_mask = (y > up) | (y < low)
                good = ~out_mask
                if out_mask.any() and good.any():
                    global_med = float(np.nanmedian(y[good]))
                    for bi in np.where(out_mask)[0]:
                        lo, hi = max(0, bi - 20), min(n, bi + 21)
                        local_good = good[lo:hi]
                        y[bi] = float(np.nanmedian(y[lo:hi][local_good])) if local_good.any() else global_med
                tsys_t[ch_i, i_start:i_end] = y

    write_antab(antab_file, header, index_line, scan_line, tsys_t.T, block_data, time_data,
        tsys_log, setup_time, dpfu_lines, polyelev_lines, log_f.station())
    print(f"ANTAB file written: {antab_file}")

    if args.vlbeer:
        _upload_to_vlbeer(antab_file, log_f)


def _upload_to_vlbeer(antab_file: str, log_f: LogFile) -> None:
    """Upload *antab_file* to vlbeer.ira.inaf.it under vlbi_arch/mmmYY/."""
    obs_date = log_f.observation_date()
    if obs_date is None:
        print("Could not determine observation date from log file; skipping upload.", file=sys.stderr)
        return

    mmm_yy = obs_date.strftime("%b%y").lower()   # e.g. "feb26"
    remote_dir = f"vlbi_arch/{mmm_yy}"
    host = "evn@vlbeer.ira.inaf.it"

    print(f"Creating remote directory {remote_dir} on vlbeer...")
    subprocess.run(["ssh", host, f"mkdir -p {remote_dir}"], check=False)

    print(f"Uploading {antab_file} to {host}:{remote_dir}/...")
    result = subprocess.run(["scp", antab_file, f"{host}:{remote_dir}/"])
    if result.returncode == 0:
        print(f"Upload successful: {host}:{remote_dir}/{os.path.basename(antab_file)}")
    else:
        print(f"Upload failed (scp exit code {result.returncode}).", file=sys.stderr)


def _find_index(time_list: list[float], target: float) -> int | None:
    """Return index of first element >= target, or None (bisect O(log n))."""
    i = bisect.bisect_left(time_list, target)
    return i if i < len(time_list) else None


if __name__ == "__main__":
    main()

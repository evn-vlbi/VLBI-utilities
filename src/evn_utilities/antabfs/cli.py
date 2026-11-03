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
from .header import AntabHeader, parse_sefd, sefd_gain_lines
from .logfile import LogFile
from .plotting import MultiChannelSelection
from .smoothing import SMOOTH_SPANS, smooth_tsys_matrix
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
    parser.add_argument("-f", "--rxg-files", default=None,
        help="Comma-separated list of RXG filenames to use (e.g. -f Ys_6cm.rxg,Ys_13cm.rxg).")
    parser.add_argument("--sefd", default=None, metavar="JY[,JY]",
        help="Nominal mode: SEFD (Jy) of the antenna at the observed band, either one value for both "
             "polarizations or RCP,LCP (e.g. --sefd 350 or --sefd 340,360). No RXG file is searched or read and "
             "/caltemp/ lines are ignored: Tsys is computed with Tcal = 1 for every channel, and DPFU is set to "
             "median(Tsys)/SEFD with no gain curve (POLY=1.0). Cannot be combined with -r/-f.")
    parser.add_argument("-o", "--output", default=None,
        help="Output ANTAB filename (default: <experiment><station>.antabfs).")
    parser.add_argument("-i", "--interactive", action="store_true", default=False,
        help="Open the interactive multi-channel editor to flag outliers and choose the smoothing method/window "
             "(initialised from --smooth-method, --smooth-window and --no-smooth) before the file is written. "
             "By default antabfs runs non-interactively (auto-smooth, no window).")
    parser.add_argument("--no-smooth", action="store_true", default=False,
        help="Disable Tsys trend smoothing (smoothing is ON by default). Outliers are still "
             "rejected per scan and replaced with the local median.")
    parser.add_argument("--plot", action="store_true", default=False,
        help="Show a read-only plot of the smoothed Tsys result before writing the file.")
    parser.add_argument("--sigma", type=float, default=3.0,
        help="Sigma-clipping threshold for smoothing (default: %(default)s).")
    parser.add_argument("--smooth-span", choices=list(SMOOTH_SPANS), default="source",
        help="source: smooth Tsys across all scans of the same source (never mixing sources) with a smoothing "
             "spline through the per-scan robust medians; scan: smooth each scan on its own with "
             "--smooth-method/--smooth-window (default: %(default)s).")
    parser.add_argument("--source-gap", type=float, default=60.0, metavar="MIN",
        help="With --smooth-span source: scans of the same source more than MIN minutes apart are smoothed "
             "as separate segments (default: %(default)s).")
    parser.add_argument("--smooth-window", type=int, default=0,
        help="Per-scan smoothing window size, --smooth-span scan with median only (0 = the full length of "
             "each scan, default: %(default)s).")
    parser.add_argument("--smooth-method", choices=["savgol", "median"], default="median",
        help="Per-scan smoothing method, --smooth-span scan only: savgol or median (default: %(default)s).")
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

    sefd = None
    if args.sefd is not None:
        if args.rxg_dir or args.rxg_files:
            parser.error("--sefd cannot be combined with -r/--rxg-dir or -f/--rxg-files (no RXG file is used).")
        try:
            sefd = parse_sefd(args.sefd)
        except ValueError as ex:
            parser.error(f"invalid --sefd value: {ex}")

    log_path = args.logfile
    if not os.path.isfile(log_path):
        _fatal(f"log file '{log_path}' not found.")

    # --- Parse the log file ---
    print(f"Reading log file: {log_path}")
    if sefd:
        print(f"SEFD mode: SEFD(R) = {sefd['R']:g} Jy, SEFD(L) = {sefd['L']:g} Jy. "
              "RXG files and /caltemp/ are not used (Tcal = 1).")
    rxg_files = args.rxg_files.split(",") if args.rxg_files else None
    log_f = LogFile(log_path, rxg_dir=args.rxg_dir, rxg_files=rxg_files, debug=args.debug, use_rxg=sefd is None)

    if sefd is None:
        _check_rxg(log_f, rxg_files, user_rxg=bool(args.rxg_dir or args.rxg_files))

    header, index_line, scan_line, tsys_data, block_data, time_data, tsys_log, setup_time = \
        log_f.get_log_data()

    if not tsys_data:
        _fatal("no Tsys data found in the log file.")

    # --- Build output filename ---
    antab_file = args.output or f"{log_f.experiment()}{log_f.station().lower()}.antabfs"
    antab_header = AntabHeader(log_f, sefd=sefd)

    # --- Prepare Tsys matrix ---
    tsys_t = prefilter(np.array(tsys_data, dtype=float).T, np.array(block_data))
    block_array = np.array(block_data)
    time_array = np.array(time_data)

    # --- Process each setup/part: (start index, end index, index into setup_time/header) ---
    part_boundaries = [(0, len(block_data), 0)]
    if len(setup_time) > 1:
        part_boundaries = []
        for bp_i in range(len(setup_time)):
            t_start = setup_time[bp_i][0]
            t_end = setup_time[bp_i + 1][0] if bp_i + 1 < len(setup_time) else time_data[-1] + 1
            idx_start = _find_index(time_data, t_start)
            idx_end = _find_index(time_data, t_end)
            if idx_start is not None and idx_end is not None:
                part_boundaries.append((idx_start, idx_end, bp_i))
            elif idx_start is not None:
                part_boundaries.append((idx_start, len(time_data), bp_i))

    # Smoothing settings; the interactive editor may change them, and the choice carries over to the next part.
    smooth_method = "none" if args.no_smooth else args.smooth_method
    smooth_window = args.smooth_window
    sigma = args.sigma
    smooth_span = args.smooth_span
    if args.source_gap <= 0:
        parser.error("--source-gap must be > 0 minutes.")
    max_gap = args.source_gap * 60.0
    scan_sources = log_f.scan_sources()
    if smooth_span == "source":
        n_scans = len(set(block_data))
        n_known = sum(1 for b in set(block_data) if b in scan_sources)
        print(f"  Sources: {len(set(scan_sources.values()))} in the log; {n_known}/{n_scans} scans with Tsys "
              f"have a source (scans without one are smoothed on their own).")

    for part_num, (i_start, i_end, setup_idx) in enumerate(part_boundaries):
        if i_start >= i_end:
            continue

        n_channels = tsys_t.shape[0]
        block_part = block_array[i_start:i_end]
        time_part = time_array[i_start:i_end]

        # Optional interactive editing: the multi-channel editor lets the user
        # flag outliers and replaces them in-place inside tsys_t.  Off by default.
        if args.interactive:
            setup_name = setup_time[setup_idx][1] if setup_idx < len(setup_time) else ""
            labels = _channel_labels(header[setup_idx]) if setup_idx < len(header) else []
            labels += [f"ch{i + 1}" for i in range(len(labels), n_channels)]
            try:
                editor = MultiChannelSelection(tsys_t[:, i_start:i_end], time_part, block_part, labels[:n_channels],
                    sigma=sigma, part_num=part_num, smooth_method=smooth_method, smooth_window=smooth_window,
                    smooth_span=smooth_span, scan_sources=scan_sources, max_gap=max_gap,
                    calib_info=_calib_info(log_f, sefd, setup_name))
            except KeyboardInterrupt:
                _fatal("aborted by user (Ctrl+C). No ANTAB file was written.")
            smooth_method, smooth_window, sigma = editor.smooth_method, editor.smooth_window, editor.sigma
            smooth_span = editor.smooth_span
            print(f"  Editor closed (part {part_num + 1}): smooth method = {smooth_method}, span = {smooth_span}, "
                  f"window = {smooth_window}, sigma = {sigma:g}")

        # Smoothing pass (ON by default): rejects outliers and fits a smooth Tsys trend, per scan or across the
        # scans of each source.  Runs after interactive editing too, so the written file always follows a clean
        # trend instead of raw per-scan scatter.  Parts (setups) are smoothed separately: sources never span setups.
        if smooth_method != "none":
            print(f"  Smoothing part {part_num + 1} ({n_channels} channels, {i_end - i_start} samples, "
                  f"span = {smooth_span})...")
            tsys_t[:, i_start:i_end] = smooth_tsys_matrix(tsys_t[:, i_start:i_end], times=time_part,
                block=block_part, sigma=sigma, window=smooth_window, method=smooth_method, span=smooth_span,
                scan_sources=scan_sources, max_gap=max_gap, report=True,
                labels=_channel_labels(header[setup_idx]) if setup_idx < len(header) else None)
        elif not args.interactive:
            # --no-smooth without the editor: still reject obvious outliers and
            # replace them with the local median so the file isn't full of spikes.
            for ch_i in range(n_channels):
                y = tsys_t[ch_i, i_start:i_end].copy()
                n = len(y)
                _, low, up, _, _, _, _ = compute_outliers(block_part, time_part, y, sigma=sigma)
                out_mask = (y > up) | (y < low)
                good = ~out_mask
                if out_mask.any() and good.any():
                    global_med = float(np.nanmedian(y[good]))
                    for bi in np.where(out_mask)[0]:
                        lo, hi = max(0, bi - 20), min(n, bi + 21)
                        local_good = good[lo:hi]
                        y[bi] = float(np.nanmedian(y[lo:hi][local_good])) if local_good.any() else global_med
                tsys_t[ch_i, i_start:i_end] = y

    # --- Optional read-only plot of the final smoothed result ---
    if args.plot:
        bbc_labels = _channel_labels(header[0]) if header else []
        _show_result_plot(tsys_t, time_array, bbc_labels, [(s, e) for s, e, _ in part_boundaries])

    # --- Write the ANTAB file (only now, so aborting the editor never leaves a partial file) ---
    dpfu_lines, polyelev_lines = antab_header.write_antab_preamble(antab_file)
    # SEFD mode: DPFU from the final (smoothed) Tsys so that median(Tsys)/DPFU = SEFD
    if sefd:
        dpfu_lines, polyelev_lines = sefd_gain_lines(log_f.station(), header, setup_time, tsys_t, time_array, sefd)
        if not dpfu_lines:
            _fatal("--sefd: no valid Tsys values to derive DPFU from.")

    write_antab(antab_file, header, index_line, scan_line, tsys_t.T, block_data, time_data,
        tsys_log, setup_time, dpfu_lines, polyelev_lines, log_f.station())
    print(f"ANTAB file written: {antab_file}")

    if args.vlbeer:
        _upload_to_vlbeer(antab_file, log_f)


def _show_result_plot(
    tsys_t: np.ndarray,
    time_array: np.ndarray,
    bbc_labels: list[str],
    part_boundaries: list[tuple[int, int]],
) -> None:
    """Display a read-only summary plot of the processed Tsys data.

    Uses the same backend negotiation as the interactive editor via the
    plotting module import side-effect.
    """
    import datetime

    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    utc = datetime.timezone.utc
    date_fmt = mdates.DateFormatter("%j %H:%M")

    for part_num, (i_start, i_end) in enumerate(part_boundaries):
        if i_start >= i_end:
            continue
        x_dt = [datetime.datetime.fromtimestamp(float(t), tz=utc) for t in time_array[i_start:i_end]]
        fig, ax = plt.subplots(figsize=(16, 9))
        fig.suptitle(f"Processed Tsys — Part {part_num + 1}  (close window to continue)")
        ax.set_xlabel("Time (UT)")
        ax.set_ylabel("Tsys [K]")
        ax.grid(True, alpha=0.35)
        ax.xaxis.set_major_formatter(date_fmt)

        n_ch = min(tsys_t.shape[0], len(bbc_labels))
        for ch_i in range(n_ch):
            ax.plot(x_dt, tsys_t[ch_i, i_start:i_end], ".", ms=3, label=bbc_labels[ch_i])

        ax.legend(loc="upper right", fontsize="small", markerscale=2, framealpha=0.85)
        fig.autofmt_xdate()
        plt.tight_layout()
        plt.show()


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


def _channel_labels(header_block: str) -> list[str]:
    """Return channel labels (e.g. ``['R1', 'R2', 'L1']``) from the ``!Column N = R1: ...`` lines of one header."""
    return [line.split("=")[1].split(":")[0].strip() for line in header_block.split("\n")
            if "Column" in line and len(line.split("=")) > 1]


def _calib_info(log_f: LogFile, sefd: dict[str, float] | None, setup: str) -> list[str]:
    """Calibration source for *setup* as display lines (one item per line) for the interactive editor panel."""
    if sefd:
        return ["No RXG file (nominal).", f"SEFD(R) = {sefd['R']:g} Jy", f"SEFD(L) = {sefd['L']:g} Jy"]
    paths = log_f.rxg_files_needed().get(setup, [])
    names = sorted({os.path.basename(p) for p in paths if p.strip()})
    return [f"RXG file(s), setup {setup}:"] + names if names else [f"No RXG file for setup {setup}."]


def _check_rxg(log_f: LogFile, rxg_files: list[str] | None, user_rxg: bool) -> None:
    """Exit with an error unless every LO has an RXG file and every channel got an RXG Tcal.

    Invariant: Tcal, DPFU and gain curve all come from RXG files; /caltemp/ lines are never mixed in or used as
    a fallback (DPFU and gain curve exist only in RXG files). Without -f/-r (*user_rxg* False) a missing RXG is
    first reported as a warning, together with whether the log has /caltemp/ lines, before the final error.
    """
    missing = log_f.missing_rxg_los()
    searched = f"  Searched directories: {', '.join(log_f.rxg_dirs())}" + \
        (f" (files: {', '.join(rxg_files)})" if rxg_files else "")
    if missing:
        los = "\n".join(f"    setup {s}: LO={lo:.2f} MHz {pol}" for s, lo, pol in missing)
        if user_rxg:
            _fatal(f"the RXG file(s) given with -r/-f do not cover the following LO(s):\n{los}\n{searched}")
        print(f"Warning: no RXG file found for the following LO(s):\n{los}\n{searched}", file=sys.stderr)
        if log_f.has_caltemp():
            print("Warning: the log has /caltemp/ lines, but they only give Tcal: DPFU and the gain curve exist only "
                  "in RXG files.", file=sys.stderr)
            _fatal("RXG file not found. Provide it with -r/--rxg-dir or -f/--rxg-files, "
                   "or use --sefd for a nominal ANTAB file.")
        _fatal("RXG file not found and no /caltemp/ lines in the log. Provide the RXG file with -r/--rxg-dir or "
               "-f/--rxg-files, or use --sefd for a nominal ANTAB file.")
    no_tcal = log_f.missing_tcal_channels()
    if no_tcal:
        chans = ", ".join(f"setup {s} {ch}" for s, ch in no_tcal)
        _fatal(f"the RXG file(s) have no Tcal table entries for the polarization of: {chans}.\n{searched}")


def _fatal(message: str) -> None:
    """Print ``Error: <message>`` to stderr (red if stderr is a terminal) and exit with status 1."""
    text = f"Error: {message}"
    if sys.stderr.isatty():
        text = f"\033[1;31m{text}\033[0m"
    sys.stdout.flush()
    print(text, file=sys.stderr)
    sys.exit(1)


def _find_index(time_list: list[float], target: float) -> int | None:
    """Return index of first element >= target, or None (bisect O(log n))."""
    i = bisect.bisect_left(time_list, target)
    return i if i < len(time_list) else None


if __name__ == "__main__":
    main()

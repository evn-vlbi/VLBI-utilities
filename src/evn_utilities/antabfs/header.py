"""ANTAB file header generation.

Builds the preamble, DPFU/POLY lines, and LO information blocks
that precede the Tsys data in an ANTAB file.
"""

from __future__ import annotations

import datetime

import numpy as np

from . import __version__ as version
from .logfile import LogFile
from .rxg import RxgFile


class AntabHeader:
    """Create the header sections of an ANTAB file.

    If *sefd* (``{'R': Jy, 'L': Jy}``) is given, the header runs in SEFD mode: no RXG file is read, and
    the GAIN/POLY lines must be built afterwards with :func:`sefd_gain_lines`.
    """

    def __init__(self, log_file: LogFile, sefd: dict[str, float] | None = None) -> None:
        self.log_f = log_file
        self.exp_name = log_file.experiment()
        self.station_name = log_file.station()
        self.sefd = sefd

    # ----------------------------------------------------------------
    # RXG information lines
    # ----------------------------------------------------------------

    def rxg_lines(self) -> list[str]:
        """Lines describing the RXG files used for each LO."""
        lines: list[str] = []
        flo_array, pol_array = self.log_f.lo_p_array()

        for setup in reversed(list(flo_array.keys())):
            for i, flo in enumerate(flo_array[setup]):
                rxg_name = self.log_f.get_rxg_filename(flo)
                print(rxg_name)
                rxg = RxgFile(rxg_name)
                lines.append(
                    f"{flo:.2f} MHz {pol_array[setup][i]}: {rxg.name} {rxg.date()}"
                )
        return lines

    # ----------------------------------------------------------------
    # DPFU lines
    # ----------------------------------------------------------------

    def dpfu_lines(self) -> dict[str, list[str]]:
        """GAIN / DPFU lines per setup."""
        result: dict[str, list[str]] = {}
        rxg_files_array = self.log_f.rxg_files_needed()
        log_data = self.log_f.get_log_data()
        header = log_data[0]

        for setup in reversed(list(rxg_files_array.keys())):
            result[setup] = []
            for rxg_name in sorted(set(rxg_files_array[setup])):
                if rxg_name.strip() == "":
                    continue
                rxg = RxgFile(rxg_name)
                dpfu_str = ",".join(rxg.dpfu()).rstrip(",")
                freq_min, freq_max = rxg.freq_min_max()
                line = f"GAIN {self.station_name} ELEV DPFU={dpfu_str}"

                # Find matching header block for this setup
                h_lines: list[str] = []
                for h in header:
                    h_parts = h.split("\n")
                    if len(h_parts) > 1 and setup in h_parts[1]:
                        h_lines = h_parts[4:]
                        break

                if_freq: list[float] = []
                if_bw: list[float] = []
                for hl in h_lines:
                    tokens = hl.split()
                    if len(tokens) < 12:
                        continue
                    freq_val = float(tokens[6])
                    bw_val = float(tokens[11])
                    if freq_val >= (freq_min - bw_val) and freq_val <= (freq_max + bw_val):
                        if_freq.append(freq_val)
                        if_bw.append(bw_val)

                if not if_freq:
                    if_freq = [freq_min, freq_max]
                    if_bw = [0.0, 0.0]

                min_idx = if_freq.index(min(if_freq))
                max_idx = if_freq.index(max(if_freq))
                line += f" FREQ={min(if_freq) - if_bw[min_idx]:.2f},{max(if_freq) + if_bw[max_idx]:.2f}"
                result[setup].append(line)

        return result

    # ----------------------------------------------------------------
    # POLY/ELEV lines
    # ----------------------------------------------------------------

    def polyelev_line(self) -> dict[str, list[str]]:
        """Gain-curve polynomial lines per setup."""
        result: dict[str, list[str]] = {}
        rxg_files_array = self.log_f.rxg_files_needed()

        for setup in reversed(list(rxg_files_array.keys())):
            result[setup] = []
            for rxg_name in sorted(set(rxg_files_array[setup])):
                if rxg_name.strip() == "":
                    continue
                rxg = RxgFile(rxg_name)
                gain_list = rxg.gain()
                coeffs = ",".join(gain_list[2:])
                if coeffs.endswith(","):
                    coeffs = coeffs[:-1]
                result[setup].append(f"POLY={coeffs} /\n")

        return result

    # ----------------------------------------------------------------
    # LO information lines
    # ----------------------------------------------------------------

    def lo_lines(self) -> list[str]:
        """Comment lines listing LO frequencies, polarizations and RXG sources."""
        lines: list[str] = []
        flo_array, pol_array = self.log_f.lo_p_array()

        setups = sorted(flo_array.keys())
        for setup in setups:
            lines.append(f"!   Setup {setup}")
            lo_lines_list: list[str] = []
            freq_list: list[float] = []
            for i, flo in enumerate(flo_array[setup]):
                rxg_name = " " if self.sefd else self.log_f.get_rxg_filename(flo)
                if rxg_name.strip() == "":
                    lo_lines_list.append(f"!     LO={flo:.2f} MHz {pol_array[setup][i]}")
                else:
                    rxg = RxgFile(rxg_name)
                    lo_lines_list.append(
                        f"!     LO={flo:.2f} MHz {pol_array[setup][i]} "
                        f"{rxg.name} {rxg.date()}"
                    )
                freq_list.append(flo)

            lines.extend(x for _, x in sorted(zip(freq_list, lo_lines_list)))

        return lines

    # ----------------------------------------------------------------
    # Preamble
    # ----------------------------------------------------------------

    def first_lines(self) -> list[str]:
        """Informational comment lines at the top of the ANTAB file."""
        today = datetime.datetime.now().strftime("%Y-%m-%d")
        wavebands = self.wave_bands()
        wb_str = " ".join(f"{w}cm" for w in wavebands) + "."
        dbbc_mode = self.log_f.dbbc_mode()

        lines = [f"! Amplitude calibration data for {self.station_name} in {self.exp_name}.",
            "! For use with AIPS task ANTAB.", f"! Waveband(s) =  {wb_str}"]
        if self.sefd:
            lines += ["! NOMINAL calibration: no RXG files used.",
                f"! DPFU = median(Tsys) / SEFD per setup and polarization, with SEFD(R) = {self.sefd['R']:g} Jy, "
                f"SEFD(L) = {self.sefd['L']:g} Jy.",
                "! Tcal = 1 for all channels (Tsys in Tcal units, scaled to the SEFD by DPFU). "
                "No gain curve available (POLY=1.0).", "! LOs used:"]
        else:
            lines.append("! RXG files used for each LO:")
        lines.extend(self.lo_lines())
        lines.append(f"! DBBC used in mode {dbbc_mode}")

        version_date = datetime.datetime.strptime(str(version), "%Y%m%d").strftime("%Y-%m-%d")
        lines.append(f"! Produced on {today} using antabfs.py version: {version_date}")

        return lines

    def write_antab_preamble(self, antab_file: str) -> list[dict[str, list[str]]]:
        """Write the preamble to *antab_file* and return DPFU/POLY data.

        In SEFD mode both returned dicts are empty; use :func:`sefd_gain_lines` once Tsys is final.
        """
        with open(antab_file, "w") as fh:
            for line in self.first_lines():
                fh.write(line + "\n")

        if self.sefd:
            return [{}, {}]
        return [self.dpfu_lines(), self.polyelev_line()]

    # ----------------------------------------------------------------
    # Waveband helpers
    # ----------------------------------------------------------------

    def wave_bands(self) -> list[str]:
        """Compute waveband names (cm) from LO frequencies."""
        flo_array = self.log_f.lo_array()
        bands: list[str] = []
        for setup in reversed(list(flo_array.keys())):
            for freq_lo in flo_array[setup]:
                bands.append(f"{3e8 / (freq_lo * 1e4):.1f}")
        return sorted(set(bands), key=float, reverse=True)


def parse_sefd(text: str) -> dict[str, float]:
    """Parse ``'350'`` or ``'350,400'`` (RCP,LCP) into ``{'R': Jy, 'L': Jy}``.

    Raises ValueError if the format is wrong or any value is not a positive number.
    """
    parts = [p.strip() for p in text.split(",")]
    if len(parts) not in (1, 2):
        raise ValueError(f"expected one value or two comma-separated values (RCP,LCP), got {text!r}")
    values = [float(p) for p in parts]
    if any(not np.isfinite(v) or v <= 0 for v in values):
        raise ValueError(f"SEFD values must be positive numbers, got {text!r}")
    return {"R": values[0], "L": values[-1]}


def _setup_columns(header_block: str) -> tuple[list[str], list[float], list[float]]:
    """Return (pol letters, sky freqs MHz, bandwidths MHz) from the ``!Column`` lines of one header block.

    Assumes the column line format written by ``LogFile._fill_header``:
    ``!Column 1 = R1: ifA, bbc01, 1610.49 MHz , USB, BW= 32.00 MHz, Tcal=1.77 K``.
    """
    pols: list[str] = []
    freqs: list[float] = []
    bws: list[float] = []
    for line in header_block.split("\n"):
        tokens = line.split()
        if line.startswith("!Column") and len(tokens) >= 12:
            pols.append(tokens[3][0])
            freqs.append(float(tokens[6]))
            bws.append(float(tokens[11]))
    return pols, freqs, bws


def sefd_gain_lines(station: str, header: list[str], setup_time: list[list], tsys_t: np.ndarray, times: np.ndarray,
                    sefd: dict[str, float]) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """Build nominal GAIN/POLY lines so that median(Tsys)/DPFU equals the given SEFD.

    Parameters
    ----------
    station : two-letter station code.
    header : per-setup header blocks (same order as *setup_time*), containing the ``!Column`` lines.
    setup_time : ``[[unix_start, setup_name], ...]`` setup boundaries.
    tsys_t : final Tsys matrix, shape ``(n_channels, n_times)``, in K.
    times : Unix timestamps of the Tsys samples, shape ``(n_times,)``.
    sefd : ``{'R': Jy, 'L': Jy}``.

    Returns
    -------
    (dpfu_lines, poly_lines) keyed by setup name, in the format expected by ``writer.write_antab``.
    One GAIN line per setup; DPFU is listed as R,L (or a single value if only one pol is present).
    Setups that appear several times in *setup_time* pool all their samples.
    """
    samples: dict[str, dict[str, list[np.ndarray]]] = {}
    columns: dict[str, tuple[list[str], list[float], list[float]]] = {}
    for i, (t_start, setup) in enumerate(setup_time):
        t_end = setup_time[i + 1][0] if i + 1 < len(setup_time) else np.inf
        mask = (times >= t_start) & (times < t_end)
        pols, freqs, bws = _setup_columns(header[i]) if i < len(header) else ([], [], [])
        if not pols:
            continue
        columns.setdefault(setup, (pols, freqs, bws))
        for ch_i, pol in enumerate(pols[:tsys_t.shape[0]]):
            samples.setdefault(setup, {}).setdefault(pol, []).append(tsys_t[ch_i, mask])

    dpfu_lines: dict[str, list[str]] = {}
    poly_lines: dict[str, list[str]] = {}
    for setup, per_pol in samples.items():
        dpfu: list[str] = []
        for pol in ("R", "L"):
            vals = np.concatenate(per_pol.get(pol, [np.array([])]))
            vals = vals[np.isfinite(vals) & (vals > 0)]
            if len(vals) == 0:
                continue
            median_tsys = float(np.median(vals))
            dpfu.append(f"{median_tsys / sefd[pol]:.6f}")
            print(f"  SEFD mode: setup {setup} pol {pol}: median Tsys = {median_tsys:.1f} K, "
                  f"SEFD = {sefd[pol]:g} Jy -> DPFU = {dpfu[-1]} K/Jy")
        if not dpfu:
            continue
        _, freqs, bws = columns[setup]
        f_min = min(f - abs(b) for f, b in zip(freqs, bws))
        f_max = max(f + abs(b) for f, b in zip(freqs, bws))
        dpfu_lines[setup] = [f"GAIN {station} ELEV DPFU={','.join(dpfu)} FREQ={f_min:.2f},{f_max:.2f}"]
        poly_lines[setup] = ["POLY=1.0 /\n"]
    return dpfu_lines, poly_lines

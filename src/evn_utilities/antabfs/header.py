"""ANTAB file header generation.

Builds the preamble, DPFU/POLY lines, and LO information blocks
that precede the Tsys data in an ANTAB file.
"""

from __future__ import annotations

import datetime

from . import __version__ as version
from .logfile import LogFile
from .rxg import RxgFile


class AntabHeader:
    """Create the header sections of an ANTAB file."""

    def __init__(self, log_file: LogFile) -> None:
        self.log_f = log_file
        self.exp_name = log_file.experiment()
        self.station_name = log_file.station()

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
                rxg_name = self.log_f.get_rxg_filename(flo)
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
            "! For use with AIPS task ANTAB.", f"! Waveband(s) =  {wb_str}", "! RXG files used for each LO:"]
        lines.extend(self.lo_lines())
        lines.append(f"! DBBC used in mode {dbbc_mode}")

        version_date = datetime.datetime.strptime(str(version), "%Y%m%d").strftime("%Y-%m-%d")
        lines.append(f"! Produced on {today} using antabfs.py version: {version_date}")

        return lines

    def write_antab_preamble(self, antab_file: str) -> list[dict[str, list[str]]]:
        """Write the preamble to *antab_file* and return DPFU/POLY data."""
        with open(antab_file, "w") as fh:
            for line in self.first_lines():
                fh.write(line + "\n")

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

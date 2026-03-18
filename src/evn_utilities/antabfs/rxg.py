"""RXG file parser for EVN antenna calibration data.

RXG files contain receiver calibration information used by the Field System.
Format:
    Line 1: LO values and ranges (e.g. 'range 4000 4300' or 'fixed 4158')
    Line 2: Creation date (yyyy mm dd)
    Line 3: FWHM beamwidth (e.g. 'frequency value' or 'constant value')
    Line 4: Polarizations available (e.g. 'lcp rcp')
    Line 5: DPFU (degrees/Jansky) per polarization
    Line 6: Gain curve polynomial (e.g. 'ELEV POLY 0.95 0.002 -3.2e-05')
    Line 7+: Tcal vs frequency (e.g. 'lcp 4650.0 1.5')
    Then: Trec (receiver temperature)
    Then: Spillover vs elevation
"""

from __future__ import annotations

from pathlib import Path


class RxgFile:
    """Parse and query an RXG (receiver gain) calibration file."""

    def __init__(self, filename: str | Path) -> None:
        self._path = Path(filename)
        self.rxgname = self._path.name
        with open(self._path) as fh:
            self._lines = fh.readlines()

    # ------------------------------------------------------------------
    # Low-level line access
    # ------------------------------------------------------------------

    def _get_line_from_param(self, parameter: str) -> list[str]:
        """Return content lines for a named parameter section.

        Parameters are identified by their canonical order (skipping comment
        lines that start with ``*``).  TCAL and SPILL sections span multiple
        non-comment lines.
        """
        param = parameter.upper()
        param_numbers = {
            "LO": 1, "DATE": 2, "FWHM": 3, "POLS": 4,
            "DPFU": 5, "GAIN": 6, "TCAL": 7, "TREC": 8, "SPILL": 9,
        }
        if param not in param_numbers:
            raise ValueError(f"Unknown RXG parameter: {param!r}")

        target = param_numbers[param]
        result: list[str] = []
        counter = 0
        in_multi = False

        for line in self._lines:
            if line.startswith("*"):
                if in_multi:
                    break
                continue
            counter += 1
            if counter == target:
                if target in (7, 9):  # multi-line sections
                    in_multi = True
                    result.append(line)
                    # stay at same counter so next non-comment line also matches
                    counter -= 1
                else:
                    result.append(line)
                    break
            elif in_multi:
                result.append(line)

        return result

    # ------------------------------------------------------------------
    # Public accessors
    # ------------------------------------------------------------------

    @property
    def name(self) -> str:
        """Filename without path."""
        return self.rxgname

    def date(self) -> str:
        """Creation / modification date string from the RXG file."""
        return self._get_line_from_param("Date")[0].strip("\n")

    def pols(self) -> list[str]:
        """Available polarizations, e.g. ``['lcp', 'rcp']``."""
        return self._get_line_from_param("POLS")[0].split()

    def dpfu(self) -> list[str]:
        """DPFU values as strings, ordered by polarization."""
        return self._get_line_from_param("DPFU")[0].split()

    def gain(self) -> list[str]:
        """Gain curve tokens (type, basis, then polynomial coefficients)."""
        return self._get_line_from_param("GAIN")[0].split()

    def lo(self) -> list[str]:
        """LO frequency value(s) as strings (1 element if fixed, 2 if range)."""
        tokens = self._get_line_from_param("LO")[0].split()
        # tokens[0] is 'range' or 'fixed'; rest are frequency values
        return tokens[1:]

    def trec(self) -> list[str]:
        """Receiver temperature value(s) as strings."""
        return self._get_line_from_param("TREC")[0].split()

    def tcal(self) -> list[str]:
        """Raw Tcal lines (excluding ``end_tcal_table`` sentinel)."""
        result = []
        for line in self._get_line_from_param("TCAL"):
            if "end_tcal_table" in line:
                continue
            result.append(line.strip("\n"))
        return result

    def freq_cal(self) -> tuple[list[float], list[float], list[float], list[float]]:
        """Parse Tcal table into per-polarization frequency and Tcal arrays.

        Returns
        -------
        freq_rcp, tcal_rcp, freq_lcp, tcal_lcp
        """
        freq_rcp: list[float] = []
        tcal_rcp: list[float] = []
        freq_lcp: list[float] = []
        tcal_lcp: list[float] = []

        for entry in self.tcal():
            parts = entry.split()
            if len(parts) == 0:
                continue
            pol = parts[0].lower()
            if pol == "rcp":
                freq_rcp.append(float(parts[1]))
                tcal_rcp.append(float(parts[2]))
            elif pol == "lcp":
                freq_lcp.append(float(parts[1]))
                tcal_lcp.append(float(parts[2]))

        return freq_rcp, tcal_rcp, freq_lcp, tcal_lcp

    def freq_min_max(self) -> tuple[float, float]:
        """Min and max frequencies covered by the Tcal table."""
        fr, _tr, fl, _tl = self.freq_cal()
        all_min: list[float] = []
        all_max: list[float] = []
        if fr:
            all_min.append(min(fr))
            all_max.append(max(fr))
        if fl:
            all_min.append(min(fl))
            all_max.append(max(fl))
        return min(all_min), max(all_max)

    def cal_vs_freq(self, pol: str) -> list[list[float]]:
        """Return ``[[freq, tcal], ...]`` pairs for a given polarization."""
        fr, tr, fl, tl = self.freq_cal()
        if pol.lower() == "rcp":
            return [[f, t] for f, t in zip(fr, tr)]
        elif pol.lower() == "lcp":
            return [[f, t] for f, t in zip(fl, tl)]
        return []

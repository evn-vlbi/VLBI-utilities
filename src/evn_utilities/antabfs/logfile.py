"""Field System log file parser.

Reads an FS log file and extracts DBBC configuration, scan metadata,
and system-temperature measurements needed to produce ANTAB files.
"""

from __future__ import annotations

import datetime
import itertools
import os
from typing import Any

from .rxg import RxgFile
from .utils import get_rxg_dirs, get_tcal


_EPOCH = datetime.datetime(1970, 1, 1, tzinfo=datetime.timezone.utc)

# Channel layout for each ``/form=`` type
_FORM_CHANNELS: dict[str, list[str]] = {
    "geo": [
        "1u", "2u", "3u", "4u", "5u", "6u", "7u", "8u",
        "1l", "8l", "9u", "au", "bu", "cu", "du", "eu",
    ],
    "astro": [
        "1u", "2u", "3u", "4u", "5u", "6u", "7u", "8u",
        "1l", "2l", "3l", "4l", "5l", "6l", "7l", "8l",
    ],
    "astro2": [
        "1u", "2u", "3u", "4u", "9u", "au", "bu", "cu",
        "1l", "2l", "3l", "4l", "9l", "al", "bl", "cl",
    ],
    "astro3": [
        "1u", "3u", "5u", "7u", "9u", "bu", "du", "fu",
        "1l", "3l", "5l", "7l", "9l", "bl", "dl", "fl",
    ],
    "lba": [
        "1u", "2u", "5u", "6u", "3u", "4u", "7u", "8u",
        "1l", "2l", "5l", "6l", "3l", "4l", "7l", "8l",
    ],
    "wastro": [
        "1u", "2u", "3u", "4u", "5u", "6u", "7u", "8u",
        "1l", "2l", "3l", "4l", "5l", "6l", "7l", "8l",
        "9u", "au", "bu", "cu", "du", "eu", "fu", "gu",
        "9l", "al", "bl", "cl", "dl", "el", "fl", "gl",
    ],
}

_PFB_FREQ = [
    1040, 1008, 976, 944, 912, 880, 848, 816,
    784, 752, 720, 688, 656, 624, 592, 560,
]


class LogFile:
    """Parse an FS log file and extract DBBC / LO / Tsys information."""

    def __init__(
        self,
        filename: str,
        rxg_dir: str | None = None,
        rxg_files: list[str] | None = None,
        debug: bool = False,
    ) -> None:
        self._rxg_dirs: list[str] = [rxg_dir] if rxg_dir is not None else get_rxg_dirs()
        self._rxg_files = rxg_files
        self._debug = debug

        self.logname = os.path.basename(filename)
        exp_station = self.logname.split(".")[0]
        self.station_name = exp_station[-2:].upper()
        self.exp_name = exp_station[:-2].lower()

        with open(filename) as fh:
            self._file_content = fh.readlines()

        # Public state populated during parsing
        self.freq_lo_mhz: dict[str, dict[str, list[float]]] = {}
        self.ifd_setup: dict[str, str] = {}
        self.pol_array: dict[str, dict[str, list[str]]] = {}
        self.band_array: dict[str, dict[str, list[str]]] = {}
        self.dbbc_mode_name: str | None = None
        self.cal_mode_name: dict[str, str] = {}
        self.last_cal_mode: str | None = None
        self.log_data: list[Any] = []

        # Private parsing state
        self._bbc_info: dict[str, list[list[str]]] = {}
        self._vsi_ch: dict[str, list[list[str]]] = {}
        self._form_type: str | None = None
        self._ch_id: list[str] | None = None
        self._ch_id_index: int | None = None
        self._scan_num = 0
        self._scan_name: str | None = None
        self._data_valid = False
        self._dt_start_int: datetime.datetime | None = None
        self._int_complete = False
        self._temp_dict: list[dict[str, list[float]]] = [{}, {}, {}, {}]
        self._header_comp = False
        self._bbc_code_list: dict[str, list[str]] = {}
        self._bw: dict[str, dict[str, float]] = {}
        self._bbc_fq: dict[str, dict[str, float]] = {}
        self._which_if: dict[str, dict[str, str]] = {}
        self._int_time = datetime.timedelta(seconds=1)
        self._tsys_log_dict: dict[float, dict[str, list[float]]] = {}

        self._scan_line: list[str] = []
        self._index_line: list[str] = []
        self._header: list[str] = []
        self._tsys_log: list[list[float]] = []
        self._setup_time: list[list[Any]] = []
        self._setup_tcal: dict[str, dict[str, list[float]]] = {}
        self._caltemp_read: dict[str, bool] = {}
        self._current_setup: str | None = None
        self._last_setup: str | None = None
        self._new_setup = False
        self._fila10g_mode: dict[str, int] = {}
        self._rec_mode: dict[str, int] = {}
        self._ch_id_set: set[str] = set()

        self._read_log()

    # ----------------------------------------------------------------
    # Main parsing loop
    # ----------------------------------------------------------------

    def _read_log(self) -> None:
        """Read the entire LOG file, populating header/scan/Tsys data."""
        time_list: list[float] = []
        block_list: list[int] = []
        tsys_line_list: list[list[float]] = []

        lines = self._file_content
        n = len(lines)
        for n_line, line in enumerate(lines):
            if not line.strip():
                continue

            # --- Header variables ---
            if not self._header_comp:
                if self._read_header(line):
                    continue

            # --- General variables ---
            if self._get_gen_var(line):
                continue

            # --- Data validation ---
            if self._check_data_valid(line):
                continue

            # --- Temperature variables ---
            if self._data_valid:
                next_line = lines[n_line + 1].strip() if n_line + 1 < n else ""
                self._read_temp_var(line, next_line)

            # --- Integration complete → Tsys calculation ---
            if self._int_complete:
                if self._current_setup not in self._bbc_code_list:
                    self._int_complete = False
                    continue

                dt = self._get_datetime(line)
                if not dt:
                    self._int_complete = False
                    continue

                tsys_aux = self._get_tsys(
                    self._bbc_code_list[self._current_setup], dt
                )

                if tsys_aux and not all(t == -1 for t in tsys_aux):
                    # dt is already parsed from this line — no need to re-parse.
                    total_seconds = (
                        dt.replace(tzinfo=datetime.timezone.utc) - _EPOCH
                    ).total_seconds()
                    time_list.append(total_seconds)
                    block_list.append(self._scan_num)
                    tsys_line_list.append(tsys_aux)

                self._int_complete = False

        self._fill_header()

        self.log_data = [self._header, self._index_line, self._scan_line, tsys_line_list,
                         block_list, time_list, self._tsys_log, self._setup_time]

    # ----------------------------------------------------------------
    # Header filling
    # ----------------------------------------------------------------

    def _fill_header(self) -> None:
        """Build ANTAB-format header strings and Tsys-from-log data."""
        for bp in range(len(self._setup_time)):
            setup = self._setup_time[bp][1]
            col_num = 1
            pol_num = {"L": 1, "R": 1}
            header_str = f"!\n! Setup {setup}\n! Calibration mode: {self.cal_mode_name.get(setup, 'UNKNOWN')}\n!\n"
            index_str = "INDEX= "

            if setup not in self._bbc_code_list:
                self._header.append(header_str[:-1])
                self._index_line.append(index_str + "\n")
                continue

            for ch in self._bbc_code_list[setup]:
                lo_fq = self.freq_lo_mhz[setup][self._which_if[setup][ch]][-1]
                pol = self.pol_array[setup][self._which_if[setup][ch]][-1]
                band = self.band_array[setup][self._which_if[setup][ch]][-1]
                bw_aux = self._bw[setup][ch]

                if self.dbbc_mode_name == "DDC":
                    if ch[-1] == "l":
                        bw_aux = -bw_aux
                    if "lsb" in band:
                        bbc_freq = -self._bbc_fq[setup][ch]
                    else:
                        bbc_freq = self._bbc_fq[setup][ch]
                    fchan = lo_fq + bbc_freq + bw_aux / 2.0
                    aux_str = f"0{ch[0]}"
                    if aux_str == "0g":
                        bbc_num = 16
                    else:
                        bbc_num = int(aux_str, 16)
                    sb_letter = ch[-1].upper()

                elif self.dbbc_mode_name == "PFB":
                    if "lsb" in band:
                        bbc_freq = -self._bbc_fq[setup][ch]
                    else:
                        bbc_freq = self._bbc_fq[setup][ch]
                    fchan = lo_fq + bbc_freq - bw_aux / 2.0
                    bbc_num = int(ch[1:])
                    sb_letter = "L"
                else:
                    continue

                tcal_val = self._setup_tcal[setup][ch][0]
                pol_char = pol[0].upper()

                header_str += (f"!Column {col_num} = {pol_char}{pol_num[pol_char]}: "
                    f"if{self._which_if[setup][ch].upper()}, bbc{bbc_num:02d}, "
                    f"{fchan:.2f} MHz , {sb_letter}SB, BW= {self._bw[setup][ch]:04.2f} MHz, "
                    f"Tcal={tcal_val:.2f} K\n")
                index_str += f"'{pol_char}{pol_num[pol_char]}',"
                col_num += 1
                pol_num[pol_char] += 1

            index_str = index_str.rstrip(",")
            self._header.append(header_str[:-1])
            self._index_line.append(index_str + "\n")

        # --- Tsys from the LOG file ---
        sorted_keys = sorted(self._tsys_log_dict)
        for ts in sorted_keys:
            tsys_entry: list[float] = [ts]
            aux_dict = self._tsys_log_dict[ts]
            if self._current_setup and self._current_setup in self._bbc_code_list:
                for ch in self._bbc_code_list[self._current_setup]:
                    try:
                        vals = aux_dict[ch]
                        tsys_entry.append(sum(vals) / len(vals))
                    except (KeyError, ZeroDivisionError):
                        tsys_entry.append(-1)
            self._tsys_log.append(tsys_entry)

    # ----------------------------------------------------------------
    # Setup parameter initialization
    # ----------------------------------------------------------------

    def _set_params(self) -> None:
        """Initialize BBC code lists, frequencies, bandwidths for current setup."""
        setup = self._current_setup
        if setup is None:
            return

        if setup not in self._bbc_code_list:
            self._bbc_code_list[setup] = []
            self._bbc_fq[setup] = {}
            self._which_if[setup] = {}
            self._bw[setup] = {}

            if setup not in self.cal_mode_name:
                print(
                    f"Calibration mode not found for setup {setup}.\n"
                    f"SINGLE calibration mode will be assumed as used by setup {setup}."
                )
                self.cal_mode_name[setup] = "SINGLE"
                self.last_cal_mode = "SINGLE"
                if setup not in self._setup_tcal:
                    self._setup_tcal[setup] = {}
                    self._caltemp_read[setup] = False
                if len(self._temp_dict) < 5:
                    self._temp_dict = [{}] + self._temp_dict

            if self.dbbc_mode_name == "PFB":
                if setup in self._vsi_ch:
                    self._bbc_code_list[setup] = self._vsi_ch[setup][0] + self._vsi_ch[setup][1]
                self._bbc_code_list[setup].sort()
                for ch in self._bbc_code_list[setup]:
                    self._bw[setup][ch] = 32.0
                    self._which_if[setup][ch] = ch[0]
                    self._bbc_fq[setup][ch] = _PFB_FREQ[int(ch[1:])]

            elif self.dbbc_mode_name == "DDC":
                if setup not in self._bbc_info:
                    self._bbc_info[setup] = []
                self._bbc_info[setup].sort()
                self._bbc_info[setup] = list(k for k, _ in itertools.groupby(self._bbc_info[setup]))
                bbc_info = self._bbc_info[setup]

                for i in range(len(bbc_info)):
                    aux = int(bbc_info[i][0][-2:])
                    if (
                        setup in self._fila10g_mode
                        and self._form_type in _FORM_CHANNELS
                    ):
                        channels = _FORM_CHANNELS[self._form_type]
                        bbc_num_index = [a for a in range(len(channels))
                            if (f"{aux:x}" in channels[a]) or (aux == 16 and "g" in channels[a])]
                        for idx in bbc_num_index:
                            aux_code = 2 ** (idx * 2) + 2 ** (idx * 2 + 1)
                            if aux_code & self._fila10g_mode[setup]:
                                code = channels[idx]
                                self._bbc_fq[setup][code] = float(bbc_info[i][1])
                                self._which_if[setup][code] = bbc_info[i][2]
                                self._bw[setup][code] = float(bbc_info[i][3])
                                self._bbc_code_list[setup].append(code)

                    elif (
                        setup in self._rec_mode
                        and self._form_type in _FORM_CHANNELS
                    ):
                        channels = _FORM_CHANNELS[self._form_type]
                        bbc_num_index = [a for a in range(len(channels)) if f"{aux:x}" in channels[a]]
                        for idx in bbc_num_index:
                            aux_code = 2 ** (idx * 2) + 2 ** (idx * 2 + 1)
                            if aux_code & self._rec_mode[setup]:
                                code = channels[idx]
                                self._bbc_fq[setup][code] = float(bbc_info[i][1])
                                self._which_if[setup][code] = bbc_info[i][2]
                                self._bw[setup][code] = float(bbc_info[i][3])
                                self._bbc_code_list[setup].append(code)
                    else:
                        if self._form_type != "geo2":
                            if self._form_type != "geo" or aux in [1, 8]:
                                lcode = f"{hex(aux)[-1]}l"
                                self._bbc_fq[setup][lcode] = float(bbc_info[i][1])
                                self._which_if[setup][lcode] = bbc_info[i][2]
                                self._bw[setup][lcode] = float(bbc_info[i][3])
                                self._bbc_code_list[setup].append(lcode)
                        ucode = f"{hex(aux)[-1]}u"
                        self._bbc_fq[setup][ucode] = float(bbc_info[i][1])
                        self._which_if[setup][ucode] = bbc_info[i][2]
                        self._bw[setup][ucode] = float(bbc_info[i][3])
                        self._bbc_code_list[setup].append(ucode)

                self._bbc_code_list[setup].sort()

            # Compute Tcal for each BBC channel from RXG files
            for ch in self._bbc_code_list[setup]:
                lo_fq = self.freq_lo_mhz[setup][self._which_if[setup][ch]][-1]
                pol = self.pol_array[setup][self._which_if[setup][ch]][-1]
                band = self.band_array[setup][self._which_if[setup][ch]][-1]
                bw_val = self._bw[setup][ch]

                if self.dbbc_mode_name == "DDC" and ch[-1] == "u":
                    bw_val = -bw_val
                if "lsb" in band:
                    bbc_freq = -self._bbc_fq[setup][ch]
                else:
                    bbc_freq = self._bbc_fq[setup][ch]

                fchan = lo_fq + bbc_freq - bw_val / 2.0
                tcal_val = get_tcal(lo_fq, pol, fchan, self.station_name.lower(),
                    rxg_dirs=self._rxg_dirs, rxg_files=self._rxg_files, debug=self._debug)

                if not self._caltemp_read.get(setup, False):
                    self._temp_dict[-1][ch] = [tcal_val]
                    self._setup_tcal[setup][ch] = [tcal_val]

        else:
            # Setup already configured — restore tcal from saved values
            for ch in self._bbc_code_list[setup]:
                self._temp_dict[-1][ch] = self._setup_tcal[setup][ch]

    # ----------------------------------------------------------------
    # Temperature variable reading
    # ----------------------------------------------------------------

    def _read_temp_var(self, line: str, next_line: str) -> None:
        """Read ``tpicd`` temperature variable from the LOG line."""
        if self._current_setup is None:
            return

        cal_mode = self.cal_mode_name.get(self._current_setup, "SINGLE")
        if cal_mode == "CONT":
            temp_ref = ["#tpicd#tpcont/"]
        else:
            temp_ref = ["#tpicd#tpi/"]

        temp_ind = self._id_line(line, temp_ref)

        if temp_ind != 0:
            if not self._header_comp and self._current_setup not in self.cal_mode_name:
                self._set_params()
                self._header_comp = True

            dt = self._get_datetime(line)
            if not dt:
                return

            if next_line == "":
                dt_next = dt + datetime.timedelta(milliseconds=200)
            else:
                dt_next = self._get_datetime(next_line)
                if not dt_next:
                    dt_next = dt + datetime.timedelta(milliseconds=200)

            if self._dt_start_int is not None:
                if (
                    (dt - self._dt_start_int) >= self._int_time
                    and (dt_next - dt) > datetime.timedelta(milliseconds=100)
                ):
                    self._int_complete = True
                    self._dt_start_int = dt

            self._get_temp_line(line, temp_ind)

    # ----------------------------------------------------------------
    # Data validation
    # ----------------------------------------------------------------

    def _check_data_valid(self, line: str) -> bool:
        """Check for ``data_valid=on/off`` lines."""
        id_idx = self._id_line(line, ["data_valid=on", "data_valid=off"])
        if id_idx == 0:
            return False

        if id_idx == 1:
            if self._current_setup is not None:
                self._data_valid = True
                self._dt_start_int = self._get_datetime(line)
                if self._new_setup:
                    self._set_params()
                self._new_setup = False
            return True

        if id_idx == 2:
            if self._current_setup is not None:
                if self._data_valid:
                    self._int_complete = True
                self._data_valid = False
            return True

        return False

    # ----------------------------------------------------------------
    # General variables
    # ----------------------------------------------------------------

    def _get_gen_var(self, line: str) -> bool:
        """Read general variables that can appear anywhere in the log."""

        # --- VSI / BBC channels ---
        if self.dbbc_mode_name == "PFB":
            vsi_num = self._id_line(line, ["/vsi1=", "/vsi2="])
            if vsi_num != 0:
                val = line.split("=")[1].split("\n")[0]
                if self._current_setup not in self._vsi_ch:
                    self._vsi_ch[self._current_setup] = [[], []]
                self._vsi_ch[self._current_setup][vsi_num - 1] = val.split(",")
                return True

        elif self.dbbc_mode_name == "DDC":
            if self._id_line(line, ["&dbbc"]) and not self._id_line(line, ["/if=bbc"]):
                parts = line.split("/")[1].split("=")
                bbc_str = parts[0]
                bbc_vals = parts[1].split(",")
                if self._current_setup not in self._bbc_info:
                    self._bbc_info[self._current_setup] = []
                self._bbc_info[self._current_setup].append(
                    [bbc_str, bbc_vals[0], bbc_vals[1], bbc_vals[2]]
                )
                return True

        # --- Scan name ---
        if self._id_line(line, [":scan_name="]):
            self._scan_name = line.split("=")[1].split(",")[0]
            self._scan_num += 1
            return True

        # --- Source ---
        if self._id_line(line, [":source="]):
            dt = self._get_datetime(line)
            if dt:
                days = (dt - datetime.datetime(dt.year, 1, 1, dt.hour, dt.minute, dt.second, dt.microsecond)).days + 1
                source_name = line.split("=")[1].split(",")[0]
                minute_frac = dt.minute + dt.second / 60.0 + dt.microsecond / 60e6
                self._scan_line.append(f"\n! {days:03d} {dt.hour:02d}:{minute_frac:05.2f}: "
                    f"scanNum={self._scan_num:04d} scanName={self._scan_name} source={source_name}")
            return True

        # --- Current setup ---
        if self._id_line(line, [":setup", ";setup", "/setup"]):
            self._current_setup = line.split("p")[-1].strip()
            if self._current_setup not in self._setup_tcal:
                self._setup_tcal[self._current_setup] = {}
                self._caltemp_read[self._current_setup] = False
            if self._current_setup != self._last_setup:
                dt = self._get_datetime(line)
                if dt:
                    time_val = (
                        dt.replace(tzinfo=datetime.timezone.utc) - _EPOCH
                    ).total_seconds()
                    self._setup_time.append([time_val, self._current_setup])
                self._last_setup = self._current_setup
                self._new_setup = True
                if self._current_setup in self.cal_mode_name:
                    if self.cal_mode_name[self._current_setup] == "SINGLE":
                        self._temp_dict = [{}, {}, {}, {}]
                    elif self.cal_mode_name[self._current_setup] == "CONT":
                        self._temp_dict = [{}, {}, {}]
                return True

        # --- ifd setup ---
        if self._id_line(line, ["&setup"]):
            parts = line.split("&")[1].split("/")
            cur_setup = parts[0].split("setup")[1]
            if "ifd" in parts[1]:
                if cur_setup not in self.ifd_setup:
                    self.ifd_setup[cur_setup] = parts[1].strip()

        # --- Fila10G mode ---
        if self._id_line(line, ["/fila10g_mode="]):
            tmp = line.split("=")[1]
            mask1 = tmp.split(",")[0][2:]
            mask2 = tmp.split(",")[1][2:]
            mask = f"0x{mask2}{mask1}"
            if self._current_setup is not None:
                self._fila10g_mode[self._current_setup] = int(mask, 16)
        elif self._id_line(line, ["_mode="]):
            if self._current_setup is not None:
                self._rec_mode[self._current_setup] = int(line.split(",")[1], 16)

        # --- Calibration mode ---
        if self._id_line(line, ["/cont_cal="]) and self._id_line(line, ["&setup"]):
            val = line.split("=")[1].split(",")[0].strip()
            if self._current_setup is not None:
                if val == "off":
                    self.cal_mode_name[self._current_setup] = "SINGLE"
                elif val == "on":
                    self.cal_mode_name[self._current_setup] = "CONT"

                if self.cal_mode_name[self._current_setup] != self.last_cal_mode:
                    self.last_cal_mode = self.cal_mode_name[self._current_setup]
                    if self.cal_mode_name[self._current_setup] == "SINGLE":
                        self._temp_dict = [{}, {}, {}, {}, {}]
                    elif self.cal_mode_name[self._current_setup] == "CONT":
                        self._temp_dict = [{}, {}, {}, {}]
            return True

        # --- LO configuration ---
        lo_refs = ["/lo=loa", "/lo=lob", "/lo=loc", "/lo=lod"]
        if self._id_line(line, lo_refs):
            if_id = line.split("&")[1].split("/")[0]
            if self._current_setup and self._current_setup in self.ifd_setup:
                if if_id != self.ifd_setup[self._current_setup]:
                    return True

            parts = line.split(",")
            if_sel = parts[0][-1]

            if self._current_setup in self.freq_lo_mhz:
                freq_aux = self.freq_lo_mhz[self._current_setup].get(if_sel, [])
                pol_aux = self.pol_array[self._current_setup].get(if_sel, [])
                band_aux = self.band_array[self._current_setup].get(if_sel, [])
            else:
                self.freq_lo_mhz[self._current_setup] = {}
                self.pol_array[self._current_setup] = {}
                self.band_array[self._current_setup] = {}
                freq_aux = []
                pol_aux = []
                band_aux = []

            try:
                new_freq = float(parts[1])
            except ValueError:
                print(f"Couldn't convert '{parts[1]}' to float. Ignoring line: '{line.strip()}'")
                return True

            new_pol = parts[3]
            new_band = parts[2]

            already_exists = any(
                f == new_freq and p == new_pol and b == new_band
                for f, p, b in zip(freq_aux, pol_aux, band_aux)
            )

            if not already_exists:
                freq_aux.append(new_freq)
                pol_aux.append(new_pol)
                band_aux.append(new_band)
                self.freq_lo_mhz[self._current_setup][if_sel] = freq_aux
                self.pol_array[self._current_setup][if_sel] = pol_aux
                self.band_array[self._current_setup][if_sel] = band_aux

            return True

        # --- Single-mode temperature variables ---
        if self._current_setup in self.cal_mode_name:
            if self.cal_mode_name[self._current_setup] == "SINGLE":
                aux_ref = ["/tpi/", "/tpical", "/tpdiff/", "/caltemp/"]
                temp_ind = self._id_line(line, aux_ref)
                if temp_ind != 0:
                    temp_ind += 1
                    self._get_temp_line(line, temp_ind)
                    return True

                if self._id_line(line, ["/tsys/"]):
                    self._store_tsys_from_log(line)
                    return True
            else:
                if self._id_line(line, ["/caltemp/"]):
                    self._get_temp_line(line, len(self._temp_dict))
                    return True
                if self._id_line(line, ["#tpicd#tsys/"]):
                    self._store_tsys_from_log(line)
                    return True
        else:
            aux_ref = ["/tpi/", "/tpical", "/caltemp/"]
            temp_ind = self._id_line(line, aux_ref)
            if temp_ind != 0:
                self._get_temp_line(line, temp_ind)
                return True

        return False

    def _store_tsys_from_log(self, line: str) -> None:
        """Parse inline Tsys from the log (``/tsys/`` or ``#tpicd#tsys/``)."""
        dt = self._get_datetime(line)
        if not dt:
            return
        time_val = (
            dt.replace(tzinfo=datetime.timezone.utc) - _EPOCH
        ).total_seconds()
        parts = line.split("/")[-1].split(",")

        if time_val not in self._tsys_log_dict:
            self._tsys_log_dict[time_val] = {}

        aux_dict = self._tsys_log_dict[time_val]
        if self._ch_id is None:
            return

        for i in range(0, len(parts), 2):
            if i + 1 >= len(parts):
                break
            key = parts[i]
            if not key or key[self._ch_id_index] not in self._ch_id_set:
                continue
            try:
                val = float(parts[i + 1])
            except (ValueError, IndexError):
                val = -1.0
            if key in aux_dict:
                aux_dict[key].append(val)
            else:
                aux_dict[key] = [val]

        self._tsys_log_dict[time_val] = aux_dict

    # ----------------------------------------------------------------
    # Tsys calculation
    # ----------------------------------------------------------------

    def _get_tsys(self, bbc_codes: list[str], dt: datetime.datetime) -> list[float]:
        """Calculate Tsys for each BBC channel."""
        if not self._header_comp and self._current_setup not in self.cal_mode_name:
            self.cal_mode_name[self._current_setup] = "SINGLE"
            self._set_params()
            self._header_comp = True

        cal_mode = self.cal_mode_name.get(self._current_setup, "SINGLE")
        temp: list[float] = []
        dt_new_order = datetime.datetime(2015, 9, 17)
        tp_zero = 0.0

        for ch in bbc_codes:
            if ch not in self._temp_dict[-1]:
                temp.append(-1)
                continue

            tcal_val = self._temp_dict[-1][ch][0]
            if tcal_val == 0:
                temp.append(-1)
                continue

            if cal_mode == "CONT":
                try:
                    vsys_on_list = self._temp_dict[1 if dt < dt_new_order else 0][ch]
                    vsys_off_list = self._temp_dict[0 if dt < dt_new_order else 1][ch]
                except (KeyError, IndexError):
                    temp.append(-1)
                    continue

                if not vsys_on_list or not vsys_off_list:
                    temp.append(-1)
                    continue

                vsys_on = sum(vsys_on_list) / len(vsys_on_list)
                vsys_off = sum(vsys_off_list) / len(vsys_off_list)

                if vsys_on <= vsys_off:
                    tsys = -1.0
                else:
                    tsys = 0.5 * tcal_val * (vsys_on + vsys_off) / (vsys_on - vsys_off)

            elif cal_mode == "SINGLE":
                try:
                    tpiprime_list = self._temp_dict[1].get(ch, [])
                    tpical_list = self._temp_dict[2].get(ch, [])
                    vsys_list = self._temp_dict[0].get(ch, [])
                    tpidiff_list = self._temp_dict[3].get(ch, []) if len(self._temp_dict) > 3 else []
                except Exception:
                    temp.append(-1)
                    continue

                if tpidiff_list:
                    if not vsys_list:
                        temp.append(-1)
                        continue
                    vsys = sum(vsys_list) / len(vsys_list)
                    tpidiff = sum(tpidiff_list) / len(tpidiff_list)
                    if tpidiff <= 0:
                        tsys = -1.0
                    else:
                        tsys = tcal_val * (vsys - tp_zero) / tpidiff
                else:
                    if not tpiprime_list or not tpical_list:
                        temp.append(-1)
                        continue
                    if not vsys_list:
                        temp.append(-1)
                        continue
                    tpiprime = sum(tpiprime_list) / len(tpiprime_list)
                    tpical = sum(tpical_list) / len(tpical_list)
                    vsys = sum(vsys_list) / len(vsys_list)
                    if tpical <= tpiprime:
                        tsys = -1.0
                    else:
                        tsys = tcal_val * (vsys - tp_zero) / (tpical - tpiprime)
            else:
                tsys = -1.0

            temp.append(tsys)

        # Clear temperature dictionaries
        self._temp_dict[0] = {}
        if cal_mode == "CONT":
            self._temp_dict[1] = {}

        return temp

    # ----------------------------------------------------------------
    # Temperature line parsing
    # ----------------------------------------------------------------

    def _get_temp_line(self, line: str, temp_ind: int) -> None:
        """Store temperature variables in the proper dictionary."""
        aux_dict = self._temp_dict[temp_ind - 1]
        parts = line.split("/")[-1].split(",")

        cal_mode = self.cal_mode_name.get(self._current_setup) if self._current_setup else None
        tpcont_det = (temp_ind == 1) and (cal_mode == "CONT")

        if tpcont_det:
            aux_range = range(0, len(parts), 3)
        else:
            aux_range = range(0, len(parts), 2)

        if self._ch_id is None:
            return

        for i in aux_range:
            if i >= len(parts):
                break
            ch_key = parts[i]
            if not ch_key or ch_key[self._ch_id_index] not in self._ch_id_set:
                continue
            dict_exist = (ch_key in aux_dict) and (temp_ind == 1) and (cal_mode is not None)
            if dict_exist:
                try:
                    aux_dict[ch_key].append(float(parts[i + 1]))
                    if tpcont_det:
                        self._temp_dict[1][ch_key].append(float(parts[i + 2]))
                except (ValueError, IndexError):
                    aux_dict[ch_key].append(-1)
                    if tpcont_det:
                        self._temp_dict[1].setdefault(ch_key, []).append(-1)
            else:
                try:
                    aux_dict[ch_key] = [float(parts[i + 1])]
                    if tpcont_det:
                        self._temp_dict[1][ch_key] = [float(parts[i + 2])]
                except (ValueError, IndexError):
                    aux_dict[ch_key] = [-1]
                    if tpcont_det:
                        self._temp_dict[1][ch_key] = [-1]

        self._temp_dict[temp_ind - 1] = aux_dict

        if temp_ind == len(self._temp_dict):
            if self._current_setup:
                self._setup_tcal[self._current_setup] = self._temp_dict[-1]
                self._caltemp_read[self._current_setup] = True

    # ----------------------------------------------------------------
    # Header reading
    # ----------------------------------------------------------------

    def _read_header(self, line: str) -> bool:
        """Read header variables (DBBC configuration, format type)."""
        if self.dbbc_mode_name is None:
            ind = self._id_line(line, ["Rack=DBBC", "equip,dbbc_"])
            if ind != 0:
                if ind == 1:
                    aux1 = line.split(" ")[1]
                    aux2 = aux1.split("_")
                    if len(aux2) == 2:
                        self.dbbc_mode_name = aux2[1].strip().upper()
                    else:
                        self.dbbc_mode_name = "DDC"
                elif ind == 2:
                    aux1 = line.split("_")[1]
                    aux2 = aux1.split(",")[0]
                    aux3 = aux2.split("/")[0]
                    self.dbbc_mode_name = aux3.upper()

                if self.dbbc_mode_name == "PFB":
                    self._ch_id = ["a", "b", "c", "d"]
                    self._ch_id_index = 0
                    self._ch_id_set = set(self._ch_id)
                elif self.dbbc_mode_name == "DDC":
                    self._ch_id = ["u", "l"]
                    self._ch_id_index = -1
                    self._ch_id_set = set(self._ch_id)
                else:
                    self.dbbc_mode_name = None

                return True

        if self._id_line(line, ["/form="]):
            self._form_type = line.split("=")[1].rstrip("\n")
            return True

        return False

    # ----------------------------------------------------------------
    # Helpers
    # ----------------------------------------------------------------

    @staticmethod
    def _id_line(line: str, references: list[str]) -> int:
        """Return 1-based index of the first matching reference, or 0."""
        return next((i + 1 for i, ref in enumerate(references) if ref in line), 0)

    @staticmethod
    def _get_datetime(line: str) -> datetime.datetime | None:
        """Parse the timestamp from an FS log line."""
        try:
            parts = line.split(".")
            time_parts = parts[2].split(":")
            usec = int(int(parts[3][:2]) * 1e4)
            year = int(parts[0])
            day = int(parts[1])
            hour = int(time_parts[0])
            minute = int(time_parts[1])
            sec = int(time_parts[2])
            dt = datetime.datetime(year, 1, 1, hour, minute, sec, usec) + datetime.timedelta(days=day - 1)
            return dt
        except (IndexError, ValueError, AttributeError):
            return None

    # ----------------------------------------------------------------
    # Public accessors
    # ----------------------------------------------------------------

    def get_log_data(self) -> list[Any]:
        return self.log_data

    def observation_date(self) -> datetime.date | None:
        """Return the date of the first log timestamp (used for archive directory naming)."""
        for line in self._file_content:
            dt = self._get_datetime(line)
            if dt is not None:
                return dt.date()
        return None

    def lo_array(self) -> dict[str, list[float]]:
        """Unique LO frequencies per setup."""
        flo_raw, _ = self.lo_p_array()
        result: dict[str, list[float]] = {}
        for setup in reversed(list(flo_raw.keys())):
            result[setup] = sorted(
                set(flo_raw[setup]), key=lambda x: flo_raw[setup].index(x)
            )
        return result

    def lo_p_array(self) -> tuple[dict[str, list[float]], dict[str, list[str]]]:
        """LO frequencies and polarizations per setup."""
        flo: dict[str, list[float]] = {}
        pols: dict[str, list[str]] = {}

        for setup in reversed(list(self.pol_array.keys())):
            # dict preserves insertion order (Python 3.7+) and gives O(1)
            # membership test — replaces the former seen: list with O(n) `in`.
            seen: dict[tuple[float, str], None] = {}
            for key in self.pol_array[setup]:
                for freq, pol in zip(self.freq_lo_mhz[setup][key], self.pol_array[setup][key]):
                    seen[(freq, pol)] = None
            flo[setup] = [f for f, _ in seen]
            pols[setup] = [p for _, p in seen]

        return flo, pols

    def station(self) -> str:
        return self.station_name

    def experiment(self) -> str:
        return self.exp_name

    def dbbc_mode(self) -> str | None:
        return self.dbbc_mode_name

    def cal_mode(self) -> dict[str, str]:
        return self.cal_mode_name

    def rxg_files_needed(self) -> dict[str, list[str]]:
        """RXG filenames required for each setup."""
        flo = self.lo_array()
        result: dict[str, list[str]] = {}
        for setup in reversed(list(flo.keys())):
            result[setup] = []
            for freq in flo[setup]:
                result[setup].append(self.get_rxg_filename(freq))
        return result

    def get_rxg_filename(self, freq_lo_mhz: float) -> str:
        """Find the RXG file matching a given LO frequency."""
        st_code = self.station_name[0].upper() + self.station_name[1].lower()

        # Build a flat list of (full_path, is_explicit) pairs across all search dirs
        candidates: list[tuple[str, bool]] = []
        for rxg_dir in self._rxg_dirs:
            if self._rxg_files:
                for fname in self._rxg_files:
                    full_path = os.path.join(rxg_dir, fname)
                    if os.path.isfile(full_path):
                        candidates.append((full_path, True))
            else:
                try:
                    for fname in os.listdir(rxg_dir):
                        if fname.endswith(".rxg") and st_code in fname:
                            candidates.append((os.path.join(rxg_dir, fname), False))
                except FileNotFoundError:
                    continue

        for full_path, _explicit in candidates:
            rxg = RxgFile(full_path)
            try:
                lo_list = rxg.lo()
                if len(lo_list) == 1:
                    if freq_lo_mhz == float(lo_list[0]):
                        return full_path
                else:
                    lo_start = float(lo_list[0])
                    lo_end = float(lo_list[1])
                    if lo_start <= freq_lo_mhz <= lo_end:
                        return full_path
            except Exception as ex:
                print(f"Error getting LO freq: {ex}")
                raise

        return " "

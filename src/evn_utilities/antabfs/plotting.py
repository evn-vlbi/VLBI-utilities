"""Interactive plotting for Tsys data.

A single multi-channel figure lets the user inspect and flag all BBC
channels at once.  Detected outliers are auto-replaced before display.
The user can:

* Click a legend entry to show / hide that channel.
* Drag a rectangle with the left mouse button to define a selection.
* Press **Delete** to flag (replace with local-median estimate) all
  data points inside the selection for every *visible* channel.
"""

from __future__ import annotations

import colorsys
import datetime
import math
import re

import matplotlib

# Try interactive backends in preference order.  Must happen before any
# import of matplotlib.pyplot so the backend is locked in first.
_INTERACTIVE_BACKENDS = ["Qt5Agg", "QtAgg", "TkAgg", "GTK3Agg", "GTK4Agg", "WXAgg"]
_backend_ok = False
for _b in _INTERACTIVE_BACKENDS:
    try:
        matplotlib.use(_b, force=True)
        _backend_ok = True
        break
    except Exception:
        continue

if not _backend_ok:
    matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.widgets import RectangleSelector

from .utils import compute_outliers, modify_data

_UTC = datetime.timezone.utc
_DATE_FMT = mdates.DateFormatter("%j %H:%M")

# Visibility states for legend entries
_ALPHA = {True: 1.0, False: 0.25}

# Special group-toggle labels in the legend
_GROUP_ALL = "all"
_GROUP_LCP = "LCP"
_GROUP_RCP = "RCP"
_GROUPS = (_GROUP_ALL, _GROUP_LCP, _GROUP_RCP)

# Match tokens like L1, L2, R3 anywhere in a label (word boundary + letter + digits)
_LCP_RE = re.compile(r"\bL\d+\b", re.IGNORECASE)
_RCP_RE = re.compile(r"\bR\d+\b", re.IGNORECASE)


def _color(index: int) -> str:
    """Maximally-spread HSV hues: 0 → ½ → ¼ → ¾ → ⅛ …"""
    if index == 0:
        hue = 0.0
    else:
        hb = 1 << (index.bit_length() - 1)
        hue = ((index - hb) * 2 + 1) / (hb * 2)
    r, g, b = colorsys.hsv_to_rgb(hue, 0.9, 0.85)
    return f"#{round(r * 255):02x}{round(g * 255):02x}{round(b * 255):02x}"


def _unix_to_dt(unix_array: np.ndarray) -> list[datetime.datetime]:
    return [datetime.datetime.fromtimestamp(float(t), tz=_UTC) for t in unix_array]


class MultiChannelSelection:
    """Single interactive figure for all BBC channels.

    Parameters
    ----------
    tsys_matrix : ndarray, shape (n_channels, n_times)
        Modified **in-place** when the window closes.
    time : ndarray
        Unix timestamps, shape (n_times,).
    block : ndarray
        Scan-block indices, shape (n_times,).
    bbc_labels : list[str]
        Channel labels (one per row of *tsys_matrix*).
    sigma : float
        MAD-sigma threshold for auto-outlier replacement.
    part_num : int
        0-based setup index shown in the window title.
    """

    def __init__(
        self,
        tsys_matrix: np.ndarray,
        time: np.ndarray,
        block: np.ndarray,
        bbc_labels: list[str],
        sigma: float = 3.0,
        part_num: int = 0,
    ) -> None:
        self.sigma = sigma
        self.block = np.asarray(block)
        self.x = np.asarray(time, dtype=float)
        self.n_ch = tsys_matrix.shape[0]
        self.labels = bbc_labels[: self.n_ch]

        # Save originals before any replacement so we can show flagged points
        # at their original values with reduced opacity.
        self.y_orig: list[np.ndarray] = [tsys_matrix[i].copy() for i in range(self.n_ch)]
        self.y: list[np.ndarray] = [tsys_matrix[i].copy() for i in range(self.n_ch)]
        self.flagged: list[np.ndarray] = [
            np.zeros(len(self.x), dtype=bool) for _ in range(self.n_ch)
        ]

        # --- Auto-replace outliers before display ---
        for ch_i in range(self.n_ch):
            y = self.y[ch_i]
            _, low, up, _, _, _, _ = compute_outliers(self.block, self.x, y, self.sigma)
            out_mask = (y > up) | (y < low)
            self.flagged[ch_i] |= out_mask
            good = ~out_mask
            if out_mask.any() and good.any():
                global_med = float(np.nanmedian(y[good]))
                good_pos = np.where(good)[0]
                for bi in np.where(out_mask)[0]:
                    lo = int(np.searchsorted(good_pos, bi - 20, "left"))
                    hi = int(np.searchsorted(good_pos, bi + 21, "right"))
                    nearby = good_pos[lo:hi]
                    y[bi] = float(np.median(y[nearby])) if len(nearby) else global_med

        # --- Build figure ---
        self._x_dt = _unix_to_dt(self.x)
        self.fig, self.ax = plt.subplots(figsize=(16, 9))
        self.fig.suptitle("Drag to select region  |  Delete = flag selection  "
                          "|  click legend to show/hide")
        self.ax.set_xlabel("Time (UT)")
        self.ax.set_ylabel("Tsys [K]")
        self.ax.grid(True, alpha=0.35)
        self.ax.xaxis.set_major_formatter(_DATE_FMT)

        # Two artists per channel:
        #   _lines[lbl]      – current (replacement) values, full opacity
        #   _lines_orig[lbl] – original values at flagged positions, alpha=0.2
        self._lines: dict[str, plt.Line2D] = {}
        self._lines_orig: dict[str, plt.Line2D] = {}
        self._visible: dict[str, bool] = {}
        for i, lbl in enumerate(self.labels):
            color = _color(i)
            (line,) = self.ax.plot(
                self._x_dt, self.y[i], ".",
                color=color, label=lbl, ms=4, picker=5,
            )
            self._lines[lbl] = line
            self._visible[lbl] = True

            # Flagged-originals artist (no legend entry)
            flag_idx = np.where(self.flagged[i])[0]
            x_flag = [self._x_dt[j] for j in flag_idx]
            y_flag = self.y_orig[i][flag_idx]
            (line_orig,) = self.ax.plot(
                x_flag, y_flag, "x",
                color=color, ms=5, alpha=0.2, zorder=1,
            )
            self._lines_orig[lbl] = line_orig

        # Tight axis limits based on actual data, with 5% padding on each side
        if self._x_dt:
            x_min_num = mdates.date2num(min(self._x_dt))
            x_max_num = mdates.date2num(max(self._x_dt))
            x_pad = (x_max_num - x_min_num) * 0.05
            self.ax.set_xlim(
                mdates.num2date(x_min_num - x_pad),
                mdates.num2date(x_max_num + x_pad),
            )
        all_y_vals = np.concatenate(self.y)
        valid_y = all_y_vals[np.isfinite(all_y_vals) & (all_y_vals > 0)]
        if len(valid_y) > 0:
            y_lo, y_hi = float(valid_y.min()), float(valid_y.max())
            y_pad = (y_hi - y_lo) * 0.05
            self.ax.set_ylim(y_lo - y_pad, y_hi + y_pad)

        # Clickable legend — group toggles first, then per-channel entries
        proxy_all = Line2D([0], [0], linestyle="none", marker="s", color="#444444", ms=7)
        proxy_lcp = Line2D([0], [0], linestyle="none", marker="s", color="#2266cc", ms=7)
        proxy_rcp = Line2D([0], [0], linestyle="none", marker="s", color="#cc2222", ms=7)
        ch_handles, ch_labels = self.ax.get_legend_handles_labels()
        legend = self.ax.legend(
            [proxy_all, proxy_lcp, proxy_rcp] + ch_handles,
            [_GROUP_ALL, _GROUP_LCP, _GROUP_RCP] + ch_labels,
            loc="upper right", fontsize="small", markerscale=2, framealpha=0.85,
        )
        self._leg_handles: dict[str, plt.Artist] = {}
        self._leg_texts: dict[str, plt.Text] = {}
        for handle, text in zip(legend.legend_handles, legend.get_texts()):
            lbl = text.get_text()
            handle.set_picker(5)
            text.set_picker(5)
            self._leg_handles[lbl] = handle
            self._leg_texts[lbl] = text

        self.fig.canvas.mpl_connect("pick_event", self._on_pick)

        # Rectangle selector (left-button drag)
        self._sel_extents: tuple[float, float, float, float] | None = None

        def _on_select(eclick, erelease) -> None:
            self._sel_extents = (
                min(eclick.xdata, erelease.xdata),
                max(eclick.xdata, erelease.xdata),
                min(eclick.ydata, erelease.ydata),
                max(eclick.ydata, erelease.ydata),
            )

        self._rect = RectangleSelector(
            self.ax,
            _on_select,
            useblit=True,
            button=[1],
            minspanx=5,
            minspany=5,
            spancoords="pixels",
            interactive=True,
        )
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)

        self.fig.autofmt_xdate()
        plt.tight_layout()
        plt.show()  # blocks until window is closed

        # --- Write results back into tsys_matrix ---
        for i in range(self.n_ch):
            tsys_matrix[i] = self.y[i]

    # ------------------------------------------------------------------ #
    # Event handlers
    # ------------------------------------------------------------------ #

    def _set_channel_visible(self, lbl: str, visible: bool) -> None:
        self._visible[lbl] = visible
        self._lines[lbl].set_visible(visible)
        self._lines_orig[lbl].set_visible(visible)
        alpha = _ALPHA[visible]
        self._leg_handles[lbl].set_alpha(alpha)
        self._leg_texts[lbl].set_alpha(alpha)

    def _on_pick(self, event) -> None:
        """Toggle channel visibility when a legend handle or text is clicked."""
        artist = event.artist
        lbl = None
        for l, h in self._leg_handles.items():
            if artist is h:
                lbl = l
                break
        if lbl is None:
            for l, t in self._leg_texts.items():
                if artist is t:
                    lbl = l
                    break
        if lbl is None:
            return

        if lbl == _GROUP_ALL:
            targets = self.labels
        elif lbl == _GROUP_LCP:
            targets = [l for l in self.labels if _LCP_RE.search(l)]
        elif lbl == _GROUP_RCP:
            targets = [l for l in self.labels if _RCP_RE.search(l)]
        else:
            # Individual channel toggle
            self._set_channel_visible(lbl, not self._visible[lbl])
            self.fig.canvas.draw_idle()
            return

        # Group toggle: hide all if any visible, else show all
        any_visible = any(self._visible[l] for l in targets)
        new_vis = not any_visible
        for l in targets:
            self._set_channel_visible(l, new_vis)
        self.fig.canvas.draw_idle()

    def _on_key(self, event) -> None:
        """Delete key: flag selected points in all visible channels."""
        if event.key != "delete" or self._sel_extents is None:
            return

        x0_mpl, x1_mpl, ymin, ymax = self._sel_extents
        x0_unix = mdates.num2date(x0_mpl).timestamp()
        x1_unix = mdates.num2date(x1_mpl).timestamp()

        for i, lbl in enumerate(self.labels):
            if not self._visible[lbl]:
                continue

            # Determine which points fall inside the selection rectangle
            # before modifying y, so we can record them as newly flagged.
            in_rect = (
                (self.x >= x0_unix) & (self.x <= x1_unix)
                & (self.y[i] >= ymin) & (self.y[i] <= ymax)
                & (self.y[i] >= 0)
            )
            self.flagged[i] |= in_rect

            self.y[i] = modify_data(
                self.x, self.y[i], self.block,
                x1_unix, x0_unix, ymax, ymin,
            )

            # Update the main (replacement) line
            self._lines[lbl].set_ydata(self.y[i])

            # Update the flagged-originals artist
            flag_idx = np.where(self.flagged[i])[0]
            self._lines_orig[lbl].set_xdata([self._x_dt[j] for j in flag_idx])
            self._lines_orig[lbl].set_ydata(self.y_orig[i][flag_idx])

        self._sel_extents = None
        self.fig.canvas.draw_idle()

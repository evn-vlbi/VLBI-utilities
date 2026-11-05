"""Interactive plotting for Tsys data.

A single multi-channel figure lets the user inspect and flag all BBC
channels at once.  Detected outliers are auto-replaced before display.
The user can:

* Click a legend entry to show / hide that channel.
* Drag a rectangle with the left mouse button to define a selection.
* Click inside the selection or press **d** (or Delete / Backspace) to flag
  all data points inside it for every *visible* channel.
* Press **Ctrl+Z** to undo the last flag / unflag action.
* Press **Escape** to cancel the current selection.
* The y axis follows the *unflagged* data: after a flag / unflag / undo it is re-fitted when a limit moves by
  more than 10% of the current span (e.g. flagging a large spike zooms in on the remaining data).
* Use the control panel next to the legend:
  - **Unflag all**: restore every flagged value (undoable).
  - **Smooth span**: scan / source (same as ``--smooth-span``); "source" smooths across all scans of the same
    source (spline through scan medians), never mixing sources.
  - **Smooth method**: savgol / median / none (same as ``--smooth-method`` / ``--no-smooth``); savgol and median
    only matter for the "scan" span, "none" disables smoothing for both.
  - **Smooth window**: integer >= 0, 0 = the whole scan (same as ``--smooth-window``).
  - **Flag threshold (sigma)**: number > 0 (same as ``--sigma``); re-runs the automatic flagging, keeping the
    manual flags, and is also used for smoothing.  Undoable.
  - **Run smoothing**: replace the data with the smoothed trend of the current *unflagged* values (flags kept).
    Each click smooths the data as it is then, so flag more and click again to smooth the new set.  Undoable.
  The thin line on each channel previews the smoothing (computed from unflagged values only).  The CLI still
  smooths once more with the final settings when the file is written.
* **Close the window** to accept the result and continue; Ctrl+C in the terminal aborts.
"""

from __future__ import annotations

import colorsys
import datetime
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

import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.widgets import Button, RadioButtons, RectangleSelector, TextBox  # noqa: E402

from .smoothing import SMOOTH_SPANS, smooth_tsys_matrix  # noqa: E402
from .utils import compute_outliers, modify_data  # noqa: E402

# Disable LaTeX rendering and mathtext parsing — many target machines do
# not have a LaTeX install and underscores / dollar signs in channel labels
# would otherwise be interpreted as math.
plt.rcParams["text.usetex"] = False
plt.rcParams["text.parse_math"] = False

class _ResizeSafeTextBox(TextBox):
    """TextBox that survives figure resizes on matplotlib 3.11.0.

    matplotlib 3.11.0 wraps ``TextBox._resize`` with ``_call_with_reparented_event``, which reads
    ``event.inaxes``; a ``ResizeEvent`` has no such attribute, so every window resize raised
    ``AttributeError``. This override keeps the upstream behaviour (stop typing) without touching the event.
    Inputs/outputs are identical to ``matplotlib.widgets.TextBox``.
    """

    def _resize(self, event) -> None:
        """Resize callback: end any in-progress text entry. ``event`` (ResizeEvent) is ignored."""
        self.stop_typing()


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

# Smoothing choices shown in the control panel; "none" is the GUI equivalent of --no-smooth.
SMOOTH_METHODS = ("savgol", "median", "none")

_INSTRUCTIONS = ("Define area with mouse.", 'Click or "D" to flag data.', "Ctrl+Z to undo.",
                 '"x" points are flagged values.', "Close window to save file.", "Ctrl+C in the terminal to abort.")

# Figure-coordinate layout: plot | legend | control panel (controls, then calibration info and instructions,
# one item per line).  Invariant: the control panel starts right of the legend (legend is anchored at the
# axes' right edge).
_AXES_RECT = dict(left=0.06, right=0.78, top=0.94, bottom=0.08)
_PANEL_X, _PANEL_W = 0.86, 0.13
_LINE_STEP = 0.022  # vertical spacing (figure fraction) between panel text lines
# After flag / unflag / undo, the y axis is re-fitted only if a limit moves by more than this fraction of the span.
_YLIM_REL_CHANGE = 0.10
# Selection rectangle: small corner/edge handles so they don't cover small boxes; a press/release closer than
# _CLICK_MAX_PX pixels is a click (flag), anything longer is a drag (select / move / resize).
_HANDLE_SIZE = 4
_HANDLE_GRAB_PX = 4
_CLICK_MAX_PX = 3


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


def parse_window(text: str) -> int:
    """Parse a smoothing-window string; return an int >= 0 or raise ValueError."""
    value = int(text.strip())
    if value < 0:
        raise ValueError(f"smooth window must be >= 0, got {value}")
    return value


def parse_sigma(text: str) -> float:
    """Parse a flagging-threshold string; return a finite float > 0 or raise ValueError."""
    value = float(text.strip())
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"sigma must be a positive number, got {value}")
    return value


class MultiChannelSelection:
    """Single interactive figure for all BBC channels.

    Parameters
    ----------
    tsys_matrix : ndarray, shape (n_channels, n_times)
        Modified **in-place** when the window closes (flagged values replaced, *not* smoothed).
    time : ndarray
        Unix timestamps, shape (n_times,).
    block : ndarray
        Scan-block indices, shape (n_times,).
    bbc_labels : list[str]
        Channel labels (one per row of *tsys_matrix*).
    sigma : float
        Initial MAD-sigma threshold for auto-outlier flagging and smoothing (editable in the panel).
    part_num : int
        0-based setup index shown in the window title.
    smooth_method : str
        Initial smoothing method, one of :data:`SMOOTH_METHODS`.
    smooth_window : int
        Initial smoothing window (0 = the whole scan).
    smooth_span : str
        Initial smoothing span, one of ``SMOOTH_SPANS`` (``"scan"`` or ``"source"``).
    scan_sources : dict or None
        ``{scan number: source name}`` used by the ``"source"`` span.
    max_gap : float
        Seconds; same-source scans further apart are smoothed separately (``"source"`` span).
    calib_info : list[str] or None
        Calibration source lines shown in the right-hand panel, one per line (RXG file names or assumed SEFD).
    show : bool
        If False the figure is built but ``plt.show()`` is not called (used by tests).

    After the window closes, ``sigma``, ``smooth_method``, ``smooth_window`` and ``smooth_span`` hold the user's
    final choice.

    Invariant: ``manual[i]`` (user-flagged) is always a subset of ``flagged[i]``; the rest of ``flagged[i]`` are
    automatic sigma-clip flags, which are recomputed from scratch when sigma changes.
    """

    def __init__(self, tsys_matrix: np.ndarray, time: np.ndarray, block: np.ndarray, bbc_labels: list[str],
                 sigma: float = 3.0, part_num: int = 0, smooth_method: str = "savgol", smooth_window: int = 0,
                 smooth_span: str = "scan", scan_sources: dict[int, str] | None = None, max_gap: float = 3600.0,
                 calib_info: list[str] | None = None, show: bool = True) -> None:
        calib_info = calib_info or []
        if smooth_method not in SMOOTH_METHODS:
            raise ValueError(f"smooth_method must be one of {SMOOTH_METHODS}, got {smooth_method!r}")
        if smooth_span not in SMOOTH_SPANS:
            raise ValueError(f"smooth_span must be one of {SMOOTH_SPANS}, got {smooth_span!r}")
        self.sigma = sigma
        self.smooth_method = smooth_method
        self.smooth_window = smooth_window
        self.smooth_span = smooth_span
        self.scan_sources = scan_sources or {}
        self.max_gap = max_gap
        self.block = np.asarray(block)
        self.x = np.asarray(time, dtype=float)
        self.n_ch = tsys_matrix.shape[0]
        self.labels = bbc_labels[: self.n_ch]

        # Save originals before any replacement so we can show flagged points
        # at their original values with reduced opacity.
        self.y_orig: list[np.ndarray] = [tsys_matrix[i].copy() for i in range(self.n_ch)]
        self.y: list[np.ndarray] = [tsys_matrix[i].copy() for i in range(self.n_ch)]
        self.flagged: list[np.ndarray] = [np.zeros(len(self.x), dtype=bool) for _ in range(self.n_ch)]
        self.manual: list[np.ndarray] = [np.zeros(len(self.x), dtype=bool) for _ in range(self.n_ch)]

        # --- Auto-replace outliers before display ---
        for ch_i in range(self.n_ch):
            self._auto_clean_channel(ch_i)

        # --- Build figure ---
        self._x_dt = _unix_to_dt(self.x)
        self._x_num_breaks, self._break_pos = self._scan_break_layout()
        self.fig, self.ax = plt.subplots(figsize=(16, 9))
        self.fig.suptitle(f"Interactive Tsys Editor - Part {part_num + 1}")
        self.ax.set_xlabel("Time (UT)")
        self.ax.set_ylabel("Tsys [K]")
        self.ax.grid(True, alpha=0.35)
        self.ax.xaxis.set_major_formatter(_DATE_FMT)

        self._build_lines()
        self._set_axis_limits()
        self._build_legend()
        self._build_controls()
        self._build_info_text(calib_info)
        self._update_smooth_lines()

        # Undo stack: each entry is a snapshot of (y_copies, flagged_copies, manual_copies, sigma)
        self._undo_stack: list[tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray], float]] = []

        # Rectangle selector (left-button drag).  _sel_extents = (x0, x1, y0, y1) in data coords, or None.
        self._sel_extents: tuple[float, float, float, float] | None = None
        # Pixel position of the last left-button press, used to tell a click (flag) from a drag (select/move).
        self._press_px: tuple[float, float] | None = None
        self._sel_at_press: tuple[float, float, float, float] | None = None

        def _on_select(_eclick, _erelease) -> None:
            # Read the rectangle itself: correct for a new drag, a move and a handle resize alike.
            # (eclick/erelease are the mouse positions, which for a move/resize are NOT the rectangle corners.)
            self._sel_extents = tuple(float(v) for v in self._rect.extents)

        self._rect = RectangleSelector(self.ax, _on_select, useblit=True, button=[1], minspanx=5, minspany=5,
                                       spancoords="pixels", interactive=True, grab_range=_HANDLE_GRAB_PX,
                                       handle_props=dict(markersize=_HANDLE_SIZE))
        self.fig.canvas.mpl_connect("pick_event", self._on_pick)
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self.fig.canvas.mpl_connect("button_press_event", self._on_button_press)
        self.fig.canvas.mpl_connect("button_release_event", self._on_button_release)

        self.fig.autofmt_xdate()
        # Reserve right-hand strip for the legend and control panel instead of tight_layout.
        plt.subplots_adjust(**_AXES_RECT)
        if not show:
            return
        plt.show()  # blocks until window is closed

        # --- Write results back into tsys_matrix ---
        for i in range(self.n_ch):
            tsys_matrix[i] = self.y[i]

    # ------------------------------------------------------------------ #
    # Figure construction
    # ------------------------------------------------------------------ #

    def _scan_break_layout(self) -> tuple[np.ndarray, np.ndarray]:
        """Return (x as matplotlib date numbers with NaN at scan boundaries, insertion positions).

        NaN breaks keep the smoothed preview line from connecting different scans.
        """
        break_pos = np.where(np.diff(self.block) != 0)[0] + 1
        x_num = np.asarray(mdates.date2num(self._x_dt), dtype=float) if self._x_dt else np.array([])
        return np.insert(x_num, break_pos, np.nan), break_pos

    def _build_lines(self) -> None:
        """Create the data, flagged-original ("x") and smoothed-preview artists for every channel."""
        # Three artists per channel:
        #   _lines[lbl]        – current (replacement) values, full opacity
        #   _lines_orig[lbl]   – original values at flagged positions ("x"), alpha=0.2
        #   _lines_smooth[lbl] – smoothed preview for the current method/window
        self._lines: dict[str, plt.Line2D] = {}
        self._lines_orig: dict[str, plt.Line2D] = {}
        self._lines_smooth: dict[str, plt.Line2D] = {}
        self._visible: dict[str, bool] = {}
        for i, lbl in enumerate(self.labels):
            color = _color(i)
            (line,) = self.ax.plot(self._x_dt, self.y[i], ".", color=color, label=lbl, ms=4, picker=5)
            self._lines[lbl] = line
            self._visible[lbl] = True

            flag_idx = np.where(self.flagged[i])[0]
            (line_orig,) = self.ax.plot([self._x_dt[j] for j in flag_idx], self.y_orig[i][flag_idx], "x",
                                        color=color, ms=5, alpha=0.2, zorder=1)
            self._lines_orig[lbl] = line_orig

            (line_smooth,) = self.ax.plot([], [], "-", color=color, lw=1.2, alpha=0.9, zorder=3)
            self._lines_smooth[lbl] = line_smooth

    def _set_axis_limits(self) -> None:
        """Tight axis limits based on actual data, with 5% padding on each side (y from unflagged data only)."""
        if self._x_dt:
            x_min_num = mdates.date2num(min(self._x_dt))
            x_max_num = mdates.date2num(max(self._x_dt))
            x_pad = (x_max_num - x_min_num) * 0.05
            self.ax.set_xlim(mdates.num2date(x_min_num - x_pad), mdates.num2date(x_max_num + x_pad))
        ylim = self._data_ylim()
        if ylim is not None:
            self.ax.set_ylim(*ylim)

    def _data_ylim(self) -> tuple[float, float] | None:
        """Return padded (low, high) y limits spanning the unflagged, valid Tsys of all channels, or None if empty.

        Flagged values are excluded, so flagging a spike lets the axis shrink to the remaining data.
        """
        if not self.y:
            return None
        vals = np.concatenate([y[~f] for y, f in zip(self.y, self.flagged)])
        vals = vals[np.isfinite(vals) & (vals > 0)]
        if len(vals) == 0:
            return None
        y_lo, y_hi = float(vals.min()), float(vals.max())
        y_pad = max((y_hi - y_lo) * 0.05, 0.5)
        return y_lo - y_pad, y_hi + y_pad

    def _update_ylim(self) -> None:
        """Re-fit the y axis to the unflagged data if either limit moved by more than _YLIM_REL_CHANGE of the span.

        Small changes are ignored so the axis does not jitter on every flag. The x axis is never touched.
        """
        new = self._data_ylim()
        if new is None:
            return
        cur_lo, cur_hi = self.ax.get_ylim()
        span = max(cur_hi - cur_lo, 1e-9)
        if abs(new[0] - cur_lo) > _YLIM_REL_CHANGE * span or abs(new[1] - cur_hi) > _YLIM_REL_CHANGE * span:
            print(f"  Editor: y axis rescaled to {new[0]:.1f} .. {new[1]:.1f} K")
            self.ax.set_ylim(*new)

    def _build_legend(self) -> None:
        """Clickable legend placed OUTSIDE the axes.

        Outside so the persistent RectangleSelector (interactive=True) handles cannot intercept clicks on it.
        """
        proxy_all = Line2D([0], [0], linestyle="none", marker="s", color="#444444", ms=7)
        proxy_lcp = Line2D([0], [0], linestyle="none", marker="s", color="#2266cc", ms=7)
        proxy_rcp = Line2D([0], [0], linestyle="none", marker="s", color="#cc2222", ms=7)
        ch_handles, ch_labels = self.ax.get_legend_handles_labels()
        legend = self.ax.legend([proxy_all, proxy_lcp, proxy_rcp] + ch_handles, list(_GROUPS) + ch_labels,
                                loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize="small", markerscale=2,
                                framealpha=0.85, borderaxespad=0.0)
        self._leg_handles: dict[str, plt.Artist] = {}
        self._leg_texts: dict[str, plt.Text] = {}
        for handle, text in zip(legend.legend_handles, legend.get_texts()):
            lbl = text.get_text()
            handle.set_picker(5)
            text.set_picker(5)
            self._leg_handles[lbl] = handle
            self._leg_texts[lbl] = text

    def _build_controls(self) -> None:
        """Control panel right of the legend, top to bottom (figure fractions, running cursor *y*).

        Unflag all | Smooth span | Smooth method | Smooth window | Flag threshold | Run smoothing | status line.
        Sets ``self._panel_bottom`` (where the info text starts).
        """
        y = _AXES_RECT["top"]

        def label(text: str) -> None:
            nonlocal y
            y -= 0.028
            self.fig.text(_PANEL_X, y, text, fontsize=10, fontweight="bold")

        def hint(text: str) -> None:
            nonlocal y
            y -= 0.022
            self.fig.text(_PANEL_X, y, text, fontsize=8, color="#555555")

        def widget_axes(height: float, gap: float = 0.006):
            nonlocal y
            y -= height + gap
            return self.fig.add_axes((_PANEL_X, y, _PANEL_W, height))

        self._btn_unflag = Button(widget_axes(0.045, gap=0.0), "Unflag all")
        self._btn_unflag.on_clicked(lambda _event: self._unflag_all())

        label("Smooth span")
        self._radio_span = RadioButtons(widget_axes(0.065), SMOOTH_SPANS, active=SMOOTH_SPANS.index(self.smooth_span))
        self._radio_span.on_clicked(self._on_span)
        hint("source: across same-source scans.")

        label("Smooth method")
        self._radio_method = RadioButtons(widget_axes(0.095), SMOOTH_METHODS,
                                          active=SMOOTH_METHODS.index(self.smooth_method))
        self._radio_method.on_clicked(self._on_method)
        hint("savgol/median: scan span only.")

        label("Smooth window (0=scan)")
        self._box_window = _ResizeSafeTextBox(widget_axes(0.04), "", initial=str(self.smooth_window))
        self._box_window.on_submit(self._on_window)
        hint("Enter to apply; used by median.")

        label("Flag threshold (sigma)")
        self._box_sigma = _ResizeSafeTextBox(widget_axes(0.04), "", initial=f"{self.sigma:g}")
        self._box_sigma.on_submit(self._on_sigma)
        hint("Enter to apply; re-runs auto-flag.")

        y -= 0.01
        self._btn_smooth = Button(widget_axes(0.045), "Run smoothing")
        self._btn_smooth.on_clicked(lambda _event: self._run_smoothing())

        y -= 0.025
        self._status = self.fig.text(_PANEL_X, y, "", fontsize=9, color="#b00020", wrap=True)
        self._panel_bottom = y - 0.012

    def _build_info_text(self, calib_info: list[str]) -> None:
        """Panel text below the controls: calibration source lines, then fixed instructions, one item per line."""
        y = self._panel_bottom
        for heading, lines, color in (("Calibration", calib_info, "#333333"),
                                      ("Instructions", _INSTRUCTIONS, "#1a8a2e")):
            self.fig.text(_PANEL_X, y, heading, fontsize=10, fontweight="bold", va="top")
            y -= _LINE_STEP
            for line in lines:
                self.fig.text(_PANEL_X, y, line, fontsize=9, color=color, va="top")
                y -= _LINE_STEP
            y -= _LINE_STEP / 2

    # ------------------------------------------------------------------ #
    # Control-panel callbacks
    # ------------------------------------------------------------------ #

    def _on_span(self, label: str) -> None:
        """Radio callback: set the smoothing span (scan / source) and refresh the preview."""
        self.smooth_span = label
        print(f"  Editor: smooth span = {label}")
        self._update_smooth_lines()
        self.fig.canvas.draw_idle()

    def _on_method(self, label: str) -> None:
        """Radio callback: set the smoothing method and refresh the preview."""
        self.smooth_method = label
        print(f"  Editor: smooth method = {label}")
        self._update_smooth_lines()
        self.fig.canvas.draw_idle()

    def _on_window(self, text: str) -> None:
        """TextBox callback: validate and set the smoothing window, or show an error and restore the old value."""
        try:
            value = parse_window(text)
        except ValueError:
            # set_val re-triggers on_submit with the old valid value (which clears the status), so set status after.
            self._box_window.set_val(str(self.smooth_window))
            self._status.set_text(f"Invalid window {text!r}: use an integer >= 0.")
            self.fig.canvas.draw_idle()
            return
        self._status.set_text("")
        if value != self.smooth_window:
            self.smooth_window = value
            print(f"  Editor: smooth window = {value}")
            self._update_smooth_lines()
        self.fig.canvas.draw_idle()

    def _on_sigma(self, text: str) -> None:
        """TextBox callback: set the flag threshold and redo the automatic flags (manual flags kept). Undoable."""
        try:
            value = parse_sigma(text)
        except ValueError:
            # set_val re-triggers on_submit with the old valid value (which clears the status), so set status after.
            self._box_sigma.set_val(f"{self.sigma:g}")
            self._status.set_text(f"Invalid sigma {text!r}: use a number > 0.")
            self.fig.canvas.draw_idle()
            return
        self._status.set_text("")
        if value != self.sigma:
            self._push_undo()
            self.sigma = value
            self._reflag_with_current_sigma()
            n_auto = int(sum((f & ~m).sum() for f, m in zip(self.flagged, self.manual)))
            print(f"  Editor: sigma = {value:g} -> {n_auto} values auto-flagged")
            self._refresh_all_lines()
        self.fig.canvas.draw_idle()

    def _reflag_with_current_sigma(self) -> None:
        """Rebuild y/flagged from the originals: keep manual flags, then re-run sigma-clipping with self.sigma.

        Manually flagged values are replaced by the local median of unflagged neighbours (as auto-flags are).
        """
        self.y = [arr.copy() for arr in self.y_orig]
        self.flagged = [m.copy() for m in self.manual]
        for ch_i in range(self.n_ch):
            self._replace_with_local_median(ch_i, self.manual[ch_i])
            self._auto_clean_channel(ch_i)

    def _unflag_all(self) -> None:
        """Restore every value to its original (pre-flag) state, including auto-flagged outliers. Undoable."""
        if not any(f.any() for f in self.flagged):
            return
        self._push_undo()
        n_flagged = int(sum(f.sum() for f in self.flagged))
        self.y = [arr.copy() for arr in self.y_orig]
        self.flagged = [np.zeros(len(self.x), dtype=bool) for _ in range(self.n_ch)]
        self.manual = [np.zeros(len(self.x), dtype=bool) for _ in range(self.n_ch)]
        print(f"  Editor: unflagged {n_flagged} values")
        self._refresh_all_lines()
        self.fig.canvas.draw_idle()

    def _compute_smoothed(self) -> np.ndarray:
        """Smooth every channel with the current method/window/sigma using only the unflagged values.

        Flagged values are passed as NaN so the smoother ignores them, then the fitted trend is evaluated at every
        sample.  Where no fit is possible (e.g. a scan with every value flagged) the current value is kept.
        Returns an array of shape (n_channels, n_times).  Must not be called with method "none".
        """
        current = np.vstack(self.y)
        unflagged = np.where(np.vstack(self.flagged), np.nan, current)
        smoothed = smooth_tsys_matrix(unflagged, times=self.x, block=self.block, sigma=self.sigma,
                                      window=self.smooth_window, method=self.smooth_method, span=self.smooth_span,
                                      scan_sources=self.scan_sources, max_gap=self.max_gap)
        return np.where(np.isfinite(smoothed), smoothed, current)

    def _run_smoothing(self) -> None:
        """Button callback: replace the data with the smoothed trend of the current unflagged values. Undoable.

        Flags are kept (the "x" markers stay); each click smooths the data as it is at that moment, so flagging more
        points and clicking again smooths the new set.
        """
        if self.smooth_method == "none":
            self._status.set_text('Smooth method is "none": nothing to run.')
            self.fig.canvas.draw_idle()
            return
        self._push_undo()
        smoothed = self._compute_smoothed()
        self.y = [smoothed[i].copy() for i in range(self.n_ch)]
        self._status.set_text("")
        print(f"  Editor: smoothing applied (span = {self.smooth_span}, method = {self.smooth_method}, "
              f"window = {self.smooth_window}, sigma = {self.sigma:g})")
        self._refresh_all_lines()
        self.fig.canvas.draw_idle()

    def _update_smooth_lines(self) -> None:
        """Recompute the smoothed preview for all channels with the current method/window."""
        if self.smooth_method == "none" or self.n_ch == 0 or len(self.x) == 0:
            for line in self._lines_smooth.values():
                line.set_data([], [])
            return
        smoothed = self._compute_smoothed()
        for i, lbl in enumerate(self.labels):
            self._lines_smooth[lbl].set_data(self._x_num_breaks, np.insert(smoothed[i], self._break_pos, np.nan))

    # ------------------------------------------------------------------ #
    # Event handlers
    # ------------------------------------------------------------------ #

    def _set_channel_visible(self, lbl: str, visible: bool) -> None:
        self._visible[lbl] = visible
        self._lines[lbl].set_visible(visible)
        self._lines_orig[lbl].set_visible(visible)
        self._lines_smooth[lbl].set_visible(visible)
        alpha = _ALPHA[visible]
        self._leg_handles[lbl].set_alpha(alpha)
        self._leg_texts[lbl].set_alpha(alpha)

    def _on_pick(self, event) -> None:
        """Toggle channel visibility when a legend handle or text is clicked."""
        artist = event.artist
        lbl = next((k for k, h in self._leg_handles.items() if artist is h), None)
        if lbl is None:
            lbl = next((k for k, t in self._leg_texts.items() if artist is t), None)
        if lbl is None:
            return

        if lbl == _GROUP_ALL:
            targets = self.labels
        elif lbl == _GROUP_LCP:
            targets = [lb for lb in self.labels if _LCP_RE.search(lb)]
        elif lbl == _GROUP_RCP:
            targets = [lb for lb in self.labels if _RCP_RE.search(lb)]
        else:
            # Individual channel toggle
            self._set_channel_visible(lbl, not self._visible[lbl])
            self.fig.canvas.draw_idle()
            return

        # Group toggle: hide all if any visible, else show all
        new_vis = not any(self._visible[lb] for lb in targets)
        for lb in targets:
            self._set_channel_visible(lb, new_vis)
        self.fig.canvas.draw_idle()

    def _push_undo(self) -> None:
        """Save current y, flagged, manual and sigma state onto the undo stack."""
        self._undo_stack.append(([arr.copy() for arr in self.y], [arr.copy() for arr in self.flagged],
                                 [arr.copy() for arr in self.manual], self.sigma))

    def _pop_undo(self) -> None:
        """Restore the last saved state from the undo stack (including sigma and its text box)."""
        if not self._undo_stack:
            return
        self.y, self.flagged, self.manual, self.sigma = self._undo_stack.pop()
        self._box_sigma.set_val(f"{self.sigma:g}")  # on_submit sees an unchanged sigma -> no re-flag
        self._refresh_all_lines()
        self._clear_selection()
        self.fig.canvas.draw_idle()

    def _clear_selection(self) -> None:
        """Remove the selection rectangle from the plot and forget its extents.

        Needed after every flag / undo / escape: a leftover rectangle would otherwise stay on screen and, once
        the y axis is rescaled, could cover most of the plot so that new presses move it instead of drawing a
        new selection.
        """
        self._sel_extents = None
        self._rect.clear()

    def _refresh_all_lines(self) -> None:
        """Update all line artists (data, flagged "x", smoothed preview) from current y / flagged arrays."""
        for i, lbl in enumerate(self.labels):
            self._lines[lbl].set_ydata(self.y[i])
            flag_idx = np.where(self.flagged[i])[0]
            self._lines_orig[lbl].set_xdata([self._x_dt[j] for j in flag_idx])
            self._lines_orig[lbl].set_ydata(self.y_orig[i][flag_idx])
        self._update_smooth_lines()
        self._update_ylim()

    def _on_key(self, event) -> None:
        """Handle keyboard shortcuts (ignored while typing in a panel text box).

        d / delete / backspace : flag selected points in visible channels.
        ctrl+z                : undo last flag / unflag / sigma change.
        escape                : cancel current selection.
        """
        if any(getattr(box, "capturekeystrokes", False) for box in (self._box_window, self._box_sigma)):
            return
        if event.key in ("d", "D", "delete", "backspace"):
            self._apply_flag()
        elif event.key == "ctrl+z":
            self._pop_undo()
        elif event.key == "escape":
            self._clear_selection()
            self.fig.canvas.draw_idle()

    def _on_button_press(self, event) -> None:
        """Remember where a left-button press happened and which selection existed at that moment."""
        if event.button == 1 and event.inaxes is self.ax:
            self._press_px = (event.x, event.y)
            self._sel_at_press = self._sel_extents

    def _on_button_release(self, event) -> None:
        """Flag when the user clicks (press and release < _CLICK_MAX_PX apart) inside the current selection.

        A click outside the selection cancels it.  Drags (select / move / resize) are handled by the
        RectangleSelector and its ``_on_select`` callback, not here.
        """
        press, self._press_px = self._press_px, None
        if event.button != 1 or press is None or event.xdata is None or event.ydata is None:
            return
        if np.hypot(event.x - press[0], event.y - press[1]) > _CLICK_MAX_PX:
            return
        sel = self._sel_at_press
        if sel is not None and sel[0] <= event.xdata <= sel[1] and sel[2] <= event.ydata <= sel[3]:
            self._sel_extents = sel
            self._apply_flag()
        else:
            self._clear_selection()
            self.fig.canvas.draw_idle()

    def _apply_flag(self) -> None:
        """Flag selected points in all visible channels, then re-clean outliers.

        After manually flagging the obvious bad points the statistics on the
        remaining good points are tighter, so we re-run the sigma-clip
        outlier detection on each affected channel and replace any newly
        identified outliers as well.
        """
        if self._sel_extents is None:
            return

        self._push_undo()

        x0_mpl, x1_mpl, ymin, ymax = self._sel_extents
        x0_unix = mdates.num2date(x0_mpl).timestamp()
        x1_unix = mdates.num2date(x1_mpl).timestamp()

        for i, lbl in enumerate(self.labels):
            if not self._visible[lbl]:
                continue
            in_rect = ((self.x >= x0_unix) & (self.x <= x1_unix) & (self.y[i] >= ymin) & (self.y[i] <= ymax)
                       & (self.y[i] >= 0) & ~self.flagged[i])   # only newly-flagged points
            if not in_rect.any():
                continue
            self.flagged[i] |= in_rect
            self.manual[i] |= in_rect
            self.y[i] = modify_data(self.x, self.y[i], self.block, x1_unix, x0_unix, ymax, ymin)
            # Re-run outlier detection on the remaining unflagged points;
            # the tighter statistics may catch borderline points that were
            # previously hidden by the now-flagged extreme outliers.
            self._auto_clean_channel(i)

        self._refresh_all_lines()
        self._clear_selection()
        self.fig.canvas.draw_idle()

    def _auto_clean_channel(self, ch_i: int) -> None:
        """Auto-detect and replace outliers on a single channel.

        Operates on currently-unflagged points only.  Newly-detected
        outliers are added to ``self.flagged[ch_i]`` and their values in
        ``self.y[ch_i]`` replaced with a local-median estimate.
        """
        y = self.y[ch_i]
        _, low, up, _, _, _, _ = compute_outliers(self.block, self.x, y, self.sigma)
        out_mask = ((y > up) | (y < low)) & ~self.flagged[ch_i]
        if not out_mask.any():
            return
        self.flagged[ch_i] |= out_mask
        self._replace_with_local_median(ch_i, out_mask)

    def _replace_with_local_median(self, ch_i: int, mask: np.ndarray) -> None:
        """Replace ``self.y[ch_i][mask]`` in place with the median of unflagged values within ±20 samples.

        Falls back to the channel's global unflagged median; does nothing if no value is unflagged.
        """
        y = self.y[ch_i]
        good = ~self.flagged[ch_i]
        if not mask.any() or not good.any():
            return
        global_med = float(np.nanmedian(y[good]))
        good_pos = np.where(good)[0]
        for bi in np.where(mask)[0]:
            lo = int(np.searchsorted(good_pos, bi - 20, "left"))
            hi = int(np.searchsorted(good_pos, bi + 21, "right"))
            nearby = good_pos[lo:hi]
            y[bi] = float(np.median(y[nearby])) if len(nearby) else global_med

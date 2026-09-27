from __future__ import annotations

from collections.abc import Sequence
import math

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QColor, QLinearGradient, QPainter, QPaintEvent, QPen
from PySide6.QtWidgets import (
    QDoubleSpinBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from plasic_app.protocol.pstream import Frame
from plasic_app.ui.theme import ACCENT, PLOT_BACKGROUND, TEXT, TEXT_MUTED
from plasic_app.ui.widgets import SegmentedControl
from plasic_app.visualization.coastlines import flattened_coastlines
from plasic_app.visualization.data import (
    RegionBounds,
    finite_color_range,
    frame_time_label,
    gaussian_grid,
    region_subset,
    resample_latitude_for_image,
)
from plasic_app.visualization.globe_view import ClimateGlobeWidget
from plasic_app.visualization.slice_view import LayerSliceView


SEA_ICE_VARIABLES = frozenset({"sea_ice_fraction", "sea_ice_thickness"})


class ColorScaleWidget(QWidget):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedHeight(48)
        self._lookup = np.asarray(pg.colormap.get("CET-D1").getLookupTable(nPts=32), dtype=np.uint8)
        self._minimum, self._maximum, self._unit = 0.0, 1.0, ""

    def set_colormap(self, name: str) -> None:
        self._lookup = np.asarray(pg.colormap.get(name).getLookupTable(nPts=32), dtype=np.uint8)
        self.update()

    def set_range(self, minimum: float, maximum: float, unit: str) -> None:
        self._minimum, self._maximum, self._unit = minimum, maximum, unit
        self.update()

    @staticmethod
    def _format(value: float) -> str:
        return f"{value:.3e}" if value and (abs(value) < 0.01 or abs(value) >= 1.0e5) else f"{value:.3f}"

    def paintEvent(self, _event: QPaintEvent) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        bar = QRectF(18, 3, max(40, self.width() - 36), 12)
        gradient = QLinearGradient(bar.left(), 0, bar.right(), 0)
        for index, color in enumerate(self._lookup):
            gradient.setColorAt(index / max(len(self._lookup) - 1, 1), QColor(*map(int, color[:3])))
        painter.setBrush(gradient)
        painter.setPen(QPen(QColor("#CBD3D9"), 1))
        painter.drawRoundedRect(bar, 5, 5)
        painter.setPen(QColor(TEXT_MUTED))
        text = QRectF(18, 19, self.width() - 36, 24)
        painter.drawText(text, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, self._format(self._minimum))
        painter.drawText(text, Qt.AlignmentFlag.AlignCenter, self._unit)
        painter.drawText(text, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, self._format(self._maximum))


class LiveFieldView(QWidget):
    """Coordinates global, regional/3-D, and global-mean visualizations."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._frames: tuple[Frame, ...] = ()
        self._variable_definition: dict[str, object] = {}
        self._frame: Frame | None = None
        self._level = 0
        self._view_mode = "map"
        self._primary_mode = "global"
        self._region_view = "2d"
        self._colormap = "CET-D1"
        self._color_mode = "dynamic"
        self._manual_levels: tuple[float, float] | None = None
        self._fixed_cache: dict[tuple[object, ...], tuple[float, float]] = {}
        # 无帧可显示时的状态文本；show_unavailable() 会把它换成具体原因，
        # 使模式切换等操作重绘时不会丢失解释。
        self._empty_state_message = "Waiting for model output…"
        self.max_points = 128
        self.steps: list[int] = []
        self.means: list[float] = []
        self.current_name: str | None = None
        self.current_level: int | None = None

        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(7)

        mode_row = QHBoxLayout()
        mode_row.addWidget(QLabel("Mode"))
        self.primary_control = SegmentedControl([
            ("Global field", "global"),
            ("Region & 3-D", "region"),
            ("Global mean", "series"),
        ])
        self.primary_control.set_option_tooltip(
            "series", "Requires at least two model-time frames for this variable."
        )
        mode_row.addWidget(self.primary_control, 1)
        self.status_label = QLabel("Waiting for model output…")
        self.status_label.setProperty("role", "secondary")
        self.status_label.setWordWrap(True)
        mode_row.addWidget(self.status_label, 2)
        root.addLayout(mode_row)

        self.region_controls = QWidget()
        region_row = QHBoxLayout(self.region_controls)
        region_row.setContentsMargins(0, 0, 0, 0)
        region_row.setSpacing(5)
        self.region_view_control = SegmentedControl([("2-D field", "2d"), ("3-D structure", "3d")])
        self.region_view_control.set_option_tooltip(
            "3d", "Requires a catalogued vertical coordinate and at least two runtime layers."
        )
        region_row.addWidget(self.region_view_control)
        self.bound_spins: dict[str, QDoubleSpinBox] = {}
        for label, key, value in (("W", "west", -60.0), ("E", "east", 60.0), ("S", "south", -30.0), ("N", "north", 30.0)):
            region_row.addWidget(QLabel(label))
            spin = QDoubleSpinBox()
            spin.setRange(-180.0 if key in {"west", "east"} else -90.0, 180.0 if key in {"west", "east"} else 90.0)
            spin.setDecimals(1)
            spin.setSuffix("°")
            spin.setValue(value)
            spin.setMaximumWidth(88)
            self.bound_spins[key] = spin
            region_row.addWidget(spin)
        self.apply_region_button = QPushButton("Apply region")
        region_row.addWidget(self.apply_region_button)
        self.layers_label = QLabel("Layers")
        region_row.addWidget(self.layers_label)
        self.layer_from = QSpinBox()
        self.layer_to = QSpinBox()
        self.layer_from.setPrefix("from ")
        self.layer_to.setPrefix("to ")
        region_row.addWidget(self.layer_from)
        region_row.addWidget(self.layer_to)
        self.slice_count = QSpinBox()
        self.slice_count.setRange(2, 12)
        self.slice_count.setValue(8)
        self.slice_count.setPrefix("show ")
        region_row.addWidget(self.slice_count)
        root.addWidget(self.region_controls)
        self.region_controls.hide()

        self.stack = QStackedWidget()
        root.addWidget(self.stack, 1)

        self.global_stack = QStackedWidget()
        self.map_plot, self.image, self.map_coasts = self._make_map()
        self.globe = ClimateGlobeWidget()
        self.global_stack.addWidget(self.map_plot)
        self.global_stack.addWidget(self.globe)
        self.stack.addWidget(self.global_stack)

        self.region_stack = QStackedWidget()
        self.region_plot, self.region_image, self.region_coasts = self._make_map()
        self.slice_view = LayerSliceView()
        self.region_stack.addWidget(self.region_plot)
        self.region_stack.addWidget(self.slice_view)
        self.stack.addWidget(self.region_stack)

        self.series_page = QWidget()
        series_layout = QVBoxLayout(self.series_page)
        series_layout.setContentsMargins(0, 0, 0, 0)
        self.series_plot = pg.PlotWidget(axisItems={"bottom": pg.AxisItem(orientation="bottom")})
        self._style_plot(self.series_plot)
        self.series_plot.setLabel("bottom", "Model time")
        self.series_plot.setLabel("left", "Global mean")
        self.series_plot.showGrid(x=True, y=True, alpha=0.12)
        self.series_curve = self.series_plot.plot(pen=pg.mkPen(color=ACCENT, width=2.2), symbol="o", symbolSize=5)
        self.coverage_label = QLabel()
        self.coverage_label.setProperty("role", "secondary")
        self.coverage_label.setWordWrap(True)
        series_layout.addWidget(self.series_plot, 1)
        series_layout.addWidget(self.coverage_label)
        self.stack.addWidget(self.series_page)

        self.color_controls = QWidget()
        color_layout = QVBoxLayout(self.color_controls)
        color_layout.setContentsMargins(0, 0, 0, 0)
        color_row = QHBoxLayout()
        color_row.addWidget(QLabel("Colorbar"))
        self.color_mode_control = SegmentedControl([("Fixed", "fixed"), ("Dynamic", "dynamic")])
        self.color_mode_control.set_value("dynamic")
        color_row.addWidget(self.color_mode_control)
        self.vmin_edit = QLineEdit()
        self.vmax_edit = QLineEdit()
        self.vmin_edit.setPlaceholderText("vmin")
        self.vmax_edit.setPlaceholderText("vmax")
        self.vmin_edit.setMaximumWidth(105)
        self.vmax_edit.setMaximumWidth(105)
        color_row.addWidget(self.vmin_edit)
        color_row.addWidget(self.vmax_edit)
        self.apply_levels_button = QPushButton("Apply")
        self.auto_levels_button = QPushButton("Use automatic")
        color_row.addWidget(self.apply_levels_button)
        color_row.addWidget(self.auto_levels_button)
        self.color_note = QLabel()
        self.color_note.setProperty("role", "secondary")
        color_row.addWidget(self.color_note, 1)
        color_layout.addLayout(color_row)
        self.color_scale = ColorScaleWidget()
        color_layout.addWidget(self.color_scale)
        root.addWidget(self.color_controls)

        self.primary_control.valueChanged.connect(self.set_primary_mode)
        self.region_view_control.valueChanged.connect(self._set_region_view)
        self.color_mode_control.valueChanged.connect(self.set_color_mode)
        self.apply_region_button.clicked.connect(self._render)
        self.apply_levels_button.clicked.connect(self._apply_manual_levels)
        self.auto_levels_button.clicked.connect(self._clear_manual_levels)
        self.layer_from.valueChanged.connect(self._vertical_range_changed)
        self.layer_to.valueChanged.connect(self._vertical_range_changed)
        self.slice_count.valueChanged.connect(self._render)
        self.slice_view.renderModeChanged.connect(self._volume_render_mode_changed)
        self._sync_capabilities()
        self.set_color_mode("dynamic")

    @staticmethod
    def _style_plot(plot: pg.PlotWidget) -> None:
        plot.setBackground(PLOT_BACKGROUND)
        plot.hideButtons()
        plot.setMenuEnabled(False)
        for name in ("left", "bottom"):
            plot.getAxis(name).setPen(pg.mkPen(TEXT_MUTED))
            plot.getAxis(name).setTextPen(pg.mkPen(TEXT_MUTED))

    def _make_map(self) -> tuple[pg.PlotWidget, pg.ImageItem, tuple[pg.PlotDataItem, pg.PlotDataItem]]:
        plot = pg.PlotWidget()
        self._style_plot(plot)
        plot.setLabel("bottom", "Longitude", units="deg")
        plot.setLabel("left", "Latitude", units="deg")
        plot.setAspectLocked(True, ratio=1.0)
        image = pg.ImageItem(axisOrder="row-major")
        plot.addItem(image)
        x, y = flattened_coastlines()
        halo = pg.PlotDataItem(x, y, connect="finite", pen=pg.mkPen((255, 255, 255, 220), width=2.8))
        coast = pg.PlotDataItem(x, y, connect="finite", pen=pg.mkPen((32, 43, 52, 205), width=0.9))
        halo.setZValue(10)
        coast.setZValue(11)
        plot.addItem(halo)
        plot.addItem(coast)
        return plot, image, (halo, coast)

    @property
    def view_mode(self) -> str:
        return self._view_mode

    @property
    def primary_mode(self) -> str:
        return self._primary_mode

    @property
    def auto_levels(self) -> bool:
        return self._color_mode == "dynamic"

    def set_history_limit(self, max_points: int) -> None:
        if max_points < 2:
            raise ValueError("history limit must be at least two")
        self.max_points = max_points

    def set_frames(self, frames: Sequence[Frame]) -> None:
        frames = tuple(frames)[-self.max_points :]
        old_name = self._frames[-1].name if self._frames else None
        new_name = frames[-1].name if frames else None
        if old_name != new_name:
            self._manual_levels = None
            self._fixed_cache.clear()
            self.vmin_edit.clear()
            self.vmax_edit.clear()
        self._frames = frames
        if self._frame is not None:
            matches = [item for item in frames if item.step == self._frame.step]
            self._frame = matches[-1] if matches else (frames[-1] if frames else None)
        elif frames:
            self._frame = frames[-1]
        self._sync_capabilities()

    def reset_data_context(self) -> None:
        """Forget data-source-dependent state before a new run or replay."""
        self._frames = ()
        self._frame = None
        self._manual_levels = None
        self._fixed_cache.clear()
        self.vmin_edit.clear()
        self.vmax_edit.clear()
        self.series_curve.setData([], [])
        self._empty_state_message = "Waiting for model output…"
        self.status_label.setText(self._empty_state_message)
        self._sync_capabilities()

    def set_variable_definition(self, definition: dict[str, object]) -> None:
        """Set catalog metadata used together with runtime shape capability checks."""
        self._variable_definition = dict(definition)
        self._sync_capabilities()

    def set_view_mode(self, mode: str) -> None:
        if mode not in {"map", "globe"}:
            raise ValueError(f"unknown global field view: {mode}")
        self._view_mode = mode
        self.global_stack.setCurrentIndex(0 if mode == "map" else 1)
        self._render()

    def set_primary_mode(self, mode: str) -> None:
        if mode not in {"global", "region", "series"}:
            raise ValueError(f"unknown visualization mode: {mode}")
        self._primary_mode = mode
        self.stack.setCurrentIndex({"global": 0, "region": 1, "series": 2}[mode])
        self.region_controls.setVisible(mode == "region")
        self.color_controls.setVisible(mode != "series")
        self._render()

    def set_color_mode(self, mode: str) -> None:
        if mode not in {"fixed", "dynamic"}:
            raise ValueError(f"unknown colorbar mode: {mode}")
        self._color_mode = mode
        if self.color_mode_control.value() != mode:
            self.color_mode_control.set_value(mode)
        fixed = mode == "fixed"
        self.vmin_edit.setEnabled(fixed)
        self.vmax_edit.setEnabled(fixed)
        self.apply_levels_button.setEnabled(fixed)
        self.auto_levels_button.setEnabled(fixed)
        self._render()

    def set_colormap(self, name: str) -> None:
        pg.colormap.get(name)
        self._colormap = name
        color_map = pg.colormap.get(name)
        self.image.setColorMap(color_map)
        self.region_image.setColorMap(color_map)
        self.globe.set_colormap(name)
        self.slice_view.set_colormap(name)
        self.color_scale.set_colormap(name)
        self._render()

    def _set_region_view(self, mode: str) -> None:
        dimensions = tuple(self._variable_definition.get("dimensions", ()))
        if mode == "3d" and (
            self._frame is None or self._frame.layers < 2 or "level" not in dimensions
        ):
            self.status_label.setText("3-D structure unavailable: this variable has no vertical coordinate.")
            self.region_view_control.set_value("2d")
            return
        self._region_view = mode
        if self.region_view_control.value() != mode:
            self.region_view_control.set_value(mode)
        self.region_stack.setCurrentIndex(0 if mode == "2d" else 1)
        self._sync_region_control_visibility()
        self._render()

    def _sync_capabilities(self) -> None:
        runtime_layers = self._frames[-1].layers if self._frames else 1
        dimensions = tuple(self._variable_definition.get("dimensions", ()))
        layers = runtime_layers if "level" in dimensions else 1
        self.region_view_control.set_option_enabled("3d", layers > 1)
        self.primary_control.set_option_enabled("series", len(self._frames) > 1)
        if layers > 1:
            self.slice_count.blockSignals(True)
            self.slice_count.setRange(2, layers)
            self.slice_count.setValue(min(self.slice_count.value(), layers))
            self.slice_count.blockSignals(False)
        maximum = max(0, layers - 1)
        for spin in (self.layer_from, self.layer_to):
            spin.blockSignals(True)
            spin.setRange(0, maximum)
        self.layer_from.setValue(min(self.layer_from.value(), maximum))
        self.layer_to.setValue(max(self.layer_from.value(), maximum))
        self.layer_from.blockSignals(False)
        self.layer_to.blockSignals(False)
        if layers < 2 and self._region_view == "3d":
            self._region_view = "2d"
            self.region_view_control.set_value("2d")
        self._sync_region_control_visibility()

    def _sync_region_control_visibility(self) -> None:
        vertical = self.region_view_control.option_enabled("3d")
        visible = vertical and self._region_view == "3d"
        for widget in (self.layers_label, self.layer_from, self.layer_to, self.slice_count):
            widget.setVisible(visible)

    def show_frame(self, frame: Frame, level: int, auto_levels: bool | None = None) -> None:
        if not self._frames or frame.name != self._frames[-1].name:
            self.set_frames((frame,))
        self._frame = frame
        self._level = min(max(int(level), 0), frame.layers - 1)
        self.current_name = frame.name
        self.current_level = self._level
        self._empty_state_message = "Waiting for model output…"
        self._render()

    def clear_history(self, name: str, level: int) -> None:
        self.current_name, self.current_level = name, level
        self.steps.clear()
        self.means.clear()
        self.series_curve.setData([], [])

    def _bounds(self) -> RegionBounds:
        return RegionBounds(**{name: spin.value() for name, spin in self.bound_spins.items()})

    def _vertical_range_changed(self) -> None:
        if self.layer_from.value() > self.layer_to.value():
            if self.sender() is self.layer_from:
                self.layer_to.setValue(self.layer_from.value())
            else:
                self.layer_from.setValue(self.layer_to.value())
        self._fixed_cache.clear()
        self._render()

    def _volume_render_mode_changed(self, mode: str) -> None:
        if self._frame is None or self._primary_mode != "region" or self._region_view != "3d":
            return
        selected_layers = self.layer_to.value() - self.layer_from.value() + 1
        self.status_label.setText(
            f"3-D {mode} · {self.slice_view.layer_count} of {selected_layers} sigma model "
            f"levels sampled · {frame_time_label(self._frame)}"
        )

    def _scope_values(self, frame: Frame) -> np.ndarray:
        if self._primary_mode == "region":
            subset = region_subset(frame.data, self._bounds()).values
            if self._region_view == "3d" and frame.layers > 1:
                return subset[self.layer_from.value() : self.layer_to.value() + 1]
            return subset[min(self._level, frame.layers - 1)]
        return frame.layer(min(self._level, frame.layers - 1))

    def _color_levels(self) -> tuple[float, float]:
        assert self._frame is not None
        nonnegative = self._frame.name in SEA_ICE_VARIABLES
        if self._color_mode == "fixed" and self._manual_levels is not None:
            self.color_note.setText("Manual range · retained across time steps")
            return self._manual_levels
        if self._color_mode == "dynamic":
            self.color_note.setText("Current time step and active spatial/vertical scope")
            return finite_color_range(self._scope_values(self._frame), nonnegative=nonnegative)
        bounds = self._bounds() if self._primary_mode == "region" else None
        signature = tuple((frame.step, frame.year, frame.month, frame.day, frame.hour, frame.minute, frame.layers) for frame in self._frames)
        key = (self._frame.name, signature, self._primary_mode, self._region_view, self._level, bounds, self.layer_from.value(), self.layer_to.value())
        cached = self._fixed_cache.get(key)
        if cached is None:
            self.status_label.setText("Computing fixed color range from loaded time steps…")
            self.status_label.repaint()
            cached = finite_color_range(
                (self._scope_values(frame) for frame in self._frames),
                nonnegative=nonnegative,
            )
            if len(self._fixed_cache) >= 16:
                self._fixed_cache.pop(next(iter(self._fixed_cache)))
            self._fixed_cache[key] = cached
        self.color_note.setText(f"Fixed automatic range · all {len(self._frames)} loaded time step(s)")
        return cached

    def _apply_manual_levels(self) -> None:
        try:
            lower, upper = float(self.vmin_edit.text()), float(self.vmax_edit.text())
        except ValueError:
            self.status_label.setText("Invalid colorbar: enter numeric vmin and vmax.")
            return
        if not math.isfinite(lower) or not math.isfinite(upper) or lower >= upper:
            self.status_label.setText("Invalid colorbar: vmin and vmax must be finite and vmin < vmax.")
            return
        self._manual_levels = (lower, upper)
        self.set_color_mode("fixed")

    def _clear_manual_levels(self) -> None:
        self._manual_levels = None
        self.vmin_edit.clear()
        self.vmax_edit.clear()
        self._fixed_cache.clear()
        self._render()

    def _render(self) -> None:
        if self._frame is None:
            self._clear_field_display(self._empty_state_message)
            self.status_label.setText(self._empty_state_message)
            return
        try:
            if self._primary_mode == "series":
                self._render_series()
            else:
                levels = self._color_levels()
                if self._primary_mode == "global":
                    self._render_global(levels)
                else:
                    self._render_region(levels)
            self._show_sea_ice_status()
        except ValueError as error:
            self.status_label.setText(str(error))

    def show_unavailable(self, message: str) -> None:
        """Drop stale field imagery when the selected variable has no frames.

        Called when the current source (live run or loaded stream) does not
        contain the selected variable: keeping the previous map, title and
        colorbar on screen would silently mislabel another variable's data.
        """
        self._frames = ()
        self._frame = None
        self._level = 0
        self._sync_capabilities()
        self._empty_state_message = message
        self._clear_field_display(message)
        self.status_label.setText(message)

    def _clear_field_display(self, message: str = "") -> None:
        self.image.clear()
        self.region_image.clear()
        self.map_plot.setTitle("")
        self.region_plot.setTitle("")
        self.globe.clear_frame(message or None)
        self.slice_view.clear()
        self.series_curve.setData([], [])
        self.coverage_label.clear()
        self.color_scale.set_range(0.0, 1.0, "")

    def _show_sea_ice_status(self) -> None:
        assert self._frame is not None
        if self._frame.name not in SEA_ICE_VARIABLES:
            return
        field = np.ma.masked_invalid(self._scope_values(self._frame)).compressed()
        if not field.size:
            note = "No valid sea-ice data in the current field."
        elif np.all(field == 0.0):
            note = "No sea ice in the current field (all valid values are 0)."
        else:
            return
        self.status_label.setText(f"{self.status_label.text()} · {note}")

    def _render_global(self, levels: tuple[float, float]) -> None:
        assert self._frame is not None
        field = np.asanyarray(self._frame.layer(self._level))
        if np.ma.isMaskedArray(field):
            field = field.filled(np.nan)
        latitudes, longitudes = gaussian_grid(*field.shape)
        order = np.argsort(longitudes)
        raster = resample_latitude_for_image(np.asarray(field)[:, order], latitudes)
        self.image.setImage(raster, autoLevels=False, levels=levels)
        self.image.setRect(QRectF(-180.0, -90.0, 360.0, 180.0))
        self.map_plot.setXRange(-180, 180, padding=0)
        self.map_plot.setYRange(-90, 90, padding=0)
        self.map_plot.setTitle(f"{self._frame.name} · {self._frame.unit} · {frame_time_label(self._frame)}", color=TEXT, size="11pt")
        self.globe.set_frame(self._frame, self._level, False, levels=levels)
        self.color_scale.set_range(*levels, self._frame.unit)
        finite = int(np.count_nonzero(np.isfinite(field)))
        total = field.size
        suffix = "" if finite == total else f" · finite coverage {100.0 * finite / total:.1f}%"
        self.status_label.setText(f"Global Gaussian grid {field.shape[1]}×{field.shape[0]} · level {self._level}{suffix}")

    def _render_region(self, levels: tuple[float, float]) -> None:
        assert self._frame is not None
        bounds = self._bounds()
        subset = region_subset(self._frame.data, bounds)
        if self._region_view == "3d":
            start, stop = self.layer_from.value(), self.layer_to.value() + 1
            values = np.asarray(subset.values[start:stop])
            if values.shape[0] < 2:
                raise ValueError("3-D structure needs at least two selected vertical layers")
            count = min(self.slice_count.value(), values.shape[0])
            indices = np.unique(np.rint(np.linspace(0, values.shape[0] - 1, count)).astype(int))
            shown = values[indices]
            self.slice_view.set_data(
                shown,
                levels=levels,
                name=self._frame.name,
                unit=self._frame.unit,
                bounds=(bounds.west, bounds.east, bounds.south, bounds.north),
                layer_indices=indices + start,
                total_layers=self._frame.layers,
            )
            self.status_label.setText(
                f"3-D {self.slice_view.render_mode} · {shown.shape[0]} of {values.shape[0]} sigma model levels sampled · {frame_time_label(self._frame)}"
            )
        else:
            field = np.asarray(subset.values[self._level])
            raster = resample_latitude_for_image(field, subset.latitudes)
            dx = 360.0 / self._frame.data.shape[-1]
            x0 = float(subset.display_longitudes[0] - dx / 2)
            x1 = float(subset.display_longitudes[-1] + dx / 2)
            y0 = float(subset.latitudes[-1])
            y1 = float(subset.latitudes[0])
            self.region_image.setImage(raster, autoLevels=False, levels=levels)
            self.region_image.setRect(QRectF(x0, y0, x1 - x0, y1 - y0))
            self.region_plot.setXRange(x0, x1, padding=0)
            self.region_plot.setYRange(bounds.south, bounds.north, padding=0)
            self._update_region_coastlines(bounds)
            self.region_plot.setTitle(f"{self._frame.name} · {self._frame.unit} · level {self._level} · {frame_time_label(self._frame)}", color=TEXT, size="11pt")
            crossing = " · crosses date line" if subset.crosses_dateline else ""
            self.status_label.setText(f"Regional 2-D field · {field.shape[1]}×{field.shape[0]} cells{crossing}")
        self.color_scale.set_range(*levels, self._frame.unit)

    def _update_region_coastlines(self, bounds: RegionBounds) -> None:
        x, y = flattened_coastlines()
        adjusted = np.asarray(x, dtype=np.float64).copy()
        if bounds.crosses_dateline:
            finite = np.isfinite(adjusted)
            adjusted[finite & (adjusted < bounds.west)] += 360.0
        for item in self.region_coasts:
            item.setData(adjusted, y, connect="finite")

    def _render_series(self) -> None:
        assert self._frame is not None
        if len(self._frames) < 2:
            self.series_curve.setData([], [])
            self.status_label.setText("Global mean time series unavailable: fewer than two time coordinates are loaded.")
            self.coverage_label.setText("The pstream format supplies time through frame timestamps; load or run multiple frames.")
            return
        means: list[float] = []
        coverage: list[float] = []
        labels: list[str] = []
        frames = self._frames[-self.max_points :]
        for frame in frames:
            level = min(self._level, frame.layers - 1)
            mean, valid = frame.global_mean_with_coverage(level)
            means.append(mean)
            coverage.append(valid)
            labels.append(frame_time_label(frame).split("  ·", 1)[0])
        x = np.arange(len(frames), dtype=np.float64)
        if np.isfinite(means).any():
            self.series_curve.setData(x, means)
        else:
            self.series_curve.setData([], [])
        tick_step = max(1, len(frames) // 6)
        ticks = [(float(index), labels[index]) for index in range(0, len(frames), tick_step)]
        if ticks[-1][0] != len(frames) - 1:
            ticks.append((float(len(frames) - 1), labels[-1]))
        self.series_plot.getAxis("bottom").setTicks([ticks])
        self.series_plot.setLabel("left", self._frame.name, units=self._frame.unit)
        self.series_plot.setTitle(f"Area-weighted global mean · sigma model level {self._level}", color=TEXT, size="11pt")
        self.series_plot.autoRange()
        self.steps = [frame.step for frame in frames]
        self.means = means
        minimum_coverage = min(coverage)
        if minimum_coverage >= 1.0 - 1.0e-12:
            detail = "Full global Gaussian grid coverage at every time step."
        else:
            detail = f"Finite-data coverage varies from {100 * minimum_coverage:.1f}% to {100 * max(coverage):.1f}%; each mean renormalizes valid Gaussian cell weights."
        self.coverage_label.setText("Gaussian quadrature area weights; computed from full global frames and unaffected by regional bounds. " + detail)
        self.status_label.setText(f"{len(frames)} time coordinates · selected vertical level {self._level}")

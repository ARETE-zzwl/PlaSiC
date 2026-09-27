from __future__ import annotations

import math
import os

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import QRectF, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QMatrix4x4, QPainter, QVector3D
from PySide6.QtWidgets import (
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QSlider,
    QSpinBox,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from plasic_app.ui.widgets import SegmentedControl
from plasic_app.visualization.data import finite_color_range

try:
    import pyqtgraph.opengl as gl
except (ImportError, ModuleNotFoundError) as error:  # pragma: no cover - installation dependent
    gl = None
    _OPENGL_IMPORT_ERROR = str(error)
else:
    _OPENGL_IMPORT_ERROR = ""


if gl is not None:
    class _LabeledGLViewWidget(gl.GLViewWidget):
        """GL view that paints all projected axis text in one safe overlay pass."""

        def __init__(self, parent: QWidget | None = None) -> None:
            super().__init__(parent=parent)
            self._tick_annotations: list[tuple[tuple[float, float, float], str]] = []
            self._title_annotations: list[tuple[tuple[float, float, float], str]] = []

        def set_axis_annotations(
            self,
            ticks: list[tuple[tuple[float, float, float], str]],
            titles: list[tuple[tuple[float, float, float], str]],
        ) -> None:
            self._tick_annotations = ticks
            self._title_annotations = titles
            self.update()

        def paintGL(self) -> None:
            super().paintGL()
            if not self._tick_annotations and not self._title_annotations:
                return
            viewport = self.getViewport()
            view_rect = QRectF(self.rect())
            viewport_matrix = QMatrix4x4()
            viewport_matrix.viewport(
                view_rect.left(),
                view_rect.bottom(),
                view_rect.width(),
                -view_rect.height(),
            )
            projection = (
                viewport_matrix
                * self.projectionMatrix(viewport, viewport)
                * self.viewMatrix()
            )
            painter = QPainter(self)
            painter.setRenderHints(
                QPainter.RenderHint.Antialiasing
                | QPainter.RenderHint.TextAntialiasing
            )
            self._paint_annotations(
                painter,
                projection,
                self._tick_annotations,
                QFont("Helvetica", 9),
                QColor(225, 232, 240, 230),
            )
            self._paint_annotations(
                painter,
                projection,
                self._title_annotations,
                QFont("Helvetica", 10, QFont.Weight.DemiBold),
                QColor(245, 247, 250, 245),
            )
            painter.end()

        @staticmethod
        def _paint_annotations(
            painter: QPainter,
            projection: QMatrix4x4,
            annotations: list[tuple[tuple[float, float, float], str]],
            font: QFont,
            color: QColor,
        ) -> None:
            painter.setFont(font)
            painter.setPen(color)
            metrics = painter.fontMetrics()
            for position, text in annotations:
                point = projection.map(QVector3D(*position)).toPointF()
                bounds = metrics.boundingRect(text)
                target = QRectF(
                    point.x() - bounds.width() / 2.0 - 3.0,
                    point.y() - bounds.height() / 2.0 - 2.0,
                    bounds.width() + 6.0,
                    bounds.height() + 4.0,
                )
                painter.drawText(target, Qt.AlignmentFlag.AlignCenter, text)
else:  # pragma: no cover - used only when the optional GL import is unavailable
    _LabeledGLViewWidget = None


def _float_volume(values: np.ndarray) -> np.ndarray:
    array = np.asanyarray(values)
    if np.ma.isMaskedArray(array):
        array = array.filled(np.nan)
    return np.asarray(array, dtype=np.float32)


def _normalize(values: np.ndarray, levels: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    lower, upper = levels
    if not math.isfinite(lower) or not math.isfinite(upper) or lower >= upper:
        lower, upper = finite_color_range(values)
    finite = np.isfinite(values)
    normalized = np.zeros(values.shape, dtype=np.float32)
    normalized[finite] = np.clip((values[finite] - lower) / (upper - lower), 0.0, 1.0)
    return normalized, finite


def colorize_field(
    values: np.ndarray,
    *,
    levels: tuple[float, float],
    lookup: np.ndarray,
    alpha: int = 230,
) -> np.ndarray:
    """Map a 2-D or 3-D field to RGBA while keeping invalid cells transparent."""
    array = _float_volume(values)
    normalized, finite = _normalize(array, levels)
    indices = np.rint(normalized * (len(lookup) - 1)).astype(np.int32)
    rgba = np.asarray(lookup[indices], dtype=np.uint8).copy()
    rgba[..., 3] = np.where(finite, np.clip(alpha, 0, 255), 0).astype(np.uint8)
    return np.ascontiguousarray(rgba)


def volume_rgba(
    values: np.ndarray,
    *,
    levels: tuple[float, float],
    lookup: np.ndarray,
    opacity: int,
    threshold: int,
) -> np.ndarray:
    """Create an RGBA transfer function suited to meteorological volumes."""
    array = _float_volume(values)
    normalized, finite = _normalize(array, levels)
    indices = np.rint(normalized * (len(lookup) - 1)).astype(np.int32)
    rgba = np.asarray(lookup[indices], dtype=np.uint8).copy()
    lower, upper = levels
    if lower < 0.0 < upper:
        zero_position = float(np.clip((0.0 - lower) / (upper - lower), 0.0, 1.0))
        denominator = np.where(normalized >= zero_position, 1.0 - zero_position, zero_position)
        strength = np.abs(normalized - zero_position) / np.maximum(denominator, 1.0e-6)
    else:
        strength = normalized
    cutoff = float(np.clip(threshold, 0, 95)) / 100.0
    ramp = np.clip((strength - cutoff) / max(1.0 - cutoff, 1.0e-6), 0.0, 1.0)
    alpha = 255.0 * (float(np.clip(opacity, 1, 100)) / 100.0) * np.power(ramp, 1.25)
    rgba[..., 3] = np.where(finite, np.rint(alpha), 0).astype(np.uint8)
    return np.ascontiguousarray(rgba)


def display_volume(values: np.ndarray, maximum_xy: int = 112, maximum_z: int = 48) -> np.ndarray:
    """Return a bounded, display-only volume with vertical interpolation."""
    array = _float_volume(values)
    y_stride = max(1, math.ceil(array.shape[1] / maximum_xy))
    x_stride = max(1, math.ceil(array.shape[2] / maximum_xy))
    array = array[:, ::y_stride, ::x_stride]
    layers = array.shape[0]
    target_layers = min(maximum_z, max(layers, layers * 4))
    if target_layers == layers:
        return np.ascontiguousarray(array)
    position = np.linspace(0.0, layers - 1.0, target_layers, dtype=np.float32)
    lower = np.floor(position).astype(np.int32)
    upper = np.minimum(lower + 1, layers - 1)
    weight = (position - lower)[:, None, None]
    low_values, high_values = array[lower], array[upper]
    result = low_values * (1.0 - weight) + high_values * weight
    result[~np.isfinite(low_values) | ~np.isfinite(high_values)] = np.nan
    return np.ascontiguousarray(result, dtype=np.float32)


def vertical_scene_fraction(index: int, count: int) -> float:
    """Map ascending sigma layer indices from the lower to the upper box face."""
    if count < 2:
        return 0.5
    return float(np.clip(index, 0, count - 1) / (count - 1))


def latitude_scene_fraction(index: int, count: int) -> float:
    """Map north-to-south source rows onto the far-to-near scene direction."""
    if count < 2:
        return 0.5
    return 1.0 - float(np.clip(index, 0, count - 1) / (count - 1))


class _ThreeSliceFallback(QWidget):
    """Software fallback that still exposes the three orthogonal sections."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QGridLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)
        self.images: list[pg.ImageItem] = []
        self.plots: list[pg.PlotWidget] = []
        labels = (
            ("Horizontal · longitude × latitude", "Longitude", "Latitude"),
            ("Zonal section · longitude × vertical", "Longitude", "Vertical sample"),
            ("Meridional section · latitude × vertical", "Latitude", "Vertical sample"),
        )
        for column, (title, bottom, left) in enumerate(labels):
            plot = pg.PlotWidget()
            plot.setBackground("#000000")
            plot.setTitle(title, color="#E8EDF3", size="9pt")
            plot.setLabel("bottom", bottom)
            plot.setLabel("left", left)
            plot.hideButtons()
            plot.setMenuEnabled(False)
            image = pg.ImageItem(axisOrder="row-major")
            plot.addItem(image)
            layout.addWidget(plot, 0, column)
            self.plots.append(plot)
            self.images.append(image)

    def update_slices(
        self,
        values: np.ndarray,
        levels: tuple[float, float],
        positions: tuple[int, int, int],
    ) -> None:
        level, latitude, longitude = positions
        sections = (
            values[level, ::-1, :],
            values[:, latitude, :],
            values[:, ::-1, longitude],
        )
        for index, (image, section) in enumerate(zip(self.images, sections, strict=True)):
            image.setImage(section, autoLevels=False, levels=levels)
            image.setRect(QRectF(0, 0, section.shape[1], section.shape[0]))
            self.plots[index].getViewBox().invertY(False)
            self.plots[index].autoRange()


class LayerSliceView(QWidget):
    """Interactive regional 3-D volume with orthogonal meteorological slices."""

    renderModeChanged = Signal(str)

    _BOX_EDGES = (
        (0, 1), (1, 3), (3, 2), (2, 0),
        (4, 5), (5, 7), (7, 6), (6, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    )

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumSize(560, 380)
        self._values: np.ndarray | None = None
        self._levels = (0.0, 1.0)
        self._unit = ""
        self._name = ""
        self._west, self._east = -180.0, 180.0
        self._south, self._north = -90.0, 90.0
        self._layer_indices = np.empty(0, dtype=np.int32)
        self._total_layers = 0
        self._colormap = "CET-D1"
        self._lookup = self._make_lookup(self._colormap)
        self._display_voxel_count = 0
        self._volume_item = None
        self._slice_items: list[object] = []
        self._box_item = None
        self._transfer_timer = QTimer(self)
        self._transfer_timer.setSingleShot(True)
        self._transfer_timer.setInterval(60)
        self._transfer_timer.timeout.connect(self._refresh_scene)

        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(6)

        controls = QHBoxLayout()
        controls.setSpacing(6)
        controls.addWidget(QLabel("Render"))
        self.render_control = SegmentedControl(
            [("Volume", "volume"), ("Slices", "slices"), ("Hybrid", "hybrid")]
        )
        controls.addWidget(self.render_control)
        self.opacity_slider = self._make_slider(5, 100, 42)
        self.threshold_slider = self._make_slider(0, 95, 30)
        controls.addWidget(QLabel("Opacity"))
        controls.addWidget(self.opacity_slider)
        self.opacity_value = QLabel("42%")
        controls.addWidget(self.opacity_value)
        controls.addWidget(QLabel("Threshold"))
        controls.addWidget(self.threshold_slider)
        self.threshold_value = QLabel("30%")
        controls.addWidget(self.threshold_value)
        controls.addStretch(1)
        root.addLayout(controls)

        slice_controls = QHBoxLayout()
        slice_controls.setSpacing(6)
        slice_controls.addWidget(QLabel("Slice position"))
        self.level_slice = QSpinBox()
        self.latitude_slice = QSpinBox()
        self.longitude_slice = QSpinBox()
        self.level_slice.setPrefix("vertical #")
        self.latitude_slice.setPrefix("latitude #")
        self.longitude_slice.setPrefix("longitude #")
        for spin in (self.level_slice, self.latitude_slice, self.longitude_slice):
            spin.setMinimumWidth(112)
            slice_controls.addWidget(spin)
        self.position_label = QLabel()
        self.position_label.setProperty("role", "secondary")
        slice_controls.addWidget(self.position_label, 1)
        root.addLayout(slice_controls)

        stage = QWidget()
        stage.setStyleSheet("background: #000000; color: #E8EDF3;")
        stage_layout = QVBoxLayout(stage)
        stage_layout.setContentsMargins(10, 8, 10, 6)
        stage_layout.setSpacing(3)
        self.title_label = QLabel("No vertical structure available")
        self.title_label.setStyleSheet("font-weight: 600; color: #F5F7FA;")
        self.detail_label = QLabel()
        self.detail_label.setStyleSheet("color: #B4BECA;")
        self.detail_label.setWordWrap(True)
        stage_layout.addWidget(self.title_label)
        stage_layout.addWidget(self.detail_label)

        self.scene_stack = QStackedWidget()
        self.fallback = _ThreeSliceFallback()
        self._opengl_error = _OPENGL_IMPORT_ERROR
        self.gl_view = None
        platform_name = os.environ.get("QT_QPA_PLATFORM", "").lower()
        if (
            gl is not None
            and os.environ.get("PLASIC_FORCE_SOFTWARE_3D") != "1"
            and platform_name not in {"offscreen", "minimal"}
        ):
            try:
                assert _LabeledGLViewWidget is not None
                self.gl_view = _LabeledGLViewWidget()
                self.gl_view.setBackgroundColor(QColor("#000000"))
                self.gl_view.setCameraPosition(distance=42.0, elevation=22.0, azimuth=-42.0)
                self.scene_stack.addWidget(self.gl_view)
            except Exception as error:  # pragma: no cover - driver/platform dependent
                self._opengl_error = str(error)
                self.gl_view = None
        self.scene_stack.addWidget(self.fallback)
        self.scene_stack.setCurrentWidget(self.gl_view or self.fallback)
        stage_layout.addWidget(self.scene_stack, 1)
        self.help_label = QLabel(
            "Left-drag to rotate · Wheel to zoom · Middle-drag or Ctrl+left-drag to pan"
            if self.gl_view is not None
            else "OpenGL unavailable · showing linked orthogonal 2-D sections"
        )
        self.help_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.help_label.setStyleSheet("color: #8F9BAA;")
        stage_layout.addWidget(self.help_label)
        root.addWidget(stage, 1)

        if self.gl_view is None:
            self.render_control.set_option_enabled("volume", False)
            self.render_control.set_option_enabled("hybrid", False)
            self.render_control.set_value("slices")
            if self._opengl_error:
                self.help_label.setToolTip(self._opengl_error)

        self.render_control.valueChanged.connect(self._render_mode_changed)
        self.opacity_slider.valueChanged.connect(self._transfer_changed)
        self.threshold_slider.valueChanged.connect(self._transfer_changed)
        for spin in (self.level_slice, self.latitude_slice, self.longitude_slice):
            spin.valueChanged.connect(self._slice_changed)
        self._sync_control_visibility()

    @staticmethod
    def _make_slider(minimum: int, maximum: int, value: int) -> QSlider:
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setRange(minimum, maximum)
        slider.setValue(value)
        slider.setMaximumWidth(120)
        return slider

    @staticmethod
    def _make_lookup(name: str) -> np.ndarray:
        return np.asarray(
            pg.colormap.get(name).getLookupTable(nPts=256, alpha=True), dtype=np.uint8
        )

    @property
    def layer_count(self) -> int:
        return 0 if self._values is None else int(self._values.shape[0])

    @property
    def using_opengl(self) -> bool:
        return self.gl_view is not None

    @property
    def render_mode(self) -> str:
        return self.render_control.value()

    def set_colormap(self, name: str) -> None:
        self._colormap = name
        self._lookup = self._make_lookup(name)
        self._refresh_scene()

    def set_data(
        self,
        values: np.ndarray,
        *,
        levels: tuple[float, float] | None,
        name: str,
        unit: str,
        bounds: tuple[float, float, float, float],
        layer_indices: np.ndarray | None = None,
        total_layers: int | None = None,
    ) -> None:
        array = _float_volume(values)
        if array.ndim != 3 or array.shape[0] < 2:
            raise ValueError("3-D structure needs at least two vertical layers")
        previous_shape = None if self._values is None else self._values.shape
        previous_positions = (
            self.level_slice.value(),
            self.latitude_slice.value(),
            self.longitude_slice.value(),
        )
        self._values = array
        self._levels = levels or finite_color_range(array)
        self._name, self._unit = name, unit
        self._west, self._east, self._south, self._north = bounds
        if layer_indices is None:
            self._layer_indices = np.arange(array.shape[0], dtype=np.int32)
        else:
            indices = np.asarray(layer_indices, dtype=np.int32)
            if indices.shape != (array.shape[0],):
                raise ValueError("layer indices must match the displayed vertical layers")
            self._layer_indices = indices
        self._total_layers = int(total_layers or array.shape[0])
        centers = (array.shape[0] // 2, array.shape[1] // 2, array.shape[2] // 2)
        requested_positions = previous_positions if previous_shape == array.shape else centers
        for spin, maximum, value in zip(
            (self.level_slice, self.latitude_slice, self.longitude_slice),
            (array.shape[0] - 1, array.shape[1] - 1, array.shape[2] - 1),
            requested_positions,
            strict=True,
        ):
            spin.blockSignals(True)
            spin.setRange(0, maximum)
            spin.setValue(value)
            spin.blockSignals(False)
        self._refresh_scene()

    def clear(self) -> None:
        """Reset the 3-D stage when the selected variable has no frames."""
        self._values = None
        self._levels = (0.0, 1.0)
        self._name = ""
        self._unit = ""
        self._layer_indices = np.empty(0, dtype=np.int32)
        self._total_layers = 0
        self._display_voxel_count = 0
        self._volume_item = None
        self._slice_items = []
        self.position_label.clear()
        self.title_label.setText("No vertical structure available")
        self.detail_label.clear()
        for image in self.fallback.images:
            image.clear()
        if self.gl_view is not None:
            self.gl_view.clear()

    def _horizontal_aspect(self) -> float:
        longitude_span = (self._east - self._west) % 360.0
        if longitude_span == 0.0:
            longitude_span = 360.0
        latitude_span = max(self._north - self._south, 1.0)
        mid_latitude = math.radians((self._north + self._south) * 0.5)
        ratio = longitude_span * max(math.cos(mid_latitude), 0.2) / latitude_span
        return float(np.clip(ratio, 0.75, 2.2))

    def _scene_dimensions(self) -> tuple[float, float, float]:
        return 11.0 * self._horizontal_aspect(), 11.0, 6.5

    def _render_mode_changed(self, mode: str) -> None:
        self._sync_control_visibility()
        self._refresh_scene()
        self.renderModeChanged.emit(mode)

    def _transfer_changed(self, _value: int) -> None:
        self.opacity_value.setText(f"{self.opacity_slider.value()}%")
        self.threshold_value.setText(f"{self.threshold_slider.value()}%")
        self._transfer_timer.start()

    def _slice_changed(self, _value: int) -> None:
        self._refresh_scene()

    def _sync_control_visibility(self) -> None:
        mode = self.render_mode
        volume_visible = mode in {"volume", "hybrid"}
        slice_visible = mode in {"slices", "hybrid"}
        for widget in (
            self.opacity_slider,
            self.opacity_value,
            self.threshold_slider,
            self.threshold_value,
        ):
            widget.setEnabled(volume_visible)
        for widget in (self.level_slice, self.latitude_slice, self.longitude_slice):
            widget.setEnabled(slice_visible)

    def _refresh_scene(self) -> None:
        if self._values is None:
            return
        level = self.level_slice.value()
        latitude = self.latitude_slice.value()
        longitude = self.longitude_slice.value()
        model_level = int(self._layer_indices[level])
        latitude_value = self._north - (self._north - self._south) * latitude / max(self._values.shape[1] - 1, 1)
        longitude_span = (self._east - self._west) % 360.0 or 360.0
        longitude_value = self._west + longitude_span * longitude / max(self._values.shape[2] - 1, 1)
        if longitude_value > 180.0:
            longitude_value -= 360.0
        self.position_label.setText(
            f"L{model_level} · {latitude_value:.1f}° lat · {longitude_value:.1f}° lon"
        )
        self.title_label.setText(f"{self._name} · {self._unit} · 3-D {self.render_mode}")
        start, end = int(self._layer_indices[0]), int(self._layer_indices[-1])
        backend = "OpenGL volume texture" if self.gl_view is not None else "software slice fallback"
        self.detail_label.setText(
            f"Longitude {self._west:g}° → {self._east:g}° · Latitude {self._south:g}° → {self._north:g}° · "
            f"Vertical: sigma model-level index [1], L{start} at lower face → "
            f"L{end} at upper face · numerical sigma values are not stored in pstream · {backend}"
        )
        positions = (level, latitude, longitude)
        self.fallback.update_slices(self._values, self._levels, positions)
        if self.gl_view is not None:
            self._update_gl_scene(positions)

    def _update_gl_scene(self, positions: tuple[int, int, int]) -> None:
        assert gl is not None and self.gl_view is not None and self._values is not None
        self.gl_view.clear()
        self._volume_item = None
        self._slice_items = []
        mode = self.render_mode
        box_x, box_y, box_z = self._scene_dimensions()
        if mode in {"volume", "hybrid"}:
            volume = display_volume(self._values)
            rgba = volume_rgba(
                volume[:, ::-1, :].transpose(2, 1, 0),
                levels=self._levels,
                lookup=self._lookup,
                opacity=self.opacity_slider.value(),
                threshold=self.threshold_slider.value(),
            )
            self._display_voxel_count = int(np.count_nonzero(rgba[..., 3]))
            # Additive blending prevents transparent volume samples from writing a
            # foreground depth wall that would hide the later slice planes.
            item = gl.GLVolumeItem(rgba, sliceDensity=1, smooth=True, glOptions="additive")
            nx, ny, nz = rgba.shape[:3]
            item.scale(box_x / nx, box_y / ny, box_z / nz)
            item.translate(-box_x / 2.0, -box_y / 2.0, -box_z / 2.0)
            self.gl_view.addItem(item)
            self._volume_item = item
        else:
            self._display_voxel_count = 0
        if mode in {"slices", "hybrid"}:
            self._add_slice_planes(positions, box_x, box_y, box_z)
        self._add_box(box_x, box_y, box_z)

    def _add_slice_planes(
        self,
        positions: tuple[int, int, int],
        box_x: float,
        box_y: float,
        box_z: float,
    ) -> None:
        assert gl is not None and self.gl_view is not None and self._values is not None
        level, latitude, longitude = positions
        nlev, nlat, nlon = self._values.shape
        slice_alpha = 150 if self.render_mode == "hybrid" else 245

        horizontal = colorize_field(
            self._values[level, ::-1, :],
            levels=self._levels,
            lookup=self._lookup,
            alpha=slice_alpha,
        ).transpose(1, 0, 2)
        horizontal_item = gl.GLImageItem(horizontal, smooth=True, glOptions="translucent")
        horizontal_item.scale(box_x / nlon, box_y / nlat, 1.0)
        horizontal_item.translate(
            -box_x / 2.0,
            -box_y / 2.0,
            -box_z / 2.0 + box_z * vertical_scene_fraction(level, nlev),
        )

        zonal = colorize_field(
            self._values[:, latitude, :],
            levels=self._levels,
            lookup=self._lookup,
            alpha=slice_alpha,
        ).transpose(1, 0, 2)
        zonal_item = gl.GLImageItem(zonal, smooth=True, glOptions="translucent")
        zonal_item.scale(box_x / nlon, box_z / nlev, 1.0)
        zonal_item.rotate(90.0, 1.0, 0.0, 0.0)
        zonal_item.translate(
            -box_x / 2.0,
            -box_y / 2.0 + box_y * latitude_scene_fraction(latitude, nlat),
            -box_z / 2.0,
        )

        meridional = colorize_field(
            self._values[:, ::-1, longitude],
            levels=self._levels,
            lookup=self._lookup,
            alpha=slice_alpha,
        ).transpose(1, 0, 2)
        meridional_item = gl.GLImageItem(meridional, smooth=True, glOptions="translucent")
        meridional_item.scale(box_y / nlat, box_z / nlev, 1.0)
        meridional_item.rotate(90.0, 0.0, 0.0, 1.0)
        meridional_item.rotate(90.0, 0.0, 1.0, 0.0)
        meridional_item.translate(
            -box_x / 2.0 + box_x * longitude / max(nlon - 1, 1),
            -box_y / 2.0,
            -box_z / 2.0,
        )
        for item in (horizontal_item, zonal_item, meridional_item):
            item.setDepthValue(10)
            self.gl_view.addItem(item)
            self._slice_items.append(item)

    def _add_box(self, box_x: float, box_y: float, box_z: float) -> None:
        assert gl is not None and self.gl_view is not None
        corners = np.asarray(
            [
                (-box_x / 2, -box_y / 2, -box_z / 2),
                (box_x / 2, -box_y / 2, -box_z / 2),
                (-box_x / 2, box_y / 2, -box_z / 2),
                (box_x / 2, box_y / 2, -box_z / 2),
                (-box_x / 2, -box_y / 2, box_z / 2),
                (box_x / 2, -box_y / 2, box_z / 2),
                (-box_x / 2, box_y / 2, box_z / 2),
                (box_x / 2, box_y / 2, box_z / 2),
            ],
            dtype=np.float32,
        )
        positions = np.asarray(
            [corners[index] for edge in self._BOX_EDGES for index in edge], dtype=np.float32
        )
        item = gl.GLLinePlotItem(
            pos=positions,
            color=(0.78, 0.84, 0.90, 0.72),
            width=1.1,
            mode="lines",
            antialias=True,
            glOptions="translucent",
        )
        item.setDepthValue(20)
        self.gl_view.addItem(item)
        self._box_item = item
        self._add_axis_ticks(box_x, box_y, box_z)

    @staticmethod
    def _tick_indices(count: int, maximum: int = 6) -> np.ndarray:
        samples = min(maximum, count)
        return np.unique(np.rint(np.linspace(0, count - 1, samples)).astype(np.int32))

    @staticmethod
    def _longitude_text(value: float) -> str:
        wrapped = ((value + 180.0) % 360.0) - 180.0
        if abs(wrapped) < 1.0e-8:
            return "0°"
        if abs(abs(wrapped) - 180.0) < 1.0e-8:
            return "180°"
        return f"{abs(wrapped):g}°{'E' if wrapped > 0 else 'W'}"

    @staticmethod
    def _latitude_text(value: float) -> str:
        if abs(value) < 1.0e-8:
            return "0°"
        return f"{abs(value):g}°{'N' if value > 0 else 'S'}"

    def _add_axis_ticks(self, box_x: float, box_y: float, box_z: float) -> None:
        """Add box-attached sigma, longitude, and latitude ticks and labels."""
        assert gl is not None and self.gl_view is not None and self._values is not None
        nlev, nlat, nlon = self._values.shape
        tick = min(box_x, box_y) * 0.035
        text_gap = tick * 1.8
        line_points: list[tuple[float, float, float]] = []
        labels: list[tuple[tuple[float, float, float], str]] = []

        # Sigma layer index: the smallest displayed index is explicitly on the
        # lower face and the largest on the upper face. Use the right-back edge
        # so its labels do not collide with the horizontal axes at the origin.
        for index in self._tick_indices(nlev):
            z = -box_z / 2.0 + box_z * vertical_scene_fraction(int(index), nlev)
            line_points.extend(
                [(box_x / 2.0, box_y / 2.0, z), (box_x / 2.0 + tick, box_y / 2.0, z)]
            )
            labels.append(
                ((box_x / 2.0 + text_gap, box_y / 2.0, z), f"σ L{int(self._layer_indices[index])}")
            )

        longitude_span = (self._east - self._west) % 360.0 or 360.0
        for index in self._tick_indices(nlon, maximum=4):
            fraction = index / max(nlon - 1, 1)
            x = -box_x / 2.0 + box_x * fraction
            line_points.extend(
                [(x, -box_y / 2.0, -box_z / 2.0), (x, -box_y / 2.0 - tick, -box_z / 2.0)]
            )
            labels.append(
                (
                    (x, -box_y / 2.0 - text_gap, -box_z / 2.0),
                    self._longitude_text(self._west + longitude_span * fraction),
                )
            )

        # Put latitude on the upper-left edge, with south on the near side and
        # north on the far side, separate from longitude on the lower-front edge.
        for index in self._tick_indices(nlat, maximum=4):
            fraction = index / max(nlat - 1, 1)
            y = -box_y / 2.0 + box_y * fraction
            latitude = self._south + (self._north - self._south) * fraction
            line_points.extend(
                [(-box_x / 2.0, y, box_z / 2.0), (-box_x / 2.0 - tick, y, box_z / 2.0)]
            )
            labels.append(
                ((-box_x / 2.0 - text_gap, y, box_z / 2.0), self._latitude_text(latitude))
            )

        tick_item = gl.GLLinePlotItem(
            pos=np.asarray(line_points, dtype=np.float32),
            color=(0.82, 0.87, 0.93, 0.9),
            width=1.2,
            mode="lines",
            antialias=True,
            glOptions="translucent",
        )
        tick_item.setDepthValue(21)
        self.gl_view.addItem(tick_item)

        axis_titles = [
            ((0.0, -box_y / 2.0 - text_gap * 2.2, -box_z / 2.0), "Longitude"),
            ((-box_x / 2.0 - text_gap * 2.1, 0.0, box_z / 2.0), "Latitude"),
            (
                (box_x / 2.0 + text_gap * 5.0, box_y / 2.0, box_z / 2.0 + tick * 2.0),
                "Sigma layer ↑",
            ),
        ]
        self.gl_view.set_axis_annotations(labels, axis_titles)

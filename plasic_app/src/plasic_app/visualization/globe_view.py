# 外部业务代码
#     │
#     ├── set_frame(frame, level, auto_levels)
#     │       ├── frame.layer(level)
#     │       └── _invalidate()
#     │                │
#     │                └── update()
#     │                     ↓ Qt事件循环
#     │
#     ├── set_colormap()
#     │       ├── _make_lookup_table()
#     │       └── _invalidate()
#     │
#     └── reset_orientation()
#             └── _invalidate()

#                     ↓

#                paintEvent()
#                     │
#         ┌───────────┴────────────┐
#         ↓                        ↓
#  _render_globe()          _coastline_path()
#         │                        │
#         ├── frame.layer()        └── coastline_segments()
#         ├── 坐标反投影
#         ├── _sample_field()
#         │      └── 双线性插值
#         ├── colormap映射
#         └── QImage
#                     ↓
#              QPainter绘制


from __future__ import annotations

import math

import numpy as np
from PySide6.QtCore import QPoint, QRectF, Qt
from PySide6.QtGui import (
    QColor,
    QImage,
    QMouseEvent,
    QPainter,
    QPainterPath,
    QPaintEvent,
    QPen,
    QRadialGradient,
    QWheelEvent,
)

from PySide6.QtWidgets import QToolTip, QWidget

# 使用 pyqtgraph 的 colormap API 生成 256 级颜色查找表（lookup table）。
import pyqtgraph as pg

# Frame 封装 PlaSiC 单个输出时刻/变量的数据，并提供 layer(level) 访问二维场。
from plasic_app.protocol.pstream import Frame
from plasic_app.ui.theme import PLOT_BACKGROUND, TEXT, TEXT_MUTED
from plasic_app.visualization.coastlines import coastline_segments


class ClimateGlobeWidget(QWidget):
    """Software-rendered orthographic globe for a PlaSiC latitude/longitude field."""

    # 初始化控件状态。
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

        # 保证控件不会缩得过小，否则球体及标题/提示文字难以辨认。
        self.setMinimumSize(360, 320)

        # 即使没有按住鼠标，也让 QWidget 持续接收 mouseMoveEvent。
        # 当前实现主要依赖拖拽时的 move，但开启 tracking 不影响逻辑。
        self.setMouseTracking(True)

        # 空闲时显示“张开手”光标，视觉上提示球体可以拖动旋转。
        self.setCursor(Qt.CursorShape.OpenHandCursor)

        # 当前显示的模型帧；None 表示还没有收到模型输出。
        self._frame: Frame | None = None

        # 无数据时在球面上方显示的提示文本。
        self._empty_message = "Waiting for model output…"

        # 当前显示的垂直层索引。
        self._level = 0

        # 是否根据当前 field 的最小/最大值自动更新色标范围。
        self._auto_levels = True

        # 实际颜色归一化范围 (lower, upper)。
        # 首次拿到 frame 之前为 None。
        self._levels: tuple[float, float] | None = None

        # pyqtgraph 颜色表名称。CET-D1 属于科学可视化色表。
        self._colormap_name = "CET-D1"

        # 预先把颜色表离散成 256×3 的 uint8 查找表，
        # 后面可通过整数索引快速把归一化场值映射为 RGB。
        self._lookup_table = self._make_lookup_table(self._colormap_name)

        # 球体朝向：yaw 控制水平旋转，pitch 控制上下倾斜。
        # 单位均为“度”，真正计算三角函数时再转成弧度。
        self._yaw = 0.0
        self._pitch = 12.0

        # 球体相对于可用绘制区域的缩放系数。
        self._zoom = 0.94

        # 保存上一次拖拽位置；None 表示当前没有进行左键拖拽。
        self._drag_position: QPoint | None = None

        self._cache_key: tuple[object, ...] | None = None
        self._cache_image: QImage | None = None

    # 该函数只依赖传入的颜色表名称，与 widget 实例状态无关，因此定义为 staticmethod。
    @staticmethod
    def _make_lookup_table(name: str) -> np.ndarray:
        # 从 pyqtgraph 获取指定 colormap，并在 [0, 1] 上均匀采样 256 个颜色。
        # alpha=False 表示这里只生成 RGB；球体透明度稍后单独构造。
        return np.asarray(
            pg.colormap.get(name).getLookupTable(
                start=0.0, stop=1.0, nPts=256, alpha=False
            ),
            dtype=np.uint8,
        )

    # 对外只读暴露当前朝向，调用方无需直接访问私有成员。
    @property
    def orientation(self) -> tuple[float, float]:
        return self._yaw, self._pitch

    def set_colormap(self, name: str) -> None:
        # 如果颜色表没有变化，不重复构建 lookup table，也不触发重绘。
        if name == self._colormap_name:
            return

        # 更新名称并重新生成 256 级 RGB 查找表。
        self._colormap_name = name
        self._lookup_table = self._make_lookup_table(name)

        # 颜色映射发生变化后，旧的缓存图像已经失效。
        self._invalidate()

    def set_frame(
        self,
        frame: Frame,
        level: int,
        auto_levels: bool,
        *,
        levels: tuple[float, float] | None = None,
    ) -> None:
        # 保存新的模型输出、层号以及自动色阶选项。
        self._frame = frame
        self._level = level
        self._auto_levels = auto_levels

        # 取出当前层对应的二维纬经度场，用于确定颜色范围。
        field = frame.layer(level)

        # 自动色阶开启时每次都重新计算；
        # 即使 auto_levels=False，第一次没有现成 _levels 时也必须初始化。
        if levels is not None:
            self._levels = levels
        elif auto_levels or self._levels is None:
            finite = np.asarray(field)[np.isfinite(field)]
            if finite.size:
                lower, upper = float(np.min(finite)), float(np.max(finite))
                if lower == upper:
                    padding = max(abs(lower) * 0.01, 0.5)
                    lower, upper = lower - padding, upper + padding
            else:
                lower, upper = 0.0, 1.0
            self._levels = (lower, upper)

        # frame/level/levels 变化后必须清除旧渲染缓存并请求 Qt 重绘。
        self._invalidate()

    def clear_frame(self, message: str | None = None) -> None:
        """Drop the current field and show the empty-state message."""
        self._frame = None
        self._levels = None
        if message:
            self._empty_message = message
        self._invalidate()

    def reset_orientation(self) -> None:
        # 恢复与初始化时完全相同的默认视角和缩放。
        self._yaw = 0.0
        self._pitch = 12.0
        self._zoom = 0.94
        self._invalidate()

    def _invalidate(self) -> None:
        # 清空像素缓存。下一个 paintEvent 会根据当前状态重新渲染。
        self._cache_key = None
        self._cache_image = None

        # update() 不会立即同步绘制，而是通知 Qt 在合适的事件循环时机触发 paintEvent。
        self.update()

    def _sample_field(
        self, field: np.ndarray, latitude: np.ndarray, longitude: np.ndarray
    ) -> np.ndarray:
        # field 约定为 [latitude_index, longitude_index] 的二维模型网格。
        nlat, nlon = field.shape

        # pstream uses Gaussian rather than equally spaced latitude rows.
        nodes, _weights = np.polynomial.legendre.leggauss(nlat)
        grid_latitude = np.arcsin(nodes)[::-1]
        row = np.interp(
            latitude,
            grid_latitude[::-1],
            np.arange(nlat, dtype=np.float64)[::-1],
        )

        # 模型经度列从 0°E 向东递增，
        # 但 arctan2 得到的是 (-pi, pi]。通过 % nlon 将负经度自然回绕到
        # [0, nlon) 的列索引空间，同时保留浮点位置供经向插值使用。
        column = (longitude / (2.0 * math.pi) * nlon) % nlon

        # 取纬向相邻的两个整数行索引。
        row0 = np.floor(row).astype(np.int32)
        # 纬度方向不周期，因此下边界必须限制在最后一行。
        row1 = np.minimum(row0 + 1, nlat - 1)

        # 取经向相邻的两个整数列索引。
        column0 = np.floor(column).astype(np.int32)
        # 经度方向是周期边界，因此最后一列的右邻居应回到第 0 列。
        column1 = (column0 + 1) % nlon

        # 浮点坐标减去左/上整数索引，得到 [0,1) 内的插值权重。
        row_weight = row - row0
        column_weight = column - column0

        # 先在 row0 上沿经度方向线性插值。
        upper = (
            field[row0, column0] * (1.0 - column_weight)
            + field[row0, column1] * column_weight
        )

        # 再在 row1 上沿经度方向线性插值。
        lower = (
            field[row1, column0] * (1.0 - column_weight)
            + field[row1, column1] * column_weight
        )

        # 最后沿纬度方向在 upper/lower 之间插值，得到标准双线性插值结果。
        return upper * (1.0 - row_weight) + lower * row_weight

    def _render_globe(self, diameter: int) -> QImage:
        assert self._frame is not None
        assert self._levels is not None

        # 构造渲染缓存 key。
        # 只有会影响“球面内部像素图像”的状态才放入 key；
        key = (
            self._frame.step,
            self._frame.name,
            self._level,
            round(self._yaw, 3),
            round(self._pitch, 3),
            diameter,
            self._colormap_name,
            self._levels,
        )

        # 当前状态与缓存完全一致时，直接复用上一张 QImage。
        # 这可避免窗口局部重绘时重复执行昂贵的 NumPy 投影与插值。
        if key == self._cache_key and self._cache_image is not None:
            return self._cache_image

        # 为直径为 diameter 的正方形图像构造像素中心坐标。
        # 最终 coordinate 大致均匀分布于 (-1, 1)，以球心为原点。
        coordinate = (
            np.arange(diameter, dtype=np.float32) + 0.5
        ) / diameter * 2.0 - 1.0

        # meshgrid 将一维坐标扩展为每个屏幕像素对应的二维 (x, y) 坐标。
        view_x, screen_y = np.meshgrid(coordinate, coordinate)

        # 图像 y 轴向下增加，而三维数学坐标通常 y 轴向上增加，所以这里翻转符号。
        view_y = -screen_y

        # 单位圆内的像素属于可见球面圆盘；圆外像素最终会设置为透明。
        radius_squared = view_x * view_x + view_y * view_y
        inside = radius_squared <= 1.0

        # 正射投影（orthographic projection）下，可见半球满足：
        # x^2 + y^2 + z^2 = 1 且 z >= 0。
        # 因此可由屏幕平面上的 x/y 反推出前半球 z。
        # maximum(0, ...) 防止圆外像素因负数开平方产生 NaN。
        view_z = np.sqrt(
            np.maximum(0.0, 1.0 - radius_squared),
            dtype=np.float32,
        )

        # 交互状态保存的是度，这里统一转换为弧度供 math.sin/cos 使用。
        yaw = math.radians(self._yaw)
        pitch = math.radians(self._pitch)
        cos_yaw, sin_yaw = math.cos(yaw), math.sin(yaw)
        cos_pitch, sin_pitch = math.cos(pitch), math.sin(pitch)

        # 从“观察者坐标/屏幕坐标”反变换回“世界球体坐标”。
        # 先恢复 pitch，再恢复 yaw。
        # 这样每个屏幕像素都能对应到球面上的一个三维单位向量。
        pitched_y = cos_pitch * view_y + sin_pitch * view_z
        pitched_z = -sin_pitch * view_y + cos_pitch * view_z
        world_x = cos_yaw * view_x + sin_yaw * pitched_z
        world_y = pitched_y
        world_z = -sin_yaw * view_x + cos_yaw * pitched_z

        # 三维单位球坐标 -> 地理纬经度：
        latitude = np.arcsin(np.clip(world_y, -1.0, 1.0))
        longitude = np.arctan2(world_x, world_z)

        # 取出当前层，并显式转成 float32，降低后续大数组运算的内存压力。
        field = np.asarray(self._frame.layer(self._level), dtype=np.float32)

        # 在球面每一个可见像素对应的纬经度位置，对模型规则网格做双线性采样。
        sampled = self._sample_field(field, latitude, longitude)

        # 按当前色阶范围归一化到 [0, 1]。
        lower, upper = self._levels
        finite = np.isfinite(sampled)
        normalized = np.zeros(sampled.shape, dtype=np.float32)
        normalized[finite] = np.clip(
            (sampled[finite] - lower) / (upper - lower), 0.0, 1.0
        )

        # [0,1] -> [0,255] 的整数索引，再查 256 级 RGB lookup table。
        # 查表后转为 float32，是因为下面还要乘阴影系数和做颜色混合。
        colors = self._lookup_table[
            np.asarray(normalized * 255.0, dtype=np.uint8)
        ].astype(np.float32)

        # limb shading：球心附近 view_z≈1，亮度系数≈1；
        # 靠近球体边缘 view_z≈0，亮度降到 0.84。
        # 这样提供轻微立体感，同时不显著改变气候场本身的色彩编码。
        shade = 0.84 + 0.16 * view_z[..., np.newaxis]
        colors *= shade

        # 经纬网每隔 30° 绘制一条线。
        grid_spacing = math.radians(30.0)

        # 计算每个像素纬度到最近 30° 倍数纬线的角距离。
        # “平移半个 spacing -> 取模 -> 再移回”的写法可以获得
        # 到最近周期网格线的对称距离。
        latitude_distance = np.abs(
            (latitude + 0.5 * grid_spacing) % grid_spacing - 0.5 * grid_spacing
        )

        # 同理计算到最近 30° 倍数经线的角距离。
        longitude_distance = np.abs(
            (longitude + 0.5 * grid_spacing) % grid_spacing - 0.5 * grid_spacing
        )

        # 距离任意经线/纬线小于 0.38° 即认为属于网格线；
        # 同时要求像素位于球面圆盘内部。
        grid = (
            (latitude_distance < math.radians(0.38))
            | (longitude_distance < math.radians(0.38))
        ) & inside

        # 网格像素与白色按 72% : 28% 混合，使网格变亮但不过分抢眼。
        colors[grid] = colors[grid] * 0.72 + 255.0 * 0.28

        # 构造最终 RGBA 图像缓冲区。
        rgba = np.zeros((diameter, diameter, 4), dtype=np.uint8)

        # RGB 写入前三个通道，并裁剪到合法的 uint8 范围。
        rgba[..., :3] = np.asarray(np.clip(colors, 0.0, 255.0), dtype=np.uint8)

        # 球面内部完全不透明，球面外部完全透明，形成圆形球体边界。
        rgba[..., 3] = np.where(inside & finite, 255, 0).astype(np.uint8)

        # QImage 需要稳定、连续的底层内存布局；这里确保 C contiguous。
        rgba = np.ascontiguousarray(rgba)

        # 用 NumPy 数组底层 buffer 构造 QImage。
        # rgba.strides[0] 是每一行占用的字节数（bytesPerLine）。
        image = QImage(
            rgba.data,
            diameter,
            diameter,
            rgba.strides[0],
            QImage.Format.Format_RGBA8888,
        ).copy()

        # 关键的 .copy() 会让 QImage 拥有自己的像素数据；
        # 否则局部变量 rgba 生命周期结束后，QImage 可能仍引用已经失效的 NumPy 内存。

        # 保存本次结果供下一次相同状态直接复用。
        self._cache_key = key
        self._cache_image = image
        return image

    def _coastline_path(self, target: QRectF) -> QPainterPath:
        """
        把地球表面的海岸线从“经纬度坐标”变成“三维球面坐标”，旋转地球，再把可见半球正射投影到二维屏幕。
        """
        # 海岸线与气候场必须使用完全一致的 yaw/pitch，才能与球面数据正确套准。
        yaw = math.radians(self._yaw)
        pitch = math.radians(self._pitch)
        cos_yaw, sin_yaw = math.cos(yaw), math.sin(yaw)
        cos_pitch, sin_pitch = math.cos(pitch), math.sin(pitch)

        # 把所有可见海岸线累积进一个 QPainterPath，最终一次性绘制。
        path = QPainterPath()

        # 每个 segment 是一条独立的经纬度折线，通常 shape 为 (N, 2)。
        for segment in coastline_segments():
            # Natural Earth 数据以“度”为单位，这里转成弧度。
            longitude = np.radians(segment[:, 0])
            latitude = np.radians(segment[:, 1])

            # 地理纬经度 -> 三维单位球笛卡尔坐标。
            # 本坐标约定中：y 指向北极，z 对应 0° 经线方向。
            cos_latitude = np.cos(latitude)
            world_x = cos_latitude * np.sin(longitude)
            world_y = np.sin(latitude)
            world_z = cos_latitude * np.cos(longitude)

            # 对世界坐标施加 yaw 旋转。
            # 这里是“世界 -> 观察者”的正向变换，因此符号与 _render_globe 中
            # “屏幕 -> 世界”的逆变换形式相对应。
            rotated_x = cos_yaw * world_x - sin_yaw * world_z
            rotated_z = sin_yaw * world_x + cos_yaw * world_z

            # 再施加 pitch，得到观察者坐标中的 y/z。
            view_y = cos_pitch * world_y - sin_pitch * rotated_z
            view_z = sin_pitch * world_y + cos_pitch * rotated_z

            # 正射投影下直接丢弃 z 深度，只使用 x/y 映射到 target 矩形。
            # rotated_x/view_y 的理论范围是 [-1, 1]。
            screen_x = target.left() + (rotated_x + 1.0) * 0.5 * target.width()
            # Qt 屏幕 y 轴向下，因此用 (1 - view_y) 完成翻转。
            screen_y = target.top() + (1.0 - view_y) * 0.5 * target.height()

            # drawing 表示上一点是否在可见半球且当前 path 正在连续绘制。
            drawing = False
            previous_x = 0.0
            previous_y = 0.0

            # view_z > 0 表示点位于面向观察者的前半球；
            # 背面的海岸线不应绘制。
            for x_value, y_value, visible in zip(screen_x, screen_y, view_z > 0.0):
                # NumPy scalar 转 Python float，交给 QPainterPath API。
                x = float(x_value)
                y = float(y_value)

                # 只有满足以下条件才把当前点与上一点连成线：
                # 1) 上一个点处于绘制状态；
                # 2) 当前点仍在可见半球；
                # 3) x/y 没出现异常的大跨度跳跃。
                continuous = (
                    drawing
                    and bool(visible)
                    and abs(x - previous_x) < target.width() * 0.35
                    and abs(y - previous_y) < target.height() * 0.35
                )

                # 连续可见时扩展当前子路径。
                if continuous:
                    path.lineTo(x, y)
                # 如果当前点可见但不能与上一点连接，则从这里开始新的子路径。
                elif visible:
                    path.moveTo(x, y)

                # 下一轮是否允许延续线段取决于当前点是否可见。
                drawing = bool(visible)
                previous_x = x
                previous_y = y

        return path

    def paintEvent(self, _event: QPaintEvent) -> None:
        # Qt 要求所有 QWidget 自绘内容都在 paintEvent 中通过 QPainter 完成。
        painter = QPainter(self)

        # 折线/椭圆启用抗锯齿；缩放 QImage 时启用平滑插值。
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)

        # 每次重绘先清成绘图面板的深色背景，避免上一帧残留。
        painter.fillRect(self.rect(), QColor(PLOT_BACKGROUND))

        # 在还没有模型数据时只显示等待提示，不尝试渲染球体。
        if self._frame is None:
            painter.setPen(QColor(TEXT_MUTED))
            font = painter.font()
            font.setPointSize(13)
            painter.setFont(font)
            painter.drawText(
                self.rect(),
                Qt.AlignmentFlag.AlignCenter,
                self._empty_message,
            )
            return

        # 为标题和底部交互提示预留空间后，取宽高中的较小值作为球体基准尺寸。
        available = min(self.width() - 42, self.height() - 72)

        # 应用用户滚轮控制的缩放；同时保证球体直径至少为 96 像素。
        diameter = max(96, int(available * self._zoom))

        # 水平居中。
        left = (self.width() - diameter) // 2

        # 垂直大致居中，但至少从 y=34 开始，给顶部标题留出空间。
        top = max(34, (self.height() - diameter) // 2 - 4)

        # target 同时用于绘制球体图像、海岸线投影以及球体轮廓。
        target = QRectF(left, top, diameter, diameter)

        # 创建球体下方的径向渐变阴影，提供悬浮/立体效果。
        # 渐变中心略低于球体底边，半径约为球直径的 58%。
        shadow = QRadialGradient(
            target.center().x(),
            target.bottom() + diameter * 0.045,
            diameter * 0.58,
        )

        # 从中心向外逐步减小 alpha，形成柔和阴影。
        shadow.setColorAt(0.0, QColor(0, 0, 0, 50))
        shadow.setColorAt(0.72, QColor(0, 0, 0, 18))
        shadow.setColorAt(1.0, QColor(0, 0, 0, 0))

        # 阴影椭圆只填充，不描边。
        painter.setBrush(shadow)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawEllipse(
            QRectF(
                left - diameter * 0.09,
                target.bottom() - diameter * 0.05,
                diameter * 1.18,
                diameter * 0.20,
            )
        )

        # 先画气候场球体像素图。
        # _render_globe 内部有状态缓存，相同视角/帧/尺寸下不会重复 NumPy 渲染。
        painter.drawImage(target, self._render_globe(diameter))

        # 根据当前视角把海岸线转换成一个屏幕空间 QPainterPath。
        coastline_path = self._coastline_path(target)
        painter.setBrush(Qt.BrushStyle.NoBrush)

        # 第一遍用较粗、半透明白线画“外描边/光晕”，提升深浅底色上的可辨识度。
        painter.setPen(
            QPen(
                QColor(255, 255, 255, 205),
                2.4,
                Qt.PenStyle.SolidLine,
                Qt.PenCapStyle.RoundCap,
                Qt.PenJoinStyle.RoundJoin,
            )
        )
        painter.drawPath(coastline_path)

        # 第二遍叠加较细的深色线，形成清晰的海岸线主体。
        painter.setPen(
            QPen(
                QColor(29, 29, 31, 195),
                0.9,
                Qt.PenStyle.SolidLine,
                Qt.PenCapStyle.RoundCap,
                Qt.PenJoinStyle.RoundJoin,
            )
        )
        painter.drawPath(coastline_path)

        # 最后为球体外缘加一圈极淡的边界线。
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.setPen(QPen(QColor(255, 255, 255, 45), 1.2))
        painter.drawEllipse(target)

        # ---------- 顶部标题 ----------
        painter.setPen(QColor(TEXT))
        title_font = painter.font()
        title_font.setPointSize(12)
        title_font.setWeight(title_font.Weight.DemiBold)
        painter.setFont(title_font)
        painter.drawText(
            QRectF(18, 10, self.width() - 36, 24),
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
            f"{self._frame.name}  ·  {self._frame.unit}"
            f"  ·  step {self._frame.step}",
        )

        # ---------- 底部交互提示 ----------
        painter.setPen(QColor(TEXT_MUTED))
        hint_font = painter.font()
        hint_font.setPointSize(10)
        hint_font.setWeight(hint_font.Weight.Normal)
        painter.setFont(hint_font)
        painter.drawText(
            QRectF(18, self.height() - 28, self.width() - 36, 18),
            Qt.AlignmentFlag.AlignCenter,
            "Drag to rotate  ·  Scroll to zoom  ·  Double-click to reset",
        )

    def mousePressEvent(self, event: QMouseEvent) -> None:
        # 只处理左键拖拽；其他鼠标键交给默认 Qt 行为。
        if event.button() == Qt.MouseButton.LeftButton:
            # 记录拖拽起点，后续 mouseMoveEvent 用相邻位置差计算旋转量。
            self._drag_position = event.position().toPoint()

            # 拖动过程中切换成“握住”光标，提供直接操作反馈。
            self.setCursor(Qt.CursorShape.ClosedHandCursor)

            # 标记事件已经由本控件处理。
            event.accept()

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        # 没有按下左键开始拖拽时，不执行旋转。
        if self._drag_position is None:
            inspected = self._inspect_position(event.position().x(), event.position().y())
            if inspected is not None:
                latitude, longitude, value = inspected
                text = f"{latitude:.2f}°, {longitude:.2f}° · {value:.5g} {self._frame.unit if self._frame else ''}"
                QToolTip.showText(event.globalPosition().toPoint(), text, self)
            return

        # 当前鼠标整数像素坐标。
        position = event.position().toPoint()

        # 与上一采样点的位移，而不是与最初按下位置的总位移。
        delta = position - self._drag_position

        # 更新“上一位置”，下一次移动继续按增量计算。
        self._drag_position = position

        # 水平拖动控制 yaw。
        # x 向右增大时这里做减法，因此球体视觉旋转方向与拖拽方向相协调。
        # % 360 保证 yaw 始终归一到 [0, 360)。
        self._yaw = (self._yaw - delta.x() * 0.55) % 360.0

        # 垂直拖动控制 pitch。
        # 限制在 [-85°, 85°]，避免接近 ±90° 时出现视角翻转/交互退化。
        self._pitch = max(-85.0, min(85.0, self._pitch + delta.y() * 0.45))

        # 视角变化后缓存图像失效，请求重新绘制。
        self._invalidate()
        event.accept()

    def _inspect_position(self, x: float, y: float) -> tuple[float, float, float] | None:
        if self._frame is None:
            return None
        available = min(self.width() - 42, self.height() - 72)
        diameter = max(96, int(available * self._zoom))
        left = (self.width() - diameter) / 2.0
        top = max(34.0, (self.height() - diameter) / 2.0 - 4.0)
        view_x = (x - left) / diameter * 2.0 - 1.0
        view_y = -((y - top) / diameter * 2.0 - 1.0)
        radius_squared = view_x * view_x + view_y * view_y
        if radius_squared > 1.0:
            return None
        view_z = math.sqrt(max(0.0, 1.0 - radius_squared))
        yaw, pitch = math.radians(self._yaw), math.radians(self._pitch)
        pitched_y = math.cos(pitch) * view_y + math.sin(pitch) * view_z
        pitched_z = -math.sin(pitch) * view_y + math.cos(pitch) * view_z
        world_x = math.cos(yaw) * view_x + math.sin(yaw) * pitched_z
        world_z = -math.sin(yaw) * view_x + math.cos(yaw) * pitched_z
        latitude = math.asin(max(-1.0, min(1.0, pitched_y)))
        longitude = math.atan2(world_x, world_z)
        sample = self._sample_field(
            np.asarray(self._frame.layer(self._level), dtype=np.float32),
            np.asarray(latitude),
            np.asarray(longitude),
        )
        value = float(sample)
        if not math.isfinite(value):
            return None
        return math.degrees(latitude), math.degrees(longitude), value

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        # 左键释放即结束拖拽。
        if event.button() == Qt.MouseButton.LeftButton:
            self._drag_position = None

            # 恢复空闲状态的“张开手”光标。
            self.setCursor(Qt.CursorShape.OpenHandCursor)
            event.accept()

    def wheelEvent(self, event: QWheelEvent) -> None:
        # angleDelta().y() 是 Qt 的滚轮增量。
        # 除以 2400 把滚轮增量转换成一个较平缓的乘法缩放因子。
        factor = 1.0 + event.angleDelta().y() / 2400.0

        self._zoom = max(0.60, min(1.16, self._zoom * factor))

        # 缩放会改变实际渲染直径，因此清除缓存并重绘。
        self._invalidate()
        event.accept()

    def mouseDoubleClickEvent(self, event: QMouseEvent) -> None:
        # 左键双击恢复默认 yaw/pitch/zoom。
        if event.button() == Qt.MouseButton.LeftButton:
            self.reset_orientation()
            event.accept()

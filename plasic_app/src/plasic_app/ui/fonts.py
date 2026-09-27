from __future__ import annotations

# QFont：
#   Qt 中用于描述字体属性的类，例如字体家族、字号、粗细等。
#
# QFontDatabase：
#   Qt 提供的字体数据库接口
from PySide6.QtGui import QFont, QFontDatabase


def interface_font(point_size: int = 13) -> QFont:
    """
    返回一个当前系统已经安装的无衬线字体，用于界面 UI。

    参数：
        point_size:
            字体大小，单位为 point（磅）。
            默认值为 13。

    返回：
        QFont:
            优先返回候选字体列表中第一个已安装的字体；
            如果候选字体全部不存在，则返回 Qt 默认字体，
            并将其字号设置为 point_size。

    这里主动检查字体是否存在的目的：
        避免直接指定一个系统中不存在的字体后，
        由 Qt 自动进行字体 alias / fallback 时产生相关警告。
    """

    for family in (
        "Helvetica Neue",
        "Arial",
        "Helvetica",
        "DejaVu Sans",
    ):
        # 检查当前系统的字体数据库中是否存在指定字体 family。
        if QFontDatabase.hasFamily(family):
            return QFont(family, point_size)

    # 如果上面的候选字体全部没有安装，
    font = QFont()

    font.setPointSize(point_size)

    return font


def fixed_width_font(point_size: int = 11) -> QFont:
    """
    返回一个当前系统已经安装的等宽字体，主要用于 process log。

    等宽字体（fixed-width / monospace）的特点是：
        每个字符占据相同的水平宽度。

    因此比较适合显示：
        - 日志；
        - 终端输出；
        - 代码；
        - 对齐后的数值或文本。

    参数：
        point_size:
            字体大小，单位为 point（磅）。
            默认值为 11。

    返回：
        QFont:
            优先返回候选等宽字体中第一个已安装的字体；
            如果都不存在，则使用 Qt 提供的系统默认 FixedFont。
    """

    for family in ("Menlo", "Monaco", "DejaVu Sans Mono", "Courier New"):

        if QFontDatabase.hasFamily(family):

            return QFont(family, point_size)

    font = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)

    font.setPointSize(point_size)

    return font

from __future__ import annotations

# importlib.resources：
#   Python 标准库中读取包内资源的接口。
#   files()/as_file() 的组合让同一份代码既能从普通源码目录加载资源，
#   也能从打包后的 Python package（如 wheel/zip）中加载资源。
#
# QIcon / QPixmap：
#   QPixmap 表示一张可供 Qt 绘制的位图；
#   QIcon 表示可以同时适配窗口、Dock、任务栏等场景的图标对象。
from importlib.resources import as_file, files

from PySide6.QtGui import QIcon, QPixmap

# 品牌 Logo 在 plasic_app.resources 包中的文件名。
# pyproject.toml 的 package-data 已经把 resources/*.png 一并打包，
# 因此安装版和源码方式启动都能找到同一张 Logo。
LOGO_RESOURCE_NAME = "logo.png"


def logo_pixmap() -> QPixmap:
    """
    从包资源中加载 PlaSiC Logo，并返回 QPixmap。

    返回：
        QPixmap:
            已解码的 Logo 位图。with 语句结束后文件句柄会关闭，
            但位图数据已经完整读入内存，可以继续使用。
    """

    logo_resource = files("plasic_app.resources").joinpath(LOGO_RESOURCE_NAME)

    with as_file(logo_resource) as logo_path:
        return QPixmap(str(logo_path))


def logo_icon() -> QIcon:
    """
    返回应用级 Logo 图标。

    macOS 上把该 QIcon 交给 QApplication.setWindowIcon() 后，
    从终端启动的 Python 进程也会在 Dock 中显示 PlaSiC Logo。
    """

    return QIcon(logo_pixmap())

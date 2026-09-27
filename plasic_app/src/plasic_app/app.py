from __future__ import annotations
import sys

def main() -> int:
    """
    PlaSiC 图形界面的程序入口函数。

    返回值：
        int:
            Qt 应用程序退出时返回的状态码。
            如果 PySide6 没有安装，则提前返回 2。
    """

    # 尝试导入 Qt 应用程序的核心类 QApplication。
    #
    # QApplication 是所有 Qt Widgets GUI 应用程序的基础对象，
    # 它负责管理：
    # - GUI 应用程序的生命周期；
    # - Qt 事件循环；
    # - 鼠标、键盘等输入事件；
    # - 全局字体、样式、调色板；
    # - 窗口之间的一些公共资源。
    #
    # 将这个 import 放在 try 块中，而不是放到文件最顶部，
    # 可以在用户没有安装 PySide6 时给出更明确、更友好的错误信息。
    try:
        from PySide6.QtGui import QSurfaceFormat
        from PySide6.QtWidgets import QApplication

    except ImportError:
        print(
            "PySide6 is not installed. Install plasic_app with "
            "`python -m pip install -e plasic_app`.",
            file=sys.stderr,
        )

        return 2

    # pyqtgraph's volume shaders need a modern OpenGL context on macOS. The
    # default format has to be selected before QApplication is constructed.
    if sys.platform == "darwin":
        surface_format = QSurfaceFormat()
        surface_format.setRenderableType(QSurfaceFormat.RenderableType.OpenGL)
        surface_format.setProfile(QSurfaceFormat.OpenGLContextProfile.CoreProfile)
        surface_format.setVersion(4, 1)
        surface_format.setDepthBufferSize(24)
        QSurfaceFormat.setDefaultFormat(surface_format)

    # 导入 PlaSiC 应用程序的主窗口类。
    from plasic_app.ui.main_window import MainWindow

    # 导入应用级 Logo 图标。
    # 在 macOS 上，QApplication 的窗口图标会同时作为 Dock 栏图标显示。
    from plasic_app.ui.branding import logo_icon

    # 导入应用程序统一使用的界面字体配置函数。
    # 后面通过 application.setFont(...) 设置为整个 Qt 应用的默认字体。
    from plasic_app.ui.fonts import interface_font

    # 导入应用程序的浅色调色板配置函数。
    # 后面会把它应用到 QApplication，
    # 从而统一控制控件背景色、文本色、按钮色等 Qt Palette 属性。
    from plasic_app.ui.theme import light_palette

    # 创建 QApplication 实例。
    # 一个基于 Qt Widgets 的 GUI 程序通常只能创建一个 QApplication。
    # sys.argv 会把当前程序的命令行参数传递给 Qt，
    # Qt 自身可以解析其中部分与 GUI 平台、显示方式等相关的参数。
    application = QApplication(sys.argv)

    # 设置应用程序名称为 "PlaSiC"。
    application.setApplicationName("PlaSiC")

    # 设置组织名称为 "PlaSiC"。
    #
    # Qt 的 QSettings 等组件通常会结合
    # organizationName 和 applicationName
    # 来确定应用程序配置项的存储位置。
    application.setOrganizationName("PlaSiC")

    # 设置整个 Qt 应用的默认窗口图标。
    #
    # 该设置必须放在创建和显示窗口之前：
    # - Windows / Linux 上，它决定任务栏和窗口标题栏使用的图标；
    # - macOS 上，非 App Bundle 方式（例如从终端执行
    #   `PYTHONPATH=src python -m plasic_app`）启动时，
    #   Qt 会把这里设置的图标同步为 Dock 栏图标，
    #   从而取代 Python 解释器默认的图标。
    application.setWindowIcon(logo_icon())

    # 设置 Qt GUI 的全局 Widget 样式为 "Fusion"。
    #
    # Fusion 是 Qt 提供的一套跨平台样式，
    # 相比完全依赖操作系统原生样式，
    # 它通常能够在 Windows、Linux、macOS 等平台上
    # 提供更加统一的控件外观。
    application.setStyle("Fusion")

    # 设置整个 QApplication 的默认字体。
    # interface_font() 是项目内部定义的字体配置函数，
    # 预计返回一个 QFont 对象。
    # 设置在 QApplication 层级意味着，
    # 如果某个子控件没有显式指定自己的字体，
    # 通常会继承这里设置的字体。
    application.setFont(interface_font())

    # 设置整个 Qt 应用程序使用的颜色调色板。
    # light_palette() 是项目内部定义的浅色主题配置函数，
    # 预计返回一个 QPalette 对象。
    # 这会影响大量标准 Qt Widget 的默认颜色表现，
    application.setPalette(light_palette())

    # 实例化 PlaSiC 的主窗口。
    window = MainWindow()

    # 显示主窗口。
    #
    # Qt Widget 创建后默认并不一定可见，
    # 调用 show() 后才会请求窗口系统将其显示出来。
    window.show()

    # 启动 Qt 主事件循环。
    #
    # 从这里开始，Qt 会持续处理：
    # - 鼠标事件；
    # - 键盘事件；
    # - 窗口重绘；
    # - 信号与槽；
    # - 定时器；
    # - 系统窗口事件等。
    #
    # application.exec() 一般会一直阻塞，
    # 直到应用程序退出（例如关闭最后一个窗口或显式调用 quit）。
    #
    # exec() 最终返回 Qt 的退出状态码，
    # main() 将该状态码直接返回给调用者。
    return application.exec()

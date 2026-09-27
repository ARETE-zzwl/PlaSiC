# QColor:
#   Qt 中表示颜色的类，可以通过十六进制字符串（如 "#FFFFFF"）构造颜色。
#
# QPalette:
#   Qt 的调色板类，用来定义窗口、文本、按钮、输入框、选中区域等
#   Qt 标准控件在不同状态下使用的基础颜色。
from PySide6.QtGui import QColor, QPalette

# =============================================================================
# 全局颜色常量
# =============================================================================
# 这里集中定义整个应用的基础视觉色板。
#
# 设计思路：
# - 整体使用浅色（日间）主题。
# - 背景以低饱和灰蓝色为主。
# - 内容卡片保持纯白。
# - 蓝绿色 ACCENT 仅用于：
#       * 主按钮
#       * 输入框焦点
#       * 进度条
#       * 链接
#       * 状态提示
#   从而避免 GUI 本身和科学数据可视化中的 colormap 争夺视觉注意力。
#
# 统一使用常量还有一个重要好处：
# 后续若需要调整主题，只需修改这里，而不需要逐个修改 Python 代码中的颜色。

# 应用最底层背景色。
# 比纯白稍深，可以让白色 Card / Header 在背景上形成轻微层次感。
BACKGROUND = "#F2F5F7"

# 常规内容表面颜色。
# 主要用于输入框 Base、Card 等核心内容区域。
SURFACE = "#FFFFFF"

# 比主 Surface 略深的次级表面颜色。
# 用于 AlternateBase 或视觉层级稍低的区域。
SURFACE_RAISED = "#F7F9FA"

# 绘图区域背景。
# 单独定义是为了未来可以独立调整 plot canvas，
# 而不必与普通 UI Surface 强绑定。
PLOT_BACKGROUND = "#FFFFFF"

# 常规边框颜色。
# 使用非常浅的灰蓝色，使卡片边界可见但不过分突出。
BORDER = "#DCE2E7"

# 主文字颜色。
# 使用深蓝灰而不是纯黑，可以降低长时间查看界面的视觉疲劳。
TEXT = "#16212B"

# 次级文字颜色。
# 用于 subtitle、说明信息、辅助标签等。
TEXT_MUTED = "#697586"

# 主强调色。
# 用于焦点、高亮、链接、进度条、主要操作按钮等。
ACCENT = "#197FA3"

# ACCENT 的 hover / 强调版本。
# 相比 ACCENT 更深，用于鼠标悬停时提供明确反馈。
ACCENT_HOVER = "#126A89"


def light_palette() -> QPalette:
    """
    创建并返回一个完整的浅色 QPalette。

    这个函数的核心目的不是简单地覆盖几个颜色，而是显式设置
    Active / Inactive / Disabled 三组颜色，使应用界面尽可能不受
    Windows、macOS、Linux 系统主题以及深色模式设置影响。

    Returns
    -------
    QPalette
        已配置完整浅色主题颜色的 Qt 调色板。
    """

    # 创建一个新的调色板实例。
    # 后续通过 setColor() 为不同 ColorGroup / ColorRole 设置颜色。
    palette = QPalette()

    # Qt 会根据窗口/控件当前状态选择不同的 ColorGroup：
    #
    # Active:
    #   当前窗口拥有焦点时使用。
    #
    # Inactive:
    #   窗口存在但当前没有获得系统焦点时使用。
    #
    # 这里让 Active 和 Inactive 使用完全相同的配色，
    # 这样窗口失去焦点后不会因为系统主题而出现明显颜色变化。
    active_groups = (
        QPalette.ColorGroup.Active,
        QPalette.ColorGroup.Inactive,
    )

    # -------------------------------------------------------------------------
    # Qt 标准 ColorRole -> 项目颜色映射
    # -------------------------------------------------------------------------
    #
    # QPalette.ColorRole 表示颜色的“用途”，而不是具体某个 Widget。
    # 例如：
    #
    # Window      -> 普通窗口背景
    # WindowText  -> 窗口上的普通文字
    # Base        -> 输入框、列表等主要内容区背景
    # Text        -> 输入控件中的文字
    # Button      -> 按钮背景
    # ButtonText  -> 按钮文字
    # Highlight   -> 被选中内容的背景
    colors = {
        # 普通窗口背景。
        QPalette.ColorRole.Window: BACKGROUND,
        # 窗口级普通文本。
        QPalette.ColorRole.WindowText: TEXT,
        # 文本框、列表、表格等控件的主要内容背景。
        QPalette.ColorRole.Base: SURFACE,
        # 交替行等次级背景。
        QPalette.ColorRole.AlternateBase: SURFACE_RAISED,
        # Tooltip 使用深色背景。
        # 这是整个浅色主题中少数采用深底浅字的区域。
        QPalette.ColorRole.ToolTipBase: "#26333D",
        # Tooltip 的白色文字。
        QPalette.ColorRole.ToolTipText: "#FFFFFF",
        # 输入框、文本编辑器等内容文字。
        QPalette.ColorRole.Text: TEXT,
        # Qt 标准按钮背景。
        QPalette.ColorRole.Button: SURFACE,
        # Qt 标准按钮文字。
        QPalette.ColorRole.ButtonText: TEXT,
        # 强亮文字。
        # 通常用于需要在深色背景上显示的文本。
        QPalette.ColorRole.BrightText: "#FFFFFF",
        # 超链接颜色。
        QPalette.ColorRole.Link: ACCENT,
        # 文本选择、列表选择等选中状态背景。
        QPalette.ColorRole.Highlight: ACCENT,
        # 选中状态下的文字颜色。
        QPalette.ColorRole.HighlightedText: "#FFFFFF",
        # 输入框 placeholder / hint 文本颜色。
        QPalette.ColorRole.PlaceholderText: "#8B96A3",
        # Qt 传统 3D palette 体系中的高亮边缘。
        QPalette.ColorRole.Light: "#FFFFFF",
        # Light 与 Mid 之间的中间亮度。
        QPalette.ColorRole.Midlight: "#E8EDF1",
        # Qt 传统边框/阴影体系中的中间色。
        QPalette.ColorRole.Mid: "#D5DCE2",
        # 较深的边缘颜色。
        QPalette.ColorRole.Dark: "#97A2AE",
        # 最深层阴影颜色。
        QPalette.ColorRole.Shadow: "#1C2730",
    }

    # -------------------------------------------------------------------------
    # 设置 Active 和 Inactive 状态
    # -------------------------------------------------------------------------
    #
    # 两种状态采用同一套颜色。
    # 这样应用在失焦时不会突然切换到操作系统默认配色。
    for group in active_groups:
        for role, color in colors.items():
            # QColor 将 "#RRGGBB" 字符串转换成 Qt 颜色对象。
            palette.setColor(group, role, QColor(color))

    # -------------------------------------------------------------------------
    # Disabled 状态
    # -------------------------------------------------------------------------
    # 先完整复制正常颜色，然后针对文字、输入区、按钮等关键角色进行覆盖。
    # 这样可以确保 Disabled 状态仍然拥有所有 ColorRole，
    # 同时在视觉上明显表现出“不可操作”状态。
    disabled = QPalette.ColorGroup.Disabled

    # 首先给 Disabled 状态设置完整基础色板。
    for role, color in colors.items():
        palette.setColor(disabled, role, QColor(color))

    # Disabled 状态下：
    # Text       -> 输入控件文字
    # WindowText -> 普通标签等文字
    # ButtonText -> 按钮文字
    # 都统一改成较浅的灰色。
    for role in (
        QPalette.ColorRole.Text,
        QPalette.ColorRole.WindowText,
        QPalette.ColorRole.ButtonText,
    ):
        palette.setColor(disabled, role, QColor("#9AA4AE"))

    # Disabled 输入区域背景。
    # 比正常 SURFACE 稍暗，用于强调当前控件不可编辑。
    palette.setColor(
        disabled,
        QPalette.ColorRole.Base,
        QColor("#F1F4F6"),
    )

    # Disabled 按钮背景。
    palette.setColor(
        disabled,
        QPalette.ColorRole.Button,
        QColor("#F4F6F8"),
    )

    # 返回构造完成的 QPalette。
    return palette


# =============================================================================
# 应用级 Qt Style Sheet
# =============================================================================
#
# STUDIO_STYLESHEET 是整个 GUI 的 QSS（Qt Style Sheet）。
#
# QSS 的语法与 CSS 类似，但选择器对应的是 Qt Widget：
#
#   QWidget
#       所有 QWidget 及其派生控件。
#
#   QWidget#Root
#       objectName == "Root" 的 QWidget。
#
#   QLabel[role="secondary"]
#       Qt dynamic property "role" == "secondary" 的 QLabel。
#
#   QPushButton:hover
#       鼠标悬停状态。
#
#   QPushButton:disabled
#       disabled 状态。
#
#   QPushButton[role="primary"]
#       dynamic property role="primary" 的按钮。
#
#
# -----------------------------------------------------------------------------
# 下面 QSS 各部分的作用说明
# -----------------------------------------------------------------------------
#
# 1. QWidget
#    设置整个应用的默认文字颜色和字号。
#
# 2. QMainWindow, QWidget#Root
#    设置应用最外层灰蓝色背景。
#
# 3. ParameterForm / ParameterContent / ParameterViewport / ParameterScroll
#    参数面板中的中间容器保持 transparent，
#    避免多层 QWidget 自己绘制背景导致出现色块。
#
# 4. HeaderBar / Card
#    统一使用：
#       - 白色背景
#       - 浅灰边框
#       - 16 px 圆角
#    从而形成主要页面容器。
#
# 5. AppLogo
#    固定 Logo 为 46 × 44，并取消背景和边框。
#
# 6. AppTitle
#    应用主标题，22 px + 700 font-weight。
#
# 7. AppSubtitle / role="secondary"
#    次级说明文本，颜色更浅、字号更小。
#
# 8. QLabel
#    普通 Label 的默认颜色。
#
# 9. role="section"
#    Section 标题，15 px、semi-bold。
#
# 10. role="eyebrow"
#     小型强调标签，使用 ACCENT 色和大写标签式视觉。
#
# 11. SourceBadge
#     信息来源/数据源一类的小型 Badge。
#
# 12. StatePill
#     状态胶囊标签。
#     默认是蓝色信息状态。
#
#     state="running":
#         绿色，表示运行中。
#
#     state="building" / state="stopping":
#         黄色，表示过渡状态。
#
# 13. HelpButton
#     帮助按钮采用弱视觉层级：
#     默认透明背景，只保留轻边框。
#
# 14. QLineEdit / QComboBox / QSpinBox / QDoubleSpinBox
#     所有主要表单输入控件统一：
#         - 31 px 最小高度
#         - 8 px 圆角
#         - 浅色背景
#         - ACCENT selection
#
# 15. 输入控件 hover
#     鼠标经过时边框变深，但不改变背景。
#
# 16. 输入控件 disabled
#     使用灰色文字 + 灰色背景。
#
# 17. 输入控件 focus
#     使用 ACCENT 蓝色边框突出当前输入焦点。
#
# 18. QComboBox::drop-down
#     移除 Qt 默认下拉按钮边框，只保留 26 px 宽度区域。
#
# 19. QAbstractItemView
#     ComboBox 下拉菜单，以及其他 Item View 的基础样式。
#
# 20. QCheckBox
#     checkbox label 与 indicator 之间保持 8 px 间距。
#
# 21. QCheckBox::indicator
#     自定义 checkbox 方框尺寸、圆角、背景和边框。
#
# 22. checkbox hover
#     边框改为强调蓝。
#
# 23. checkbox checked
#     使用 ACCENT 蓝色填充。
#
# 24. checkbox disabled
#     使用灰色弱化状态。
#
# 25. QPushButton
#     普通按钮统一基础样式。
#
# 26. QPushButton:hover
#     轻微灰色背景 + 较深边框。
#
# 27. QPushButton:pressed
#     使用更明显的灰蓝色表示按下反馈。
#
# 28. QPushButton:disabled
#     灰化不可用按钮。
#
# 29. role="primary"
#     主要 CTA 按钮。
#     使用 ACCENT 蓝底白字。
#
# 30. primary:hover
#     使用更深的 ACCENT_HOVER。
#
# 31. role="danger"
#     删除、终止等危险操作按钮。
#     使用低饱和红色，而不是高饱和纯红。
#
# 32. role="quick"
#     快捷操作按钮，视觉权重低于 primary。
#
# 33. segment="true"
#     Segmented Control 内部按钮。
#     默认透明，checked 项使用白色背景。
#
# 34. SegmentedControl
#     Segment 按钮的整体容器背景和边框。
#
# 35. QGroupBox
#     参数分组容器。
#     使用浅色背景、轻边框和顶部标题。
#
# 36. QGroupBox::title
#     对 GroupBox 标题本身进行定位和样式控制。
#
# 37. QScrollArea
#     去掉 Qt 默认 ScrollArea 边框和背景。
#
# 38. vertical / horizontal QScrollBar
#     自定义滚动条，使其更轻量。
#
# 39. QProgressBar
#     6 px 高的细进度条。
#
# 40. QProgressBar::chunk
#     进度部分使用 ACCENT 蓝。
#
# 41. QTextEdit
#     文本日志/说明区域使用比 Card 稍暗的背景。
#
# 42. QSplitter::handle
#     Splitter handle 本身透明，但保留 10 px 可拖拽区域。
#
# 43. QToolTip
#     深色 Tooltip，和 light palette 中的 ToolTipBase 相呼应。
#
# 44. QMessageBox
#     强制消息框使用白色背景，避免被系统深色主题影响。
#
#
# 注意：
# 下面字符串保持原始内容不变。
# 这是刻意的，因为如果把 Python "#" 注释放进字符串，它会成为 QSS 内容；
# 如果加入 "/* ... */"，虽然 Qt 通常会把它视为注释，但仍然改变了
# STUDIO_STYLESHEET 的实际字符串值。
STUDIO_STYLESHEET = """
QWidget {
    color: #16212B;
    font-size: 13px;
}
QMainWindow, QWidget#Root {
    background: #F2F5F7;
}
QWidget#ParameterForm,
QWidget#ParameterContent,
QWidget#ParameterViewport,
QScrollArea#ParameterScroll {
    background: transparent;
}
QFrame#HeaderBar {
    background: #FFFFFF;
    border: 1px solid #DCE2E7;
    border-radius: 16px;
}
QFrame#Card {
    background: #FFFFFF;
    border: 1px solid #DCE2E7;
    border-radius: 16px;
}
QLabel#AppLogo {
    min-width: 46px;
    max-width: 46px;
    min-height: 44px;
    max-height: 44px;
    background: transparent;
    border: none;
}
QLabel#AppTitle {
    color: #14212B;
    font-size: 22px;
    font-weight: 700;
}
QLabel#AppSubtitle, QLabel[role="secondary"] {
    color: #697586;
    font-size: 12px;
}
QLabel {
    color: #24313C;
    background: transparent;
}
QLabel[role="section"] {
    color: #17232D;
    font-size: 15px;
    font-weight: 600;
}
QLabel[role="eyebrow"] {
    color: #197FA3;
    font-size: 10px;
    font-weight: 700;
}
QLabel#SourceBadge {
    color: #62707D;
    background: #F5F8F9;
    border: 1px solid #E0E6EA;
    border-radius: 10px;
    padding: 5px 9px;
    font-size: 11px;
}
QLabel#StatePill {
    background: #ECF6FA;
    color: #146E8D;
    border: 1px solid #C9E4ED;
    border-radius: 12px;
    padding: 5px 12px;
    font-size: 11px;
    font-weight: 600;
}
QLabel#StatePill[state="running"] {
    background: #EAF7EF;
    color: #257248;
    border-color: #CBE8D5;
}
QLabel#StatePill[state="building"],
QLabel#StatePill[state="stopping"] {
    background: #FFF6E5;
    color: #8B641A;
    border-color: #F0DDAF;
}
QPushButton#HelpButton {
    min-height: 28px;
    padding: 2px 14px;
    color: #53616E;
    background: transparent;
    border: 1px solid #D5DDE3;
    border-radius: 9px;
    font-weight: 500;
}
QPushButton#HelpButton:hover {
    color: #1D2B35;
    background: #F5F8FA;
    border-color: #BFC9D1;
}
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {
    min-height: 31px;
    background: #FBFCFD;
    border: 1px solid #D4DCE2;
    border-radius: 8px;
    padding: 2px 9px;
    color: #1B2731;
    selection-background-color: #197FA3;
    selection-color: #FFFFFF;
}
QLineEdit:hover, QComboBox:hover, QSpinBox:hover,
QDoubleSpinBox:hover {
    border-color: #AEBAC4;
}
QLineEdit:disabled, QComboBox:disabled,
QSpinBox:disabled, QDoubleSpinBox:disabled {
    background: #F1F4F6;
    color: #98A2AD;
    border-color: #E2E7EB;
}
QLineEdit:focus, QComboBox:focus, QSpinBox:focus,
QDoubleSpinBox:focus {
    border: 1px solid #197FA3;
}
QComboBox::drop-down {
    width: 26px;
    border: none;
}
QComboBox QAbstractItemView, QAbstractItemView {
    background: #FFFFFF;
    color: #1D2933;
    selection-background-color: #E2F1F5;
    selection-color: #163744;
    border: 1px solid #CBD4DB;
    outline: none;
    padding: 4px;
}
QCheckBox {
    spacing: 8px;
    color: #2A3742;
}
QCheckBox::indicator {
    width: 16px;
    height: 16px;
    border: 1px solid #AAB6C0;
    border-radius: 5px;
    background: #FFFFFF;
}
QCheckBox::indicator:hover {
    border-color: #328EAD;
}
QCheckBox::indicator:checked {
    background: #197FA3;
    border-color: #197FA3;
}
QCheckBox::indicator:disabled {
    background: #F0F3F5;
    border-color: #D7DEE3;
}
QPushButton {
    min-height: 31px;
    background: #FFFFFF;
    color: #24313C;
    border: 1px solid #D1D9DF;
    border-radius: 8px;
    padding: 2px 14px;
    font-weight: 500;
}
QPushButton:hover {
    background: #F5F8FA;
    border-color: #ABB8C2;
}
QPushButton:pressed {
    background: #EAF0F3;
}
QPushButton:disabled {
    color: #A1AAB3;
    background: #F4F6F8;
    border-color: #E1E6EA;
}
QPushButton[role="primary"] {
    color: #FFFFFF;
    background: #197FA3;
    border-color: #197FA3;
    font-weight: 700;
}
QPushButton[role="primary"]:hover {
    background: #126A89;
    border-color: #126A89;
}
QPushButton[role="danger"] {
    color: #B74A4A;
    background: #FFF7F7;
    border-color: #F0CFCF;
}
QPushButton[role="danger"]:hover {
    color: #A03838;
    background: #FFF0F0;
    border-color: #E5B5B5;
}
QPushButton[role="quick"] {
    color: #4F5D69;
    background: #F7F9FA;
    border-color: #DCE2E7;
    font-size: 12px;
}
QPushButton[role="quick"]:hover {
    color: #126A89;
    background: #EBF6F9;
    border-color: #B8DCE7;
}
QPushButton[segment="true"] {
    min-height: 27px;
    background: transparent;
    border: none;
    border-radius: 7px;
    padding: 1px 13px;
    color: #697586;
}
QPushButton[segment="true"]:hover {
    color: #24313C;
    background: #E7ECEF;
}
QPushButton[segment="true"]:checked {
    background: #FFFFFF;
    color: #16212B;
    font-weight: 600;
}
QFrame#SegmentedControl {
    background: #EDF1F4;
    border: 1px solid #D7DEE4;
    border-radius: 9px;
}
QGroupBox {
    margin-top: 15px;
    padding: 14px 10px 9px 10px;
    border: 1px solid #E2E7EB;
    border-radius: 10px;
    background: #FAFBFC;
    color: #26343E;
    font-weight: 600;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 6px;
    color: #5F6C78;
    background: #FFFFFF;
    font-size: 11px;
}
QScrollArea {
    background: transparent;
    border: none;
}
QScrollBar:vertical {
    background: transparent;
    width: 10px;
    margin: 2px;
}
QScrollBar::handle:vertical {
    background: #C8D0D7;
    border-radius: 4px;
    min-height: 28px;
}
QScrollBar::handle:vertical:hover {
    background: #AEB9C2;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0;
}
QScrollBar:horizontal {
    background: transparent;
    height: 10px;
    margin: 2px;
}
QScrollBar::handle:horizontal {
    background: #C8D0D7;
    border-radius: 4px;
    min-width: 28px;
}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
    width: 0;
}
QProgressBar {
    min-height: 6px;
    max-height: 6px;
    background: #E3E8EC;
    border: none;
    border-radius: 3px;
    text-align: center;
    color: transparent;
}
QProgressBar::chunk {
    background: #197FA3;
    border-radius: 3px;
}
QTextEdit {
    background: #F7F9FA;
    color: #45535F;
    border: 1px solid #E0E5E9;
    border-radius: 10px;
    padding: 9px;
    font-size: 11px;
    selection-background-color: #DDEEF3;
}
QSplitter::handle {
    background: transparent;
    width: 10px;
}
QToolTip {
    color: #FFFFFF;
    background: #26333D;
    border: 1px solid #26333D;
    border-radius: 6px;
    padding: 6px;
}
QMessageBox {
    background: #FFFFFF;
}
"""




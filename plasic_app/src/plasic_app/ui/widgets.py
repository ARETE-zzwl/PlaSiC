# 高屋建瓴地来说
# 这个脚本做了一件很简单的事：
# 它定义了一个“多选一”的自定义 UI 控件。
# 比如界面上有：
# [ 小时 ] [ 天 ] [ 月 ]
# 用户只能选一个。你点“天”，程序就把内部值设成 "day"，同时通知外部：“用户现在选的是 day”。

from __future__ import annotations

# Signal 是 Qt 的信号机制。
# 这里用于定义 SegmentedControl 对外暴露的 valueChanged 信号。
from PySide6.QtCore import Signal

from PySide6.QtWidgets import (
    # QButtonGroup 用于管理一组按钮。
    # 本组件通过它实现按钮之间的“互斥选择”，
    # 即同一时间只允许一个 segment 处于选中状态。
    QButtonGroup,
    # QFrame 作为 SegmentedControl 的基类。
    # 它本身是 QWidget 的子类，同时方便通过 QSS 对整个控件设置边框、背景等样式。
    QFrame,
    # 水平布局。
    # 所有 segment 按钮会从左到右依次排列。
    QHBoxLayout,
    # 每一个 segment 实际上使用一个 QPushButton 实现。
    QPushButton,
)


class SegmentedControl(QFrame):
    # 自定义 Qt 信号。
    #
    # 信号携带一个 str 参数，表示当前被选中的 segment 对应的 value。
    #
    # 例如：
    # options = [("小时", "hour"), ("天", "day")]
    #
    # 当用户点击“天”时，会 emit：
    #     valueChanged.emit("day")
    valueChanged = Signal(str)

    def __init__(
        self,
        # options 中每一项都是一个 (label, value) 二元组：
        options: list[tuple[str, str]],
        # Qt 标准 parent 参数。
        # 传入父 QObject/QWidget 后，Qt 会负责父子对象之间的生命周期管理。
        parent=None,
    ) -> None:
        # 初始化 QFrame 基类。
        super().__init__(parent)

        # SegmentedControl 至少必须存在一个选项。
        # 如果 options 是空列表：
        #     []
        # 那么组件不存在任何可选值，因此直接抛出 ValueError，
        # 避免后续 value() 等方法进入没有合法状态的情况。
        if not options:
            raise ValueError("segmented control requires options")

        # 设置 Qt objectName。
        #
        # 主要用途之一是通过 Qt Style Sheet（QSS）定位这个组件，例如：
        #
        #     QFrame#SegmentedControl {
        #         ...
        #     }
        #
        # objectName 不影响 Python 属性访问。
        self.setObjectName("SegmentedControl")

        # 创建按钮组，并将当前 SegmentedControl 作为 parent。
        #
        # QButtonGroup 本身不会负责按钮布局，
        # 它只负责按钮之间的逻辑分组和选择关系。
        self._group = QButtonGroup(self)

        # 设置为互斥模式。
        # 因为下面所有 QPushButton 都会设置为 checkable，
        # 所以 exclusive=True 后，一组按钮中只能有一个按钮处于 checked 状态。
        # 这正是 segmented control / 单选控件所需要的行为。
        self._group.setExclusive(True)

        # 保存：
        #     value -> QPushButton
        # 的映射关系。
        #
        # 例如：
        #     self._buttons = {
        #         "hour": <QPushButton ...>,
        #         "day": <QPushButton ...>,
        #     }
        #
        # 这样 set_value("day") 可以直接根据 value 找到对应按钮。
        self._buttons: dict[str, QPushButton] = {}

        # 创建水平布局，并直接绑定到当前 QFrame。
        # 由于传入了 self：
        #     QHBoxLayout(self)
        # 因此该 layout 会成为 SegmentedControl 的布局。
        layout = QHBoxLayout(self)

        # 设置布局四个方向的边距：
        #
        #     left   = 2
        #     top    = 2
        #     right  = 2
        #     bottom = 2
        #
        # 给外围 QFrame 与内部按钮之间留出少量空间。
        layout.setContentsMargins(2, 2, 2, 2)

        # 设置相邻按钮之间的间距为 1 像素。
        #
        # Segmented Control 通常要求按钮之间非常紧凑，
        # 因而这里使用较小的 spacing。
        layout.setSpacing(1)

        # 遍历所有选项，并为每一个选项创建一个按钮。
        #
        # enumerate() 同时提供：
        #     index：当前选项序号
        #     label：按钮显示文字
        #     value：按钮对应的内部值
        for index, (label, value) in enumerate(options):

            # 根据 label 创建按钮。
            #
            # 例如：
            #     label == "小时"
            #
            # 则按钮上显示“小时”。
            button = QPushButton(label)

            # 将普通 QPushButton 设置为可选中按钮。
            #
            # 默认 QPushButton 点击之后不会维持 checked 状态，
            # setCheckable(True) 后，它就可以拥有：
            #
            #     checked / unchecked
            #
            # 两种状态。
            button.setCheckable(True)

            # 给按钮设置一个动态 Qt Property：
            #
            #     segment = True
            #
            # 这个属性通常用于 QSS 样式选择，例如：
            #
            #     QPushButton[segment="true"] {
            #         ...
            #     }
            #
            # 因此可以只给 SegmentedControl 内的这些按钮应用特定样式。
            button.setProperty("segment", True)

            # 将按钮加入 QButtonGroup。
            #
            # 因为前面已经设置：
            #
            #     self._group.setExclusive(True)
            #
            # 所以加入同一个 group 的按钮之间会自动保持互斥。
            self._group.addButton(button)

            # 将当前 value 与对应的 QPushButton 建立映射。
            # 后续：
            #     value()
            #     set_value()
            # 都会依赖这个字典。
            self._buttons[value] = button

            # 将按钮加入水平布局。
            # 按照 options 的顺序从左到右排列。
            layout.addWidget(button)

            # 默认选中第一个选项。
            # index == 0 只对第一个按钮成立。
            # 因此 SegmentedControl 创建完成以后一定有一个初始选中值。
            if index == 0:
                button.setChecked(True)

            # 监听按钮的 clicked 信号。
            # QPushButton.clicked 信号本身会传递一个 checked 参数，
            # 所以 lambda 的第一个参数 _checked 用来接收它。
            #
            # 当前组件并不直接使用 _checked，
            # 因而以下划线开头表示“有这个参数，但这里不使用”。
            #
            # selected=value 非常关键：
            #
            # 它通过 lambda 的默认参数，在当前循环迭代时就把 value
            # 保存下来，从而避免 Python 闭包的 late binding（延迟绑定）问题。
            #
            # 如果简单写成：
            #
            #     lambda _checked: self.valueChanged.emit(value)
            #
            # 那么 lambda 中的 value 会在真正点击按钮时才查找，
            # 此时循环通常已经结束，所有 lambda 很可能都会引用最后一个 value。
            #
            # 当前写法：
            #
            #     selected=value
            #
            # 会让每个按钮分别保存自己对应的 value。
            button.clicked.connect(
                lambda _checked, selected=value: self.valueChanged.emit(selected)
            )

    def value(self) -> str:
        # 遍历所有 value -> button 映射，
        # 查找当前处于 checked 状态的按钮。
        for value, button in self._buttons.items():

            # 如果当前按钮被选中，直接返回其对应的内部 value。
            if button.isChecked():
                return value

        # 正常情况下，由于：
        #
        # 1. 至少存在一个 option；
        # 2. 第一个按钮初始化时会被 checked；
        # 3. QButtonGroup 设置为 exclusive；
        #
        # 理论上总能在上面的循环中找到一个 checked button。
        #
        # 这里仍然提供 fallback：
        # 如果由于某些外部状态修改导致没有按钮被选中，
        # 则返回字典中的第一个 value。
        #
        # next(iter(self._buttons))
        # 等价于取得字典迭代顺序中的第一个 key。
        return next(iter(self._buttons))

    def set_value(self, value: str) -> None:
        # 根据传入的 value 找到对应按钮。
        #
        # 这里直接使用字典索引，因此：
        #
        # 如果传入不存在的 value，例如：
        #
        #     self.set_value("invalid")
        #
        # 会直接抛出 KeyError。
        #
        # 当前实现没有主动捕获或转换这个异常。
        button = self._buttons[value]

        # 只有当目标按钮当前还没有被选中时，
        # 才真正执行状态切换和信号发送。
        #
        # 这样可以避免重复调用：
        #
        #     set_value("day")
        #     set_value("day")
        #
        # 时连续发送两次相同的 valueChanged 信号。
        if not button.isChecked():

            # 将目标按钮设置为 checked。
            # 因为 QButtonGroup 是 exclusive，
            # Qt 会自动取消当前其他按钮的 checked 状态。
            button.setChecked(True)

            # 主动发出 valueChanged 信号。
            #
            # 这里需要手动 emit，是因为程序调用：
            #
            #     button.setChecked(True)
            #
            # 与用户实际点击按钮是两种不同操作。
            #
            # 上面连接的是 button.clicked 信号，
            # setChecked(True) 本身不会模拟一次 clicked。
            #
            # 因此这里手动发送 valueChanged，
            # 从而保证“代码修改值”和“用户点击修改值”
            # 都能够通知 SegmentedControl 的外部使用者。
            self.valueChanged.emit(value)

    def set_option_enabled(self, value: str, enabled: bool) -> None:
        """Enable or disable one segment while keeping a valid selection."""
        button = self._buttons[value]
        button.setEnabled(enabled)
        if not enabled and button.isChecked():
            replacement = next(
                (candidate for candidate in self._buttons.values() if candidate.isEnabled()),
                None,
            )
            if replacement is not None:
                replacement.click()

    def option_enabled(self, value: str) -> bool:
        return self._buttons[value].isEnabled()

    def set_option_tooltip(self, value: str, text: str) -> None:
        self._buttons[value].setToolTip(text)

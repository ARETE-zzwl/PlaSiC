#                    ParameterCatalog
#                    参数 Schema / Metadata
#                           │
#                           │ exposed()
#                           ▼
#               ┌──────────────────────┐
#               │    ParameterForm     │
#               │                      │
#               │  Dynamic UI Builder  │
#               └──────────┬───────────┘
#                          │
#             ┌────────────┼────────────┐
#             │            │            │
#             ▼            ▼            ▼
#         boolean        number        path
#             │            │            │
#             ▼            ▼            ▼
#        QCheckBox   QDoubleSpinBox  PathEditor
#                                       │
#                                   ┌───┴────┐
#                                   ▼        ▼
#                               QLineEdit  QPushButton

# 参数变化：
# Widget
#   │
#   │ Qt Signal
#   ▼
# _emit_value(name)
#   │
#   ├── update dependencies
#   │
#   └── ParameterForm.valueChanged(name, value)
#                   │
#                   ▼
#              外部应用程序

# 而配置数据：
# GUI widgets
#     │
#     │ values()
#     ▼
# Python dict

# 反方向：
# Python dict
#     │
#     │ set_values()
#     ▼
# GUI widgets


from __future__ import annotations

from collections import defaultdict

from pathlib import Path

from typing import Any

from PySide6.QtCore import Signal

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

# 参数目录/参数元数据对象。
# ParameterForm 会根据这个 catalog 动态构建整个参数配置界面。
from plasic_app.config.catalog import ParameterCatalog


class PathEditor(QWidget):
    """
    文件路径 / 目录路径编辑控件。
    这个控件把两个功能组合在一起：
    1. 一个 QLineEdit：
       - 显示当前路径；
       - 允许用户直接手动输入路径。
    2. 一个 QPushButton：
       - 点击后打开文件或目录选择对话框；
       - 用户选中路径后自动回填到 QLineEdit。
    参数
    ----
    value:
        初始路径字符串。
    directory:
        True  表示选择目录；
        False 表示选择文件。
    parent:
        Qt 父控件。
    """

    # 自定义信号。
    # Signal(str) 表示该信号携带一个字符串参数。
    # 当 QLineEdit 内文本发生变化时，会将新的路径字符串继续向外发出。
    valueChanged = Signal(str)

    def __init__(
        self, value: str, directory: bool, parent: QWidget | None = None
    ) -> None:
        # 初始化 QWidget 基类，并建立 Qt 父子对象关系。
        super().__init__(parent)

        # 保存当前 PathEditor 的模式：
        # True  -> 浏览目录
        # False -> 浏览文件
        self.directory = directory

        # 路径文本输入框，初始内容为传入的 value。
        self.edit = QLineEdit(value)

        # 浏览按钮。
        # 使用省略号字符 “…” 表达“打开选择对话框”这一常见 UI 含义。
        self.button = QPushButton("…")

        # 限制按钮最大宽度，避免按钮占据过多横向空间。
        self.button.setMaximumWidth(34)

        # PathEditor 内部采用水平布局：
        # +--------------------------------------+
        # | [ 路径输入框................ ] [ … ] |
        # +--------------------------------------+
        layout = QHBoxLayout(self)

        # 由于 PathEditor 本身通常被作为表单里的一个 field 使用，
        # 内部布局不额外留下边距，以便与其他普通输入控件视觉对齐。
        layout.setContentsMargins(0, 0, 0, 0)

        # 输入框和浏览按钮之间留 6 px 间隔。
        layout.setSpacing(6)

        # 先添加路径输入框。
        layout.addWidget(self.edit)

        # 再添加浏览按钮。
        layout.addWidget(self.button)

        # QLineEdit.textChanged 信号本身也携带 str。
        # 这里直接把它连接到 PathEditor.valueChanged，
        # 相当于把内部输入框的文本变化“转发”为整个自定义控件的值变化。
        self.edit.textChanged.connect(self.valueChanged)

        # 点击浏览按钮时调用 _browse，
        # 打开对应的文件 / 目录选择对话框。
        self.button.clicked.connect(self._browse)

    def _browse(self) -> None:
        """
        打开路径选择对话框。
        directory=True:
            使用 QFileDialog.getExistingDirectory() 选择目录。
        directory=False:
            使用 QFileDialog.getOpenFileName() 选择文件。
        如果用户取消对话框，Qt 返回空字符串；
        此时不会修改当前输入框内容。
        """

        if self.directory:
            # 目录模式。
            # 参数：
            # self               -> 对话框父窗口
            # "Select Directory" -> 对话框标题
            # self.edit.text()   -> 初始打开位置
            value = QFileDialog.getExistingDirectory(
                self, "Select Directory", self.edit.text()
            )
        else:
            # 文件模式。
            # getOpenFileName 返回一个二元组：
            # (选中文件路径, 使用的文件过滤器)
            # 当前这里只关心路径，因此第二个返回值用 _ 丢弃。
            value, _ = QFileDialog.getOpenFileName(
                self, "Select File", self.edit.text()
            )

        # 用户真正选中了文件/目录时 value 非空。
        # setText 会触发 edit.textChanged，
        # 进而触发 PathEditor.valueChanged。
        if value:
            self.edit.setText(value)

    def value(self) -> str:
        """
        获取当前路径值。
        统一提供 value() 接口，使 ParameterForm 不需要直接访问内部 QLineEdit。
        """
        return self.edit.text()

    def set_value(self, value: str) -> None:
        """
        设置当前路径值。
        注意：
        QLineEdit.setText() 在文本实际发生变化时会触发 textChanged，
        因此这里也可能间接触发 valueChanged。
        """
        self.edit.setText(value)


class ParameterForm(QWidget):
    """
    根据 ParameterCatalog 动态生成参数配置表单。

    主要职责包括：

    1. 从 catalog.exposed() 获取需要暴露给用户的参数；
    2. 根据参数 type 自动创建合适的 Qt 控件；
    3. 按照 section 对参数进行分组；
    4. 支持参数搜索；
    5. 支持隐藏 / 显示 advanced 参数；
    6. 支持 depends_on 参数依赖关系；
    7. 对外统一提供读取和设置参数值的接口；
    8. 参数发生变化时发出 valueChanged(name, value) 信号。
    """

    # 参数变化信号。
    # 第一个参数 str：
    #     参数名，例如 "input_path"
    # 第二个参数 object：
    #     参数当前值，可以是 bool / int / float / str / enum data 等。
    valueChanged = Signal(str, object)

    def __init__(
        self,
        catalog: ParameterCatalog,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)

        # 设置 objectName。
        self.setObjectName("ParameterForm")

        # 保存参数目录对象。
        self.catalog = catalog

        # 保存：
        # 参数名 -> 对应 Qt 控件
        # 例如：
        # {
        #     "enabled": QCheckBox(...),
        #     "iterations": QSpinBox(...),
        #     "output_path": PathEditor(...),
        # }
        # 这样 value() / set_values() 可以通过参数名直接查找控件。
        self.controls: dict[str, QWidget] = {}

        # 保存所有表单行的元数据。
        # 每一项都是：
        # (
        #     definition,  # catalog 中的参数定义
        #     label,       # 左侧 QLabel
        #     control,     # 右侧参数输入控件
        # )
        #
        # 后续过滤、advanced 显示和依赖更新都会遍历该列表。
        self.rows: list[tuple[dict[str, Any], QLabel, QWidget]] = []

        # -------------------------------------------------------------
        # 整个 ParameterForm 的最外层布局
        # -------------------------------------------------------------
        outer = QVBoxLayout(self)

        # 外层不留额外边距，
        # 通常由 ParameterForm 的父布局决定整体间距。
        outer.setContentsMargins(0, 0, 0, 0)

        # 顶部过滤区与下面滚动区之间的垂直间距。
        outer.setSpacing(9)

        # -------------------------------------------------------------
        # 顶部过滤栏
        # [ Filter parameters...             ] [ Show advanced ]
        # -------------------------------------------------------------

        filter_row = QHBoxLayout()

        # 搜索框和 advanced checkbox 之间的间距。
        filter_row.setSpacing(8)

        # 参数搜索输入框。
        self.search = QLineEdit()

        # 搜索框 placeholder。
        self.search.setPlaceholderText("Filter parameters…")

        # 是否显示 advanced 参数。
        self.advanced = QCheckBox("Show advanced")

        # 搜索框 stretch=1：
        # 当横向空间增加时，搜索框优先扩展；
        # checkbox 一般保持自身 sizeHint 的宽度。
        filter_row.addWidget(self.search, 1)

        # 添加 advanced 开关。
        filter_row.addWidget(self.advanced)

        # 将整个过滤栏加入最外层布局。
        outer.addLayout(filter_row)

        # -------------------------------------------------------------
        # 可滚动参数区域
        # -------------------------------------------------------------

        scroll = QScrollArea()

        # 设置 objectName，方便 QSS / 自动测试等使用。
        scroll.setObjectName("ParameterScroll")

        # QScrollArea 的实际显示窗口是 viewport。
        # 给 viewport 单独设置 objectName，通常是为了更精确地设置背景等样式。
        scroll.viewport().setObjectName("ParameterViewport")

        # 当 ScrollArea 大小时，让内部 widget 自动调整宽度。
        # 如果为 False，内部 widget 会保留自己的固定尺寸，
        # 容易造成不必要的水平滚动。
        scroll.setWidgetResizable(True)

        # scroll 内部真正承载所有参数组的 content widget。
        content = QWidget()
        content.setObjectName("ParameterContent")

        # content 中使用垂直布局：
        # [Section A]
        # [Section B]
        # [Section C]
        # ...
        self.content_layout = QVBoxLayout(content)

        # content 自身不增加边距。
        self.content_layout.setContentsMargins(0, 0, 0, 0)

        # 各 section group box 之间的距离。
        self.content_layout.setSpacing(6)

        # 指定 QScrollArea 中真正需要滚动的 widget。
        scroll.setWidget(content)

        # 将 scroll 加入最外层布局。
        # stretch=1 表示它会占据 ParameterForm 中绝大多数剩余垂直空间。
        outer.addWidget(scroll, 1)

        # -------------------------------------------------------------
        # 按 section 对参数定义分组
        # -------------------------------------------------------------
        # sections 的结构：
        # {
        #     "General": [definition1, definition2, ...],
        #     "Input":   [definition3, definition4, ...],
        #     ...
        # }
        sections: dict[str, list[dict[str, Any]]] = defaultdict(list)

        # catalog.exposed() 应返回所有需要展示在 GUI 中的参数定义。
        # 每个 definition 预计至少包含：
        # {
        #     "name": ...,
        #     "label": ...,
        #     "type": ...,
        #     "section": ...,
        # }
        for definition in catalog.exposed():
            # 根据 definition["section"] 分类。
            sections[definition["section"]].append(definition)

        # 保存所有 section 的 QGroupBox。
        # 后面 _update_visibility() 会控制整个 section 是否显示。
        self.groups: list[QGroupBox] = []

        # -------------------------------------------------------------
        # 为每一个 section 创建 QGroupBox + QFormLayout
        # -------------------------------------------------------------

        for section, definitions in sections.items():
            # section 名作为 GroupBox 标题。
            group = QGroupBox(section)

            # 每个 section 内使用 QFormLayout。
            # 典型布局：
            # Parameter A: [ control ]
            # Parameter B: [ control ]
            # Parameter C: [ control ]
            form = QFormLayout(group)

            # form 内边距：
            # left=10, top=13, right=10, bottom=10
            # top 稍大一般是为了给 GroupBox 标题留视觉空间。
            form.setContentsMargins(10, 13, 10, 10)

            # Label 列与 Field 列之间的横向距离。
            form.setHorizontalSpacing(12)

            # 各参数行之间的垂直距离。
            form.setVerticalSpacing(7)

            # 指定 field 列的增长策略。
            # AllNonFixedFieldsGrow 表示所有非固定尺寸 field
            # 可以使用表单中的剩余水平空间。
            form.setFieldGrowthPolicy(
                QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
            )

            # 当窗口变窄、某一行长度不足时允许换行。
            # WrapLongRows：
            # 对过长的行进行 wrap，而不是强行挤压。
            form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)

            # 为 section 中每一个参数创建一行。
            for definition in definitions:
                # 根据参数类型创建实际输入控件。
                control = self._create_control(definition)

                # 左侧参数标签。
                label = QLabel(definition["label"])

                # 参数说明。
                # get(..., "") 可保证字段缺失时仍返回字符串，
                # 避免 KeyError。
                description = definition.get("description", "")

                # 参数警告信息。
                warning = definition.get("warning", "")

                # 将 description 和 warning 合并为 tooltip。
                # 只加入非空字符串，避免出现多余空行。
                tooltip = "\n".join(value for value in (description, warning) if value)

                # Label 和 Control 都设置 tooltip。
                # 这样无论鼠标停留在参数名还是输入控件上，
                # 用户都可以看到同样的说明。
                label.setToolTip(tooltip)
                control.setToolTip(tooltip)

                # 将当前参数作为一行加入表单。
                form.addRow(label, control)

                # 记录：
                # 参数名 -> 控件
                self.controls[definition["name"]] = control

                # 同时记录该行的 definition / label / control。
                # 后续过滤显示和依赖逻辑使用。
                self.rows.append((definition, label, control))

            # 将当前 section 加入滚动区域的垂直布局。
            self.content_layout.addWidget(group)

            # 保存 group 引用。
            self.groups.append(group)

        # 在所有 group 后添加 stretch。
        # 如果参数组总高度小于滚动区域，
        # stretch 会把 group 顶到顶部，而不是在整个区域里垂直分散。
        self.content_layout.addStretch(1)

        # -------------------------------------------------------------
        # 信号连接
        # -------------------------------------------------------------

        # 搜索关键字发生改变时，重新计算各参数是否可见。
        self.search.textChanged.connect(self._update_visibility)

        # Show advanced 状态变化时，也重新计算参数可见性。
        self.advanced.toggled.connect(self._update_visibility)

        # 初始化一次 UI 可见性。
        # 这样初始状态下 advanced 参数会按 checkbox 状态正确隐藏，
        # 同时依赖关系也会被更新。
        self._update_visibility()

    def _emit_value(self, name: str) -> None:
        """
        当某一个参数控件变化时执行。

        执行顺序：

        1. 更新所有依赖控件的 enabled 状态；
        2. 读取发生变化参数的当前值；
        3. 对外发出 valueChanged(name, value)。

        这里先更新 dependencies，
        可以保证外部监听 valueChanged 时看到的界面状态已经同步。
        """

        # 当前参数改变可能会影响其他参数的 depends_on 条件，
        # 因此先重新计算依赖关系。
        self._update_dependencies()

        # 发出：
        # 参数名 + 参数当前值
        self.valueChanged.emit(name, self.value(name))

    def _create_control(self, definition: dict[str, Any]) -> QWidget:
        """
        根据参数 definition["type"] 创建对应 Qt 控件。

        当前类型映射：
        boolean   -> QCheckBox
        enum      -> QComboBox
        integer   -> QSpinBox
        number    -> QDoubleSpinBox
        path      -> PathEditor(directory=False)
        directory -> PathEditor(directory=True)
        其他类型   -> QLineEdit

        所有控件的“值发生变化”信号最终都会连接到：
            self._emit_value(parameter_name)

        从而让 ParameterForm 对外统一发出 valueChanged。
        """

        # 参数数据类型。
        kind = definition["type"]

        # 参数默认值。
        default = definition.get("default")

        # 参数内部名称。
        name = definition["name"]

        # -------------------------------------------------------------
        # boolean -> checkbox
        # -------------------------------------------------------------
        if kind == "boolean":
            control = QCheckBox()

            # 将默认值转换为 bool。
            control.setChecked(bool(default))

            # toggled 信号会传出一个 bool。
            # 这里不直接使用信号中的值，
            # 而是统一由 _emit_value -> self.value(name) 获取当前值。
            # n=name 是非常关键的：
            # 它把当前循环中的 name 作为 lambda 默认参数立即绑定，
            # 避免 Python closure 的 late binding 问题。
            control.toggled.connect(lambda _value, n=name: self._emit_value(n))

            return control

        # -------------------------------------------------------------
        # enum -> combobox
        # -------------------------------------------------------------
        if kind == "enum":
            control = QComboBox()

            # choice_labels 可用于：
            # 实际值:
            #     "rk4"
            # GUI 显示:
            #     "Runge-Kutta 4"
            # 若没有映射，则直接显示 choice 自身的字符串形式。
            labels = definition.get("choice_labels", {})

            # 将所有枚举选项加入 QComboBox。
            for choice in definition["choices"]:
                # addItem(display_text, userData)
                # display_text:
                #     用户看到的字符串。
                # userData:
                #     真正对应的参数值。
                control.addItem(labels.get(str(choice), str(choice)), choice)

            # 根据 userData 查找 default 对应的下拉项。
            index = control.findData(default)

            # 如果没有找到，findData 返回 -1。
            # max(index, 0) 保证至少选择第一个选项。
            control.setCurrentIndex(max(index, 0))

            # currentIndexChanged 信号携带新的 index。
            # 这里同样忽略 index，
            # 统一使用 self.value(name) 读取 currentData。
            control.currentIndexChanged.connect(
                lambda _value, n=name: self._emit_value(n)
            )

            return control

        # -------------------------------------------------------------
        # integer -> QSpinBox
        # -------------------------------------------------------------
        if kind == "integer":
            control = QSpinBox()
            # 设置整数允许范围。
            # definition 未指定 minimum / maximum 时，
            # 使用接近 32-bit signed integer 范围的默认值。
            control.setRange(
                int(definition.get("minimum", -2147483647)),
                int(definition.get("maximum", 2147483647)),
            )

            # 设置默认值。
            control.setValue(int(default))

            # 用户修改整数后通知 ParameterForm。
            control.valueChanged.connect(lambda _value, n=name: self._emit_value(n))

            return control

        # -------------------------------------------------------------
        # number -> QDoubleSpinBox
        # -------------------------------------------------------------
        if kind == "number":
            control = QDoubleSpinBox()

            # 最多显示 / 编辑 12 位小数。
            # 注意这里只是 QDoubleSpinBox 的显示/编辑精度配置；
            # 底层返回值仍然是 Python float。
            control.setDecimals(12)

            # 设置浮点数范围。
            # 如果 definition 没有定义范围，
            # 使用非常宽的 [-1e30, 1e30]。
            control.setRange(
                float(definition.get("minimum", -1.0e30)),
                float(definition.get("maximum", 1.0e30)),
            )

            # 设置默认值。
            control.setValue(float(default))

            # 数值变化时触发统一的 _emit_value。
            control.valueChanged.connect(lambda _value, n=name: self._emit_value(n))

            return control

        # -------------------------------------------------------------
        # path / directory -> PathEditor
        # -------------------------------------------------------------
        if kind in {"path", "directory"}:
            control = PathEditor(
                # 路径统一转换成字符串。
                str(default),
                # directory 类型打开目录选择器，
                # path 类型打开文件选择器。
                directory=kind == "directory",
            )

            # PathEditor.valueChanged 携带字符串，
            # 但这里同样忽略信号值，
            # 统一通过 _emit_value -> value() 获取。
            control.valueChanged.connect(lambda _value, n=name: self._emit_value(n))

            return control

        # -------------------------------------------------------------
        # fallback -> QLineEdit
        # -------------------------------------------------------------
        # 如果 type 不是当前显式支持的类型，
        # 默认使用普通字符串输入框。
        control = QLineEdit(str(default))

        # 文本变化时通知 ParameterForm。
        control.textChanged.connect(lambda _value, n=name: self._emit_value(n))

        return control

    def value(self, name: str) -> Any:
        """
        获取指定参数当前值。

        通过 Qt 控件的实际类型决定怎样读取值，
        从而在 ParameterForm 外部提供统一接口。

        Parameters
        ----------
        name:
            参数名称。

        Returns
        -------
        Any
            参数当前值。
        """

        # 根据参数名取得对应控件。
        # 这里使用 [] 而不是 get()：
        # 如果 name 不存在，会直接抛出 KeyError，
        # 可以尽早暴露调用者传错参数名的问题。
        control = self.controls[name]

        # checkbox -> bool
        if isinstance(control, QCheckBox):
            return control.isChecked()

        # combo box -> 当前选项 userData
        # 注意不是 currentText()，
        # 因为 UI 显示文本可能经过 choice_labels 转换。
        if isinstance(control, QComboBox):
            return control.currentData()

        # spin box -> int 或 float
        # QSpinBox.value() 返回 int；
        # QDoubleSpinBox.value() 返回 float。
        if isinstance(control, (QSpinBox, QDoubleSpinBox)):
            return control.value()

        # 自定义路径编辑器 -> str
        if isinstance(control, PathEditor):
            return control.value()

        # 普通文本输入 -> str
        if isinstance(control, QLineEdit):
            return control.text()

        # 如果未来 controls 中加入了新的控件类型，
        # 但这里没有增加相应读取逻辑，则显式抛出 TypeError，
        # 避免静默返回错误结果。
        raise TypeError(type(control))

    def values(self) -> dict[str, Any]:
        """
        获取当前所有参数值。

        返回形式：
        {
            parameter_name: current_value,
            ...
        }
        """

        # 遍历 controls 中的所有参数，
        # 对每个参数调用统一的 value(name)。
        return {name: self.value(name) for name in self.controls}

    def set_values(self, values: dict[str, Any]) -> None:
        """
        批量设置参数值。

        values 中允许只包含部分参数。

        如果 values 包含 catalog 当前不存在的参数名，
        会直接跳过，而不是报错。

        这对于：
        - 加载旧版本配置文件；
        - 部分更新；
        - 配置向前/向后兼容；
        通常比较方便。
        """

        # 遍历调用者传入的参数。
        for name, value in values.items():
            # 找到对应 Qt 控件。
            control = self.controls.get(name)

            # 如果参数不存在，忽略。
            if control is None:
                continue

            # ---------------------------------------------------------
            # 根据实际控件类型设置值
            # ---------------------------------------------------------

            if isinstance(control, QCheckBox):
                # checkbox 使用 bool。
                control.setChecked(bool(value))

            elif isinstance(control, QComboBox):
                # enum 控件内部保存的是 userData，
                # 因此通过 findData 查找实际参数值对应的选项。
                index = control.findData(value)

                # 只有找到有效选项才更新。
                # 如果传入一个 catalog 中不存在的 enum value，
                # 保持原选择不变。
                if index >= 0:
                    control.setCurrentIndex(index)

            elif isinstance(control, (QSpinBox, QDoubleSpinBox)):
                # QSpinBox / QDoubleSpinBox 都提供 setValue。
                control.setValue(value)

            elif isinstance(control, PathEditor):
                # PathEditor 内部统一存储字符串。
                control.set_value(str(value))

            elif isinstance(control, QLineEdit):
                # 文本输入统一转换为 str。
                control.setText(str(value))

        # 批量修改之后重新计算一次 depends_on。
        # 虽然各个 setXxx() 本身可能已经触发 valueChanged，
        # 最后再统一执行一次可以保证最终依赖状态与最终 values 一致。
        self._update_dependencies()

    def reset_defaults(self) -> None:
        """
        将所有参数恢复到 catalog 定义的默认值。
        """

        # catalog.defaults() 应返回：
        # {
        #     parameter_name: default_value,
        #     ...
        # }
        # 直接复用 set_values()，
        # 避免重复实现不同控件的赋值逻辑。
        self.set_values(self.catalog.defaults())

    def _update_dependencies(self) -> None:
        """
        根据每个参数 definition 中的 depends_on
        更新控件 enabled / disabled 状态。

        depends_on 预计形式类似：

        {
            "depends_on": {
                "use_feature": True,
                "mode": "advanced"
            }
        }

        表示：

        当前参数只有在：
            use_feature == True
        且：
            mode == "advanced"

        时才允许编辑。

        注意：
        当前逻辑只是 disable control，
        不会隐藏 label/control，也不会修改参数本身已有的值。
        """

        # 首先获取所有参数当前值。
        # 后面所有 depends_on 判断都基于这个快照。
        values = self.values()

        # 遍历所有参数行。
        for definition, _label, control in self.rows:
            # 获取依赖定义。
            # 如果没有 depends_on，则默认空 dict。
            dependencies = definition.get("depends_on", {})

            # 所有依赖条件都满足时才 enabled。
            #
            # 例如：
            #
            # dependencies = {
            #     "foo": True,
            #     "bar": "x",
            # }
            #
            # 会检查：
            #
            # values["foo"] == True
            # AND
            # values["bar"] == "x"
            #
            # all([]) 在 Python 中为 True，
            # 所以没有 dependency 的参数自然保持 enabled。
            enabled = all(
                values.get(key) == expected for key, expected in dependencies.items()
            )

            # 设置当前参数控件是否允许交互。
            control.setEnabled(enabled)

    def _update_visibility(self) -> None:
        """
        根据以下两个条件更新参数的显示状态：

        1. 搜索框中的关键字；
        2. "Show advanced" 是否勾选。

        搜索范围包括：

        - definition["name"]
        - definition["label"]
        - definition["description"]
        - definition["component"]

        搜索使用 casefold()，
        因此大小写不敏感。

        如果一个 group 中没有任何可见参数，
        整个 QGroupBox 会被隐藏。
        """

        # 获取搜索字符串：
        # strip():
        #     去掉首尾空格。
        # casefold():
        #     用于大小写不敏感比较，
        #     比 lower() 更适合通用 Unicode 文本。
        query = self.search.text().strip().casefold()

        # 当前是否显示 advanced 参数。
        show_advanced = self.advanced.isChecked()

        # 记录至少包含一个 visible 参数的 GroupBox。
        # 最后只显示这些 group。
        visible_groups: set[QGroupBox] = set()

        # 遍历所有参数行。
        for definition, label, control in self.rows:
            # ---------------------------------------------------------
            # 构建可搜索字符串
            # ---------------------------------------------------------
            # 用空格连接后统一 casefold。
            searchable = " ".join(
                [
                    definition["name"],
                    definition["label"],
                    definition.get("description", ""),
                    definition.get("component", ""),
                ]
            ).casefold()

            # ---------------------------------------------------------
            # 判断当前参数是否应该显示
            # ---------------------------------------------------------
            #
            # 条件一：
            #
            # show_advanced == True
            # OR
            # 当前参数不是 advanced
            #
            # 条件二：
            #
            # query 为空
            # OR
            # query 是 searchable 的子串
            #
            # 两部分必须同时成立。
            visible = (show_advanced or not definition.get("advanced", False)) and (
                not query or query in searchable
            )

            # Label 与对应 Control 保持一致的 visible 状态。
            label.setVisible(visible)
            control.setVisible(visible)

            # 如果当前参数可见，
            # 则记录它所在的 GroupBox。
            if visible:
                # QFormLayout 中 label 的 parentWidget
                # 通常就是当前 section 对应的 QGroupBox。
                parent = label.parentWidget()

                # 做类型检查后加入集合。
                if isinstance(parent, QGroupBox):
                    visible_groups.add(parent)

        # -------------------------------------------------------------
        # 更新整个 section 的可见性
        # -------------------------------------------------------------

        for group in self.groups:
            # 只有至少含一个 visible 参数的 group 才显示。
            group.setVisible(group in visible_groups)

        # 可见性变化后重新计算一次依赖状态。
        # 注意：
        # visible 和 enabled 是两个独立概念：
        # visible:
        #     控件是否显示。
        # enabled:
        #     控件是否允许用户编辑。
        self._update_dependencies()

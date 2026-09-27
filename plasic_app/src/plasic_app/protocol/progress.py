from __future__ import annotations

import json
from pathlib import Path
from typing import Any

# progress.txt 中一条进度记录必须包含的字段。
# 该文件不再使用 protocol/version 或 event 事件类型。
_PROGRESS_FIELDS = frozenset(
    {
        "completed_steps",
        "total_steps",
        "date",
        "elapsed_seconds",
        "steps_per_second",
    }
)


class ProgressTail:
    """
    用于“增量读取”一个 progress 文件的辅助类。

    典型使用场景：
    - 某个外部进程持续向文件末尾追加 JSON Lines 格式的数据；
    - 当前程序周期性调用 read_new()；
    - 每次只读取“上一次读取位置之后”新增的数据；
    - 如果最后一行尚未写完整，则暂存下来，等待下一次读取时继续拼接。

    该类内部维护两个关键状态：

    1. offset
       当前已经从文件中读取到的字节位置。
       下一次读取时会从该位置继续，而不是重新从文件开头读取。

    2. remainder
       上一次读取时遇到的“不完整文本行”。
       例如文件最后只写入了一半 JSON：
           {"completed_steps": 1928
       此时不会立刻解析，而是保存到 remainder，
       等下一次读取到后半部分后再拼接成完整的一行。
    """

    def __init__(self, path: Path) -> None:
        self.path = Path(path)

        self.offset = 0

        self.remainder = ""

    def reset(self) -> None:

        # 将文件读取位置恢复到文件开头。
        self.offset = 0

        # 丢弃之前保存的未完成文本行。
        self.remainder = ""

    def read_new(self) -> list[dict[str, Any]]:
        """
        读取自上一次调用以来文件新增的完整事件。

        返回值：
            list[dict[str, Any]]
            每一个元素对应文件中的一条完整 JSON 记录。

        文件中的每条记录必须包含以下字段：
            completed_steps、total_steps、date、
            elapsed_seconds、steps_per_second

        否则会抛出：
            ValueError("unsupported progress record")

        此外，如果某一完整行不是合法 JSON，
        json.loads() 会直接抛出 JSONDecodeError。
        """

        if not self.path.exists():
            return []

        with self.path.open("rb") as stream:
            stream.seek(self.offset)

            # 从当前 offset 一直读取到当前文件末尾。
            # payload 的类型为 bytes。
            payload = stream.read()

            # 保存读取结束后的文件位置，
            # 供下一次调用 read_new() 时继续读取。
            #
            # 即使当前新增内容中包含“不完整的最后一行”，
            # offset 仍然会推进到文件末尾。
            #
            # 那部分不完整数据不会通过重新读取文件获得，
            # 而是通过下面的 self.remainder 保存在内存中。
            self.offset = stream.tell()

        # 如果从上一次 offset 到当前文件末尾没有新增任何字节，
        # 则直接返回空列表。
        # remainder 此时保持原样：
        # 如果之前存在半行数据，它仍然等待未来新增内容进行拼接。
        if not payload:
            return []

        text = self.remainder + payload.decode("utf-8")

        # 按“行”拆分文本，同时保留每一行末尾的换行符。
        #
        # keepends=True 非常关键，因为后续需要通过检查
        # 最后一行是否以 '\n' 或 '\r' 结尾，
        # 判断该行是否已经完整写入。
        # 示例：
        #   "a\nb\n" -> ["a\n", "b\n"]
        #   "a\nb"   -> ["a\n", "b"]
        # 第二种情况下，"b" 没有换行符，
        # 因而会被视为可能仍在写入中的“不完整行”。
        lines = text.splitlines(keepends=True)

        # 如果存在至少一行，并且最后一行没有以换行符结尾，
        # 则认为最后一行可能尚未写完整。
        # 不在本次调用中解析它，而是将其从 lines 中取出，
        # 保存到 self.remainder，等待下一次调用补全。
        if lines and not lines[-1].endswith(("\n", "\r")):
            # pop() 同时完成两件事：
            # 1. 从本次待解析的 lines 中移除最后一行；
            # 2. 将该不完整行保存到 remainder。
            self.remainder = lines.pop()
        else:
            # 如果最后一行已经以换行符结束，
            # 说明当前读取到的文本全部由完整行组成，
            # 不需要保留任何 remainder。
            self.remainder = ""

        # 用于收集本次成功解析并通过协议校验的事件。
        events = []

        # 逐行处理所有“确认已经完整写入”的文本行。
        for line in lines:
            # 忽略空行以及只包含空白字符的行。
            #
            # strip() 会移除：
            # - 空格
            # - tab
            # - '\n'
            # - '\r'
            # 等空白字符。
            if line.strip():
                # 将当前完整文本行解析为 JSON。
                #
                # 这里预期 JSON 顶层结构是一个 object，
                # 因此后续代码直接调用 event.get(...)。
                #
                # 如果 line 不是合法 JSON，
                # json.loads() 会抛出 json.JSONDecodeError，
                # 本函数不会捕获该异常。
                event = json.loads(line)

                # 新的 progress.txt 记录不再携带 protocol/version，
                # 因此改为检查五个固定进度字段是否齐全。
                if (
                    not isinstance(event, dict)
                    or not _PROGRESS_FIELDS.issubset(event)
                ):
                    raise ValueError("unsupported progress record")

                events.append(event)

        return events
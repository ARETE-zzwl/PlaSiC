#!/usr/bin/env python3
"""Reader/writer for the PlaSiC sequential-unformatted restart format.

中文说明：PlaSiC 顺序无格式（sequential-unformatted）重启文件的读写工具。

Format (little-endian, native byte order on the machines used here):

中文说明：文件格式（小端序，本试验所用机器上为原生字节序）：

    [u32 16][16-byte blank-padded name][u32 16]
    [u32 N][N bytes payload][u32 N]

English: every record is delimited by two equal u32 markers, so a truncated
file can be detected immediately.  Payloads are stored column-major: for an
``R x C`` record the file holds C consecutive blocks of R float32 values,
i.e. ``array[row + column * R]``.  Integer records hold a single int32 value.

中文说明：每条记录由前后两个相等的 u32 标记界定，因此可以立即发现被截断
的文件。数据按列主序存储：对于 ``R x C`` 的记录，文件中是 C 个连续的
“R 个 float32”数据块，即 ``array[row + column * R]``。整型记录只包含一个
int32 数值。
"""
from __future__ import annotations

import struct
from pathlib import Path

import numpy as np

NAME_BYTES = 16    # 记录名固定为 16 字节 / record name is always 16 bytes
FRAME = "<I"       # 记录边界标记：小端 u32 / record frame marker: little-endian u32


def read_raw_records(path: str | Path) -> "list[tuple[str, bytes]]":
    """Return ``[(name, payload_bytes), ...]`` preserving file order.

    中文说明：按文件中的顺序返回 ``[(记录名, 数据字节), ...]`` 列表。
    """
    records: list[tuple[str, bytes]] = []
    with Path(path).open("rb") as stream:
        while True:
            lead = stream.read(4)
            if lead == b"":
                break  # 文件正常结束 / clean end of file
            if len(lead) != 4:
                raise ValueError("truncated name marker")
            (name_bytes,) = struct.unpack(FRAME, lead)
            if name_bytes != NAME_BYTES:
                raise ValueError(f"unexpected name length {name_bytes}")
            name = stream.read(NAME_BYTES)
            trail = stream.read(4)
            if len(name) != NAME_BYTES or len(trail) != 4:
                raise ValueError("truncated name record")
            (trail_bytes,) = struct.unpack(FRAME, trail)
            if trail_bytes != name_bytes:
                raise ValueError("mismatched name markers")

            lead = stream.read(4)
            if len(lead) != 4:
                raise ValueError("truncated data marker")
            (payload_bytes,) = struct.unpack(FRAME, lead)
            payload = stream.read(payload_bytes)
            trail = stream.read(4)
            if len(payload) != payload_bytes or len(trail) != 4:
                raise ValueError("truncated data record")
            (trail_bytes,) = struct.unpack(FRAME, trail)
            if trail_bytes != payload_bytes:
                raise ValueError("mismatched data markers")

            key = name.decode("ascii", "replace").strip()  # 去掉右侧补空格 / strip blank padding
            records.append((key, payload))
    return records


def decode_float(payload: bytes, count: int | None = None) -> np.ndarray:
    """Decode one float32 payload into a NumPy array.

    中文说明：把一段 float32 数据解码为 NumPy 数组。

    When ``count`` is given, the payload length is checked against it.
    中文说明：若给定了 ``count``，则同时校验数据个数是否与预期一致。
    """
    values = np.frombuffer(payload, dtype="<f4")
    if count is not None and values.size != count:
        raise ValueError(f"payload holds {values.size} floats, expected {count}")
    return values.copy()


def decode_int(payload: bytes) -> int:
    """Decode a 4-byte little-endian int32 record.

    中文说明：把 4 字节小端 int32 数据解码为整数。
    """
    if len(payload) != 4:
        raise ValueError("integer record must hold exactly 4 bytes")
    return struct.unpack("<i", payload)[0]


def encode_float(values: np.ndarray) -> bytes:
    """Encode an array as a contiguous little-endian float32 payload.

    中文说明：把数组编码为连续的 float32 小端数据块。
    """
    return np.ascontiguousarray(np.asarray(values, dtype="<f4").reshape(-1)).tobytes()


def encode_int(value: int) -> bytes:
    """Encode an integer as a 4-byte little-endian int32 payload.

    中文说明：把整数编码为 4 字节小端 int32 数据块。
    """
    return struct.pack("<i", int(value))


def write_raw_records(path: str | Path, records: "list[tuple[str, bytes]]") -> None:
    """Write records in the given order; payload bytes are passed through.

    中文说明：按给定顺序写出记录，数据字节原样写入、不做转换。
    """
    with Path(path).open("wb") as stream:
        for name, payload in records:
            # 记录名右补空格到 16 字节 / pad the record name with blanks to 16 bytes
            key = name.encode("ascii")[:NAME_BYTES].ljust(NAME_BYTES, b" ")
            stream.write(struct.pack(FRAME, NAME_BYTES) + key + struct.pack(FRAME, NAME_BYTES))
            stream.write(struct.pack(FRAME, len(payload)) + payload + struct.pack(FRAME, len(payload)))


def replace_records(
    template: "list[tuple[str, bytes]]", replacements: dict[str, bytes]
) -> "list[tuple[str, bytes]]":
    """Return a copy of template with the given records replaced.

    中文说明：返回模板记录列表的副本，其中指定记录被替换为新数据。

    Unknown names in ``replacements`` are appended at the end.
    中文说明：``replacements`` 中模板里不存在的记录名会被追加到列表末尾。
    """
    seen = set()
    output: list[tuple[str, bytes]] = []
    for name, payload in template:
        if name in replacements:
            output.append((name, replacements[name]))
            seen.add(name)
        else:
            output.append((name, payload))
    for name, payload in replacements.items():
        if name not in seen:
            output.append((name, payload))
    return output

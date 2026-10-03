#!/usr/bin/env python3
"""Reader for the PlaSiC ``--frames`` binary stream.

中文说明：PlaSiC ``--frames`` 二进制帧流的读取工具。

File header (64 bytes)::

中文说明：文件头（64 字节）::

    "PLASICPS" | u16 version | u16 header_bytes | u32 0 | u32 0x01020304
    u32 nlon | u32 nlat | u32 nlev | u32 ocean_levels | u32 sizeof(float)

Frame header (128 bytes)::

中文说明：帧头（128 字节）::

    "PFRAME1\\0" | u16 version | u16 header_bytes | u32 variable_id
    u32 payload_bytes | u64 step | u32 year, month, day, hour, minute
    u32 level_index | u32 nlon | u32 nlat | u32 layers
    f32 minimum | f32 maximum | u32 payload_hash | char name[32] | char unit[16]

English: the payload is float32 ``[layer][lat][lon]`` with latitude running
north -> south and longitude fastest, exactly the model's global Gaussian grid.

中文说明：数据体为 float32 的 ``[层][纬][经]`` 数组，纬度由北向南排列，
经度方向变化最快，与模式的全球高斯网格完全一致。
"""
from __future__ import annotations

import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

FILE_HEADER_BYTES = 64   # 文件头固定长度 / fixed size of the file header
FRAME_HEADER_BYTES = 128 # 每帧帧头固定长度 / fixed size of each frame header


@dataclass
class Frame:
    """One decoded frame of the stream.

    中文说明：帧流中解码出来的单帧数据。
    """

    variable_id: int        # 变量编号 / variable identifier
    name: str               # 变量名 / variable name
    unit: str               # 单位 / physical unit
    step: int               # 模式时间步 / model time step
    date: tuple[int, int, int, int, int]  # (年, 月, 日, 时, 分) / (year, month, day, hour, minute)
    level_index: int        # 层序号 / level index
    layers: int             # 数据层层数 / number of vertical layers
    nlon: int               # 经向格点数 / number of longitudes
    nlat: int               # 纬向格点数 / number of latitudes
    values: np.ndarray      # 形状 (layers, nlat, nlon) / shape (layers, nlat, nlon)

    @property
    def valid_time(self):
        """Frame valid time as a ``datetime.datetime``.

        中文说明：把帧的日期字段转换为 ``datetime.datetime`` 有效时刻。
        """
        import datetime as dt

        year, month, day, hour, minute = self.date
        return dt.datetime(year, month, day, hour, minute)


def iter_frames(path: str | Path):
    """Yield :class:`Frame` objects from a frame stream file.

    中文说明：逐个生成帧流文件中的 :class:`Frame` 对象。
    """
    with Path(path).open("rb") as stream:
        header = stream.read(FILE_HEADER_BYTES)
        if len(header) != FILE_HEADER_BYTES or header[:8] != b"PLASICPS":
            raise ValueError(f"{path} is not a PlaSiC frame stream")
        # 文件头中记录的网格尺寸（此处读出后仅作格式校验）
        # grid dimensions recorded in the file header (read here for format checking)
        nlon, nlat, nlev = struct.unpack_from("<III", header, 20)
        while True:
            frame_header = stream.read(FRAME_HEADER_BYTES)
            if frame_header == b"":
                break  # 正常读到文件尾 / clean end of stream
            if len(frame_header) != FRAME_HEADER_BYTES or frame_header[:8] != b"PFRAME1\0":
                raise ValueError("corrupt frame header")
            # 帧头各字段的偏移量必须与模式 C 代码写出的布局保持一致。
            # Field offsets must match the layout written by the model's C code.
            (variable_id, payload_bytes) = struct.unpack_from("<II", frame_header, 12)
            (step,) = struct.unpack_from("<Q", frame_header, 20)
            date = struct.unpack_from("<IIIII", frame_header, 28)
            (level_index, frame_nlon, frame_nlat, layers) = struct.unpack_from("<IIII", frame_header, 48)
            name = frame_header[76:108].split(b"\0")[0].decode("ascii", "replace")
            unit = frame_header[108:124].split(b"\0")[0].decode("ascii", "replace")
            payload = stream.read(payload_bytes)
            if len(payload) != payload_bytes:
                raise ValueError("truncated frame payload")
            # 按 C 顺序的 [层][纬][经] 直接 reshape，无需转置。
            # Reshape directly in C order [layer][lat][lon]; no transpose needed.
            values = np.frombuffer(payload, dtype="<f4").reshape(layers, frame_nlat, frame_nlon)
            yield Frame(
                variable_id=variable_id,
                name=name,
                unit=unit,
                step=step,
                date=tuple(int(v) for v in date),
                level_index=level_index,
                layers=layers,
                nlon=frame_nlon,
                nlat=frame_nlat,
                values=values.copy(),
            )


def load_stream(path: str | Path) -> "dict[str, list[Frame]]":
    """Group every frame of one file by variable name, preserving order.

    中文说明：把单个文件中的所有帧按变量名分组，并保持原有的时间顺序。
    """
    grouped: dict[str, list[Frame]] = {}
    for frame in iter_frames(path):
        grouped.setdefault(frame.name, []).append(frame)
    return grouped


if __name__ == "__main__":
    import sys

    # 命令行用法：python frame_io.py a.frames b.frames ...
    # 逐个文件打印变量名、帧数、时刻数与数组形状。
    # CLI usage: python frame_io.py a.frames b.frames ...
    # Print variable name, frame count, number of times and array shape per file.
    for frame_path in sys.argv[1:]:
        grouped = load_stream(frame_path)
        print(frame_path)
        for name, frames in grouped.items():
            steps = sorted({frame.step for frame in frames})
            print(f"  {name:32s} {len(frames):4d} frames, {len(steps)} times, "
                  f"shape {frames[0].values.shape}, unit {frames[0].unit}")

# 数值模式 / Producer
#         │
#         │ 持续 append 二进制数据
#         ▼
# ┌──────────────────────┐
# │      .pstream 文件    │
# │                      │
# │ FILE_HEADER          │
# │ FRAME_HEADER + DATA  │
# │ FRAME_HEADER + DATA  │
# │ FRAME_HEADER + DATA  │
# │ ...                  │
# └──────────┬───────────┘
#            │
#            ▼
#     PStreamReader
#            │
#            │ read_new()
#            ▼
#       list[Frame]
#            │
#            ▼
#       FrameBuffer
#            │
#       ┌────┴────┐
#       ▼         ▼
#  latest()    history()


from __future__ import annotations

from collections import defaultdict, deque

from dataclasses import dataclass

# lru_cache 用于缓存高斯纬度权重，
# 避免每次计算 global_mean 时重复求解 Legendre 多项式根。
from functools import lru_cache

from pathlib import Path

# struct 用于按照固定的二进制布局解析文件头和帧头。
import struct

from typing import Deque

import numpy as np

from numpy.typing import NDArray

# ============================================================================
# 二进制协议结构定义
# ============================================================================

# | 格式  | 含义               |      大小 |
# | --- | ---------------- | ------: |
# | `s` | bytes 字符串        |    指定长度 |
# | `H` | unsigned short   | 2 bytes |
# | `I` | unsigned int     | 4 bytes |
# | `i` | signed int       | 4 bytes |
# | `q` | signed long long | 8 bytes |
# | `f` | float            | 4 bytes |


# 整个 pstream 文件最前面的“文件头”结构。
# struct.Struct 会预编译该二进制格式，
# 后续可以通过 FILE_HEADER.unpack(...) 高效解析。
FILE_HEADER = struct.Struct("<8sHHIIIIIII24s")


# 单个数据帧前面的“帧头”结构。
# (
#     magic,
#     version,
#     header_size,
#     variable_id,
#     payload_bytes,
#     step,
#     year,
#     month,
#     day,
#     hour,
#     minute,
#     level_index,
#     nlon,
#     nlat,
#     layers,
#     minimum,
#     maximum,
#     checksum,
#     name,
#     unit,
#     _reserved,
# ) = unpacked
FRAME_HEADER = struct.Struct("<8sHHIIqiiiiiiiiiffI32s16s4s")


# 文件级 magic number。
# 用于判断打开的文件是否真的是本协议定义的 pstream 文件。
# 固定长度为 8 字节。
FILE_MAGIC = b"PLASICPS"


# 单个 frame 的 magic number。
# 用于在流中确认当前位置确实是一个合法帧的开始。
# "\0" 是末尾的 NUL 字节，使其总长度满足协议定义。
FRAME_MAGIC = b"PFRAME1\0"

PROTOCOL_VERSION = 1


# ============================================================================
# 协议相关异常
# ============================================================================
class ProtocolError(ValueError):
    """
    pstream 二进制协议格式错误。

    继承自 ValueError，表示：
    文件本身可以被读取，但其内容不符合当前读取器所支持的协议。
    """

    pass


# ============================================================================
# 辅助函数
# ============================================================================


def _decode_fixed(value: bytes) -> str:
    """
    解码固定长度的 C 风格字符串字段。

    二进制协议中的 name / unit 等字段通常是固定长度 byte array，
    实际字符串结束位置通过 NUL（\\0）表示。

    例如：
        b"temperature\\0\\0\\0..."

    会得到：
        "temperature"

    处理流程：
    1. 以第一个 b"\\0" 为界进行切分；
    2. 只保留前面的有效字节；
    3. 按 UTF-8 解码为 Python str。
    """
    return value.split(b"\0", 1)[0].decode("utf-8")


def fnv1a(payload: bytes) -> int:
    """
    计算 payload 的 32-bit FNV-1a 哈希值。

    这里用作每个 frame payload 的 checksum，
    用于检测：
    - 文件写入过程中数据损坏；
    - 读取过程中数据截断或变化；
    - producer / consumer 对 payload 理解不一致等问题。

    FNV-1a 32-bit 算法：

        hash 初始值 = 2166136261

        对每一个 byte：
            hash ^= byte
            hash *= 16777619

        每轮结果限制为 32 bit。
    """

    # FNV-1a 32-bit offset basis。
    value = 2166136261

    # 逐字节计算哈希。
    for byte in payload:
        # 先与当前字节做 XOR。
        value ^= byte

        # 再乘 FNV prime。
        # "& 0xFFFFFFFF" 强制结果保留低 32 bit，
        # 等价于 uint32 溢出行为。
        value = (value * 16777619) & 0xFFFFFFFF

    return value


# ============================================================================
# 高斯纬度面积权重
# ============================================================================


# 函数参数 nlat
#      │
#      ▼
# ┌───────────────┐
# │ cache 中有吗？│
# └───────────────┘
#     │       │
#    有       没有
#     │       │
#     ▼       ▼
# 直接返回   执行函数
# 缓存结果   leggauss(nlat)
#             │
#             ▼
#           存入 cache
@lru_cache(maxsize=8)
def _gaussian_weights(nlat: int) -> NDArray[np.float64]:
    """
    返回 nlat 点高斯纬度对应的 Legendre-Gauss 积分权重。
    """
    # np.polynomial.legendre.leggauss 返回 (nodes, weights)，
    # nodes 是 Legendre 多项式在 [-1, 1] 上的根，
    # weights 是各纬度带在积分中的面积权重。
    _, weights = np.polynomial.legendre.leggauss(nlat)
    return weights


# ============================================================================
# 文件头数据结构
# ============================================================================


@dataclass(frozen=True)
class PStreamHeader:
    """
    pstream 文件级别的元数据。

    frozen=True 表示 dataclass 创建后不可修改，
    可以理解为一个 immutable value object。
    """

    version: int
    flags: int
    nlon: int
    nlat: int
    nlev: int
    ocean_levels: int
    float_size: int


# ============================================================================
# 单个数据帧
# ============================================================================


@dataclass(frozen=True)
class Frame:
    """
    表示 pstream 中读取出来的一个完整数据帧。
    """

    variable_id: int
    name: str
    unit: str
    step: int
    year: int
    month: int
    day: int
    hour: int
    minute: int
    level_index: int
    minimum: float
    maximum: float

    data: NDArray[np.float32]

    @property
    def layers(self) -> int:
        """
        返回当前 Frame 中包含的 layer 数量。

        因为 data 的 shape 定义为：
            (layers, nlat, nlon)

        所以 shape[0] 即 layer 数。
        """
        return int(self.data.shape[0])

    def layer(self, index: int = 0) -> NDArray[np.float32]:
        if not 0 <= index < self.layers:
            raise IndexError(index)

        return self.data[index]

    def global_mean(self, index: int = 0) -> float:
        """
        计算指定 layer 上所有网格点的纬度面积加权平均值。

        使用高斯纬度对应的 Legendre-Gauss 积分权重进行面积加权：

            global_mean = Σⱼ Σᵢ data[j,i] · wⱼ  /  Σⱼ Σᵢ wⱼ

        其中 wⱼ 是第 j 个高斯纬度的 Legendre-Gauss 权重，
        等价于气候学中标准的 cos(latitude) 面积加权全球平均。
        """

        return self.global_mean_with_coverage(index)[0]

    def global_mean_with_coverage(self, index: int = 0) -> tuple[float, float]:
        """Return the valid-area-normalized mean and represented area fraction."""
        from plasic_app.visualization.data import gaussian_global_mean

        return gaussian_global_mean(self.layer(index))


# ============================================================================
# pstream 增量读取器
# ============================================================================


class PStreamReader:
    """
    pstream 二进制文件的增量读取器。

    该类的主要目标不是一次性读取整个文件，
    而是支持某个 producer 持续向文件末尾追加 Frame 时，
    consumer 可以周期性调用 read_new() 获取新产生的数据。

    内部维护两个主要状态：

    1. header
       文件级 header。
       第一次成功读取后会缓存下来，不需要每次重新解析。

    2. offset
       已经成功消费到的文件字节位置。
       下一次 read_new() 从这里继续读取。

    一个重要设计点是：
    只有完整且校验成功的 frame 才会推进 self.offset。

    因此如果 producer 当前只写了一半 frame，
    本次 read_new() 会停止读取，
    下次调用时仍然会从该 frame 的开头重新尝试。
    """

    def __init__(self, path: Path) -> None:
        self.path = Path(path)

        self.header: PStreamHeader | None = None

        self.offset = 0

    def reset(self) -> None:
        self.header = None
        self.offset = 0

    @property
    def exhausted(self) -> bool:
        """已经消费完文件当前的全部字节（可用于判断流是否读尽）。"""
        if self.header is None:
            return False
        try:
            return self.offset >= self.path.stat().st_size
        except OSError:
            return False

    def _read_file_header(self, stream) -> bool:
        # 按 FILE_HEADER 固定长度读取文件头。
        payload = stream.read(FILE_HEADER.size)

        if len(payload) < FILE_HEADER.size:
            return False

        (
            magic,
            version,
            header_size,
            flags,
            endian,
            nlon,
            nlat,
            nlev,
            ocean_levels,
            float_size,
            _reserved,
        ) = FILE_HEADER.unpack(payload)

        if magic != FILE_MAGIC:
            raise ProtocolError("invalid pstream file magic")

        if version != PROTOCOL_VERSION:
            raise ProtocolError(f"unsupported pstream version {version}")

        if header_size != FILE_HEADER.size:
            raise ProtocolError(f"unexpected file header size {header_size}")

        if endian != 0x01020304 or float_size != 4:
            raise ProtocolError("unsupported pstream numeric layout")

        self.header = PStreamHeader(
            version=version,
            flags=flags,
            nlon=nlon,
            nlat=nlat,
            nlev=nlev,
            ocean_levels=ocean_levels,
            float_size=float_size,
        )

        # 文件头已经消费完毕。
        # 下一次 frame 读取应该从 FILE_HEADER 后面开始。
        self.offset = FILE_HEADER.size

        return True

    def read_new(self, max_frames: int | None = None) -> list[Frame]:
        """
        从 pstream 文件中读取自上次调用以来新增的完整 Frame。
        """
        if max_frames is not None and max_frames < 1:
            raise ValueError("max_frames must be positive")

        if not self.path.exists():
            return []

        frames: list[Frame] = []

        with self.path.open("rb") as stream:

            if self.header is None:
                if not self._read_file_header(stream):
                    return []

            # 跳转到上一次成功处理结束的位置。
            # 第一次读取时通常是 FILE_HEADER.size；
            # 后续则是上一帧结束的位置。
            stream.seek(self.offset)

            while True:

                frame_offset = stream.tell()

                header_payload = stream.read(FRAME_HEADER.size)

                if not header_payload:
                    break

                # 能读到一些数据，但不足一个完整 frame header。
                #
                # 通常说明 producer 正在写下一个 frame，
                # 只是 header 还没有写完。
                #
                # 不认为是协议错误，等待下次 read_new() 再处理。
                if len(header_payload) < FRAME_HEADER.size:
                    break

                # frame header 已完整，按固定协议格式解包。
                unpacked = FRAME_HEADER.unpack(header_payload)

                (
                    magic,
                    version,
                    header_size,
                    variable_id,
                    payload_bytes,
                    step,
                    year,
                    month,
                    day,
                    hour,
                    minute,
                    level_index,
                    nlon,
                    nlat,
                    layers,
                    minimum,
                    maximum,
                    checksum,
                    name,
                    unit,
                    _reserved,
                ) = unpacked

                if magic != FRAME_MAGIC:
                    raise ProtocolError(f"invalid frame magic at byte {frame_offset}")

                if version != PROTOCOL_VERSION:
                    raise ProtocolError(f"unsupported frame version {version}")

                if header_size != FRAME_HEADER.size:
                    raise ProtocolError(f"unexpected frame header size {header_size}")

                # 每个 payload element 固定为：
                expected_bytes = nlon * nlat * layers * 4

                if (
                    nlon <= 0
                    or nlat <= 0
                    or layers <= 0
                    or payload_bytes != expected_bytes
                ):
                    raise ProtocolError("invalid frame dimensions")

                # 根据 frame header 声明的 payload_bytes 读取完整的数值数据。
                payload = stream.read(payload_bytes)

                if len(payload) < payload_bytes:
                    break

                # 对完整 payload 计算 FNV-1a checksum，
                # 与 header 中存储的 checksum 比较。
                if fnv1a(payload) != checksum:
                    raise ProtocolError(f"checksum mismatch at byte {frame_offset}")

                # 将原始 payload 解释为 NumPy float32 数据。
                data = (
                    # "<f4"：
                    #   <  = little endian
                    #   f4 = 4-byte floating point，即 float32
                    np.frombuffer(payload, dtype="<f4").reshape((layers, nlat, nlon))
                    # frombuffer 默认直接引用 payload 的底层 memory。
                    # 这里显式 copy()，
                    # 让 Frame 中的 data 拥有独立 NumPy 内存，
                    # 与原始 bytes 对象解耦。
                    .copy()
                )

                # 数据完全通过验证后，
                # 构造一个 immutable Frame 对象。
                frames.append(
                    Frame(
                        variable_id=variable_id,
                        name=_decode_fixed(name),
                        unit=_decode_fixed(unit),
                        step=step,
                        year=year,
                        month=month,
                        day=day,
                        hour=hour,
                        minute=minute,
                        level_index=level_index,
                        minimum=minimum,
                        maximum=maximum,
                        data=data,
                    )
                )

                self.offset = stream.tell()

                if max_frames is not None and len(frames) >= max_frames:
                    break

        return frames


# ============================================================================
# 内存中的 Frame 历史缓存
# ============================================================================


class FrameBuffer:
    """
    按变量名称维护最近若干个 Frame 的内存缓存。

    数据结构大致为：

        {
            "temperature": deque([...]),
            "pressure": deque([...]),
            "humidity": deque([...]),
        }

    每个变量拥有独立的 deque，
    且 deque 的最大长度由 max_frames_per_variable 控制。

    当新的 Frame 持续加入并超过 maxlen 时，
    deque 会自动丢弃最旧的 Frame。
    """

    def __init__(self, max_frames_per_variable: int = 128) -> None:

        if max_frames_per_variable < 1:
            raise ValueError("buffer length must be positive")

        self.max_frames_per_variable = max_frames_per_variable

        self._frames: dict[str, Deque[Frame]] = defaultdict(
            lambda: deque(maxlen=self.max_frames_per_variable)
        )

    def append(self, frame: Frame) -> None:
        """
        向缓存中加入一个 Frame。
        Frame 按 frame.name 分类。
        """

        # 例如 frame.name == "temperature"，
        # 则该 frame 被加入 temperature 对应的 deque。
        #
        # 如果超过 maxlen，
        # deque 自动删除最旧元素。
        self._frames[frame.name].append(frame)

    def extend(self, frames: list[Frame]) -> None:

        for frame in frames:
            self.append(frame)

    def latest(self, name: str) -> Frame | None:
        """
        获取指定变量的最新一个 Frame。

        如果该变量尚无任何缓存数据，则返回 None。
        """

        values = self._frames.get(name)

        return values[-1] if values else None

    def variables(self) -> tuple[str, ...]:
        """
        返回当前缓存中已经存在的所有变量名称。

        返回值：
            排序后的 tuple[str, ...]
        """
        return tuple(sorted(self._frames))

    def history(self, name: str) -> tuple[Frame, ...]:
        """
        返回指定变量当前缓存的全部历史 Frame。

        顺序：
            最旧 -> 最新

        如果变量不存在，则返回空 tuple。
        """

        # 同样使用 get()，避免查询不存在变量时
        # 意外触发 defaultdict 创建新的 key。
        return tuple(self._frames.get(name, ()))

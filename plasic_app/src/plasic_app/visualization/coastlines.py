from __future__ import annotations

# lru_cache 用于缓存函数返回结果。
# 这里两个函数的 maxsize 都是 1，因此第一次调用完成后，
# 后续调用会直接复用缓存结果，避免重复读取资源文件和重复构造数组。
from functools import lru_cache

# importlib.resources.files 用于以“包资源”的方式访问随 Python 包一起分发的文件。
# 相比直接拼接文件系统路径，这种方式对普通安装、wheel 等打包形式更加友好。
from importlib.resources import files

import json

import numpy as np

from numpy.typing import NDArray


# 缓存 coastline_segments() 的结果。
@lru_cache(maxsize=1)
def coastline_segments() -> tuple[NDArray[np.float32], ...]:
    resource = files("plasic_app.resources").joinpath(
        "coastlines_110m.json"
    )

    payload = json.loads(resource.read_text(encoding="utf-8"))

    if payload.get("schema_version") != 1:
        raise ValueError("unsupported coastline resource version")

    # 将 JSON 中的每一条 coastline polyline 转成 NumPy 数组。
    # 返回值是 tuple，其中每个元素代表一条独立的海岸线线段：
    # segment.shape 通常类似：
    #     (N, 2)
    # 每一行的两个元素分别是：
    #     [longitude, latitude]
    # dtype 被显式转换为 float32，以降低内存占用，
    # 并让后续绘图使用统一的数据类型。
    return tuple(
        np.asarray(line, dtype=np.float32)
        for line in payload["lines"]
        if len(line) >= 2
    )


# 同样缓存扁平化后的 coastline。
# flattened_coastlines() 本身又依赖已经缓存的 coastline_segments()，
# 因而完整的“读取 JSON -> 解析 -> 展平”流程通常只执行一次。
@lru_cache(maxsize=1)
def flattened_coastlines() -> tuple[NDArray[np.float32], NDArray[np.float32]]:

    longitudes: list[float] = []

    latitudes: list[float] = []

    # coastline_segments() 返回多个彼此独立的 coastline segment。
    # 每个 segment 是一个二维 NumPy 数组，
    # 每一行是：
    #     [longitude, latitude]
    for segment in coastline_segments():

        # previous 保存当前 segment 中上一个点的经度。
        # 初始为 None，因为处理第一个点时不存在前驱点。
        # 这里主要用于检测相邻两个坐标点是否跨越了 ±180° 经线
        # （即国际日期变更线附近的经度跳变）。
        previous: float | None = None

        for longitude, latitude in segment:
            current = float(longitude)

            # 如果不是 segment 的第一个点，
            # 就检查当前经度与前一个经度的差值。
            #
            # 正常情况下，相邻 coastline 点的经度变化不会特别大。
            # 但如果一条地理线穿过 ±180° 经线，例如：
            #
            #     previous = 179°
            #     current  = -179°
            #
            # 几何意义上只跨越了约 2°，
            # 但直接做经度数值差：
            #
            #     abs(-179 - 179) = 358°
            #
            # 如果直接把这两个点交给普通二维绘图库，
            # 绘图库可能会从 +179° 横跨整张地图连接到 -179°，
            # 产生一条错误的超长水平线。
            if previous is not None and abs(current - previous) > 180.0:

                # 在经度数组中插入 NaN。
                # 很多绘图库（例如基于连续 x/y 数组绘制折线的组件）
                # 会把 NaN 当成“断开当前折线”的标记。
                longitudes.append(float("nan"))

                # 纬度数组必须同步插入 NaN，
                latitudes.append(float("nan"))

            # 添加当前点的经度。
            longitudes.append(current)

            # 添加当前点的纬度。
            latitudes.append(float(latitude))

            # 当前经度保存为 previous，
            # 供下一轮循环判断经度是否出现 >180° 的跳变。
            previous = current

        # 一个 segment 处理完成后，再插入一组 NaN。
        #
        # 目的与上面的日期变更线断开类似：
        # 告诉绘图库当前 coastline segment 已结束，
        # 下一条 segment 不应该与当前 segment 的最后一个点连起来。
        #
        # 因此最后可以把所有 coastline 都放进同一个 PlotDataItem，
        # 而无需为每个 segment 分别创建一个绘图对象。
        longitudes.append(float("nan"))
        latitudes.append(float("nan"))

    # 将两个 Python list 转换成连续的 float32 NumPy 数组。
    # 最终得到：
    # longitudes = [lon0, lon1, ..., NaN, lonN, ..., NaN, ...]
    # latitudes  = [lat0, lat1, ..., NaN, latN, ..., NaN, ...]
    # 两个数组长度完全相同。
    # NaN 同时承担两种“断线”用途：
    # 1. 分隔不同 coastline segment；
    # 2. 防止跨 ±180° 经线时产生错误的地图横跨连线。
    return (
        np.asarray(longitudes, dtype=np.float32),
        np.asarray(latitudes, dtype=np.float32),
    )
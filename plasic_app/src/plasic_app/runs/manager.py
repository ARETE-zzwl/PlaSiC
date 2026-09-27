from __future__ import annotations

from dataclasses import dataclass

from datetime import datetime, timezone

import json

# 这里主要使用 os.link() 创建硬链接。
import os

from pathlib import Path

# 当硬链接创建失败时，使用 shutil.copy2() 复制 restart 文件。
import shutil

# uuid4 用于给每次运行生成随机、低碰撞概率的唯一标识。
from uuid import uuid4

# 说明：
# 运行目录只保存实验本身产生的结果（config.json、progress.txt、
# frames.pstream、model.log、monthly/ 等）。
# C runtime 所需的运行参数一律通过命令行直接传递，
# 因此不再生成 model_config.json 和 build.json。

# 实验配置对象。
#
# ExperimentConfig 负责保存和验证一次 PlaSiC experiment/run 所需要的完整配置。
from plasic_app.config.models import ExperimentConfig

from plasic_app.paths import default_runs_root


# frozen=True 表示 dataclass 实例创建之后字段不可重新赋值，
# 使 PreparedRun 更接近一个不可变的“运行描述对象”。
@dataclass(frozen=True)
class PreparedRun:
    """
    描述一次已经准备完成的 PlaSiC 运行。
    """

    # 本次运行的唯一标识。
    # 当前实现中格式类似：
    #   20260818T123456Z-a1b2c3d4
    # 前半部分是 UTC 时间戳，
    # 后半部分是 UUID 的前 8 个十六进制字符。
    run_id: str

    # 本次运行的独立工作目录。
    directory: Path

    initial_restart: Path | None

    @property
    def progress(self) -> Path:

        # progress.txt 使用“一行一个 JSON 对象”的格式，
        # 由 C runtime 在模型运行过程中持续追加进度记录。
        return self.directory / "progress.txt"

    @property
    def frames(self) -> Path:
        """
        返回本次运行的 frame 数据流文件路径。

        文件名固定为：
            frames.pstream
        """

        return self.directory / "frames.pstream"

    @property
    def log(self) -> Path:
        """
        返回本次模型运行日志文件路径。

        文件名固定为：
            model.log
        """

        return self.directory / "model.log"


class RunManager:
    """
    管理 PlaSiC 每次独立运行所使用的工作目录和元数据文件。
    """

    def __init__(self, root: Path | None = None) -> None:

        # 如果调用者传入 root，就使用该路径；
        # 否则使用项目定义的默认 runs 目录。
        self.root = (root or default_runs_root()).resolve()

    @staticmethod
    def _write_json(path: Path, payload: object) -> None:
        """
        将一个 Python 对象以统一格式写入 JSON 文件。
        """

        # json.dumps() 首先将 Python 对象序列化成字符串，
        # 然后通过 Path.write_text() 一次性写入文件。
        path.write_text(
            json.dumps(
                payload,
                indent=2,
                sort_keys=True,
                # 不将中文等非 ASCII 字符转义成 \uXXXX。
                ensure_ascii=False
            )
            # 确保 JSON 文件最后包含换行符。
            + "\n",
            encoding="utf-8",
        )

    def prepare(
        self,
        config: ExperimentConfig,
        c_root: Path,
    ) -> PreparedRun:
        """
        为一次 PlaSiC 模型运行准备独立的运行目录。

        返回：
            PreparedRun:
                描述已经准备完成的运行目录和关键文件。

        主要流程：
            1. 校验配置；
            2. 生成唯一 run ID；
            3. 创建运行目录；
            4. 根据需要准备 restart 文件；
            5. 写入 config.json；
            6. 返回 PreparedRun。

        说明：
            C runtime 的运行参数由 run_invocation() 通过命令行传递，
            因此这里不再写入 model_config.json / build.json。
        """

        # 首先验证 ExperimentConfig。
        #
        # validate() 预期返回错误字符串列表；
        # 空列表意味着配置合法。
        errors = config.validate()

        # 如果存在任何配置错误，
        # 则不创建运行任务，直接抛出 ValueError。
        if errors:
            raise ValueError("\n".join(errors))

        # 获取当前 UTC 时间，并格式化为适合作为目录名的形式。
        # 格式：
        #   YYYYMMDDTHHMMSSZ
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

        # 为本次运行生成唯一 ID。
        #
        # UUID4 是随机 UUID；
        # .hex 得到不包含 "-" 的 32 位十六进制字符串；
        # [:8] 只取前 8 位。
        #
        # 最终 run_id 例如：
        #
        #   20260818T125830Z-4f82a12c
        #
        # 时间戳便于人工排序和识别，
        # UUID 后缀则降低同一秒内创建多个 run 时发生重名的概率。
        run_id = f"{timestamp}-{uuid4().hex[:8]}"

        # 每个 run 使用一个独立目录：
        #   <runs_root>/<run_id>
        directory = self.root / run_id

        directory.mkdir(parents=True, exist_ok=False)
        (directory / "monthly").mkdir()

        # 默认没有准备 restart 文件。
        # None 表示后续将执行 cold start。
        initial_restart: Path | None = None

        # 空的 restart path 表示 cold start。
        if config.run.restart:

            # 将配置中的 restart 路径解析成实际模型文件路径。
            # _resolve_model_path() 会结合 c_root
            # 处理相对于模型根目录的路径。
            source_restart = config._resolve_model_path(
                config.run.restart, c_root
            )

            # restart 文件必须实际存在，并且必须是普通文件。
            if not source_restart.is_file():

                # 如果 restart 不存在，
                # 直接报告缺失文件的位置。
                raise FileNotFoundError(source_restart)

            # 在本次独立运行目录中，
            # 将 restart 统一命名为：
            #
            #   initial.restart
            initial_restart = directory / "initial.restart"

            try:
                # 优先创建硬链接。
                #
                # 优点是：
                # - 不需要复制大型 restart 文件；
                # - 创建速度快；
                # - 不额外占用一份文件内容的磁盘空间。
                #
                # source_restart 和 initial_restart
                # 此时会指向同一底层 inode/文件数据。
                os.link(source_restart, initial_restart)

            except OSError:
                # 硬链接可能因为多种原因失败
                # 这种情况下退化为普通文件复制。
                #
                # copy2() 相比 copy()，
                # 还会尽可能保留文件的 metadata，
                # 例如修改时间等。
                shutil.copy2(source_restart, initial_restart)

        # 保存用户/应用层完整的 ExperimentConfig，
        # 作为这次运行的意图记录。
        # C runtime 的运行参数不在这里落盘：
        # run_invocation() 会把它们直接拼进命令行。
        self._write_json(directory / "config.json", config.to_dict())

        # 所有准备工作完成后，
        # 返回一个不可变的 PreparedRun 对象。
        return PreparedRun(
            run_id=run_id,
            directory=directory,
            initial_restart=initial_restart,
        )

    def list_runs(self) -> list[Path]:
        """
        列出当前 runs root 中所有有效的运行目录。

        判断一个目录是不是有效 run 的条件是：
        1. 它本身是一个目录；
        2. 目录中存在 config.json。

        返回：
            list[Path]:
                按路径名称逆序排序后的 run directory 列表。
        """

        if not self.root.exists():
            return []

        # 遍历 runs root 的一级子项，
        # 过滤出符合 run directory 特征的路径。
        return sorted(
            (
                path
                for path in self.root.iterdir()
                if path.is_dir()
                and (path / "config.json").is_file()
            ),

            # 按 Path 的自然排序结果逆序返回。
            #
            # 由于 run_id 以：
            #
            #   YYYYMMDDTHHMMSSZ
            #
            # 开头，因此在当前命名规则下，
            # 逆序通常也就意味着最新的 run 排在最前面。
            reverse=True,
        )

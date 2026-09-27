
#                     ProcessController
#                            │
#           ┌────────────────┴────────────────┐
#           │                                 │
#           ▼                                 ▼
#    build_process                       run_process
#      QProcess                            QProcess
#           │                                 │
#           │ start                           │ start
#           ▼                                 ▼
#        compiler                          PlaSiC
#           │                                 │
#           │ output                          │ output
#           ▼                                 ▼
# _read_build_output()              _read_run_output()
#           │                                 │
#           │                                 ├──> 写 log 文件
#           │                                 │
#           └────────────┬────────────────────┘
#                        ▼
#                 logReceived
#                        │
#                        ▼
#                      GUI


# 运行结束：

# QProcess.finished
#        │
#        ├── build → _build_finished()
#        │
#        └── run   → _run_finished()
#                         │
#                         ▼
#                 stateChanged("idle")

from __future__ import annotations

import os
from pathlib import Path
import signal

from PySide6.QtCore import (
    QObject,
    QProcess,
    QProcessEnvironment,
    QTimer,
    Signal,
)

from plasic_app.controller.commands import (
    BuildInvocation,
    RunInvocation,
)


class ProcessController(QObject):
    """
    负责统一管理“构建（build）”和“运行（run）”两个外部进程。

    这个 Controller 基于 Qt 的 QProcess，因此外部进程的启动、输出读取、
    结束通知等都是异步的，不会阻塞 Qt GUI 主线程。

    对外通过 Qt Signal 暴露以下事件：

    - logReceived:
        收到构建或运行进程输出时触发。
        参数为日志文本。

    - stateChanged:
        Controller 状态发生变化时触发。
        当前代码中可能出现的状态包括：
        "building"、"running"、"stopping"、"idle"。

    - buildFinished:
        build 结束时触发。
        第一个参数表示是否成功，
        第二个参数为错误描述；成功时为空字符串。

    - runFinished:
        run 结束时触发。
        第一个参数表示是否成功，
        第二个参数为进程退出码。
    """

    logReceived = Signal(str)

    stateChanged = Signal(str)

    buildFinished = Signal(bool, str)

    runFinished = Signal(bool, int)

    def __init__(self, parent: QObject | None = None) -> None:
        """
        初始化进程控制器。

        这里分别创建两个 QProcess：

        - build_process:
            专门执行构建命令。

        - run_process:
            专门执行实际模型/程序运行命令。

        两个 QProcess 都以 self 为 Qt parent，因此 ProcessController
        被销毁时，Qt 的对象树机制也会负责管理这些对象的生命周期。
        """
        super().__init__(parent)

        # 独立的构建进程。
        self.build_process = QProcess(self)

        # 独立的运行进程。
        self.run_process = QProcess(self)

        # 当前 run 对应的日志文件对象。
        #
        # 仅 start_run() 后有效；
        # run 结束后会在 _run_finished() 中关闭并恢复为 None。
        self._run_log = None

        # 标记当前 run 是否是由用户主动 stop 引起的。
        #
        # 这和 exit_code 配合使用：
        # 即使进程在 SIGTERM 后最终返回 0，
        # 如果是用户主动停止，也不会被判断成正常成功运行。
        self._stopping = False

        # 将 build 进程的 stdout 和 stderr 合并。
        # 这样后续只需要读取 StandardOutput，
        # stderr 中的信息也会一并读到。
        self.build_process.setProcessChannelMode(
            QProcess.ProcessChannelMode.MergedChannels
        )

        # 当 build 进程有新的输出可读取时，
        # 调用 _read_build_output()。
        self.build_process.readyReadStandardOutput.connect(
            self._read_build_output
        )

        # build 进程结束后调用 _build_finished()。
        self.build_process.finished.connect(self._build_finished)

        # run 进程也采用 stdout/stderr 合并模式。
        self.run_process.setProcessChannelMode(
            QProcess.ProcessChannelMode.MergedChannels
        )

        # run 产生新输出时进行读取。
        self.run_process.readyReadStandardOutput.connect(
            self._read_run_output
        )

        # run 结束时执行统一清理以及状态通知。
        self.run_process.finished.connect(self._run_finished)

    @staticmethod
    def _environment(values: dict[str, str]) -> QProcessEnvironment:
        """
        将 Python 字典转换成 QProcessEnvironment。

        参数
        ----------
        values:
            环境变量映射，例如：

            {
                "PATH": "...",
            }

        返回
        ----
        QProcessEnvironment
            可直接通过 QProcess.setProcessEnvironment() 设置给子进程。

        注意
        ----
        这里创建的是一个新的 QProcessEnvironment，
        并逐项插入传入的环境变量。

        也就是说，该函数本身没有调用
        QProcessEnvironment.systemEnvironment() 来自动复制当前系统环境。
        最终子进程拥有哪些环境变量，取决于 invocation.environment
        中实际传入的内容。
        """
        environment = QProcessEnvironment()

        # 将 invocation 中准备好的环境变量逐项写入 Qt 环境对象。
        for key, value in values.items():
            environment.insert(key, value)

        return environment

    @staticmethod
    def _decode(process: QProcess) -> str:
        """
        从指定 QProcess 中读取当前所有可用的标准输出，并解码成字符串。

        由于前面已经设置：
            ProcessChannelMode.MergedChannels

        所以这里虽然调用的是 readAllStandardOutput()，
        实际上可以同时获取 stdout 和 stderr 的内容。

        UTF-8 解码时使用 errors="replace"：
        如果外部程序输出了非法 UTF-8 字节，不会抛异常，
        而是使用替换字符代替，从而避免 GUI 因日志解码问题崩溃。
        """
        return bytes(process.readAllStandardOutput()).decode(
            "utf-8", errors="replace"
        )

    @property
    def build_running(self) -> bool:
        """
        返回 build 进程当前是否仍然处于活动状态。

        QProcess 状态通常包括：

        - NotRunning
        - Starting
        - Running

        因此这里不仅 Running 算作正在运行，
        Starting 阶段也会返回 True。

        这样可以避免进程刚开始启动、尚未真正进入 Running 时，
        又重复启动另一个任务。
        """
        return (
            self.build_process.state()
            != QProcess.ProcessState.NotRunning
        )

    @property
    def run_running(self) -> bool:
        """
        返回 run 进程当前是否仍然处于活动状态。
        """
        return (
            self.run_process.state()
            != QProcess.ProcessState.NotRunning
        )

    def start_build(self, invocation: BuildInvocation) -> None:
        """
        启动一次 build。

        invocation 中包含启动 build 所需的信息

        Controller 明确限制：
        build 和 run 不能同时执行。
        """

        # build/run 任一进程已经启动时，都拒绝再次启动新进程。
        # 这样保证当前 Controller 始终最多只管理一个活动任务。
        if self.build_running or self.run_running:
            raise RuntimeError("another process is already active")

        # 通知 UI 或其他监听者当前进入构建状态。
        self.stateChanged.emit("building")

        # 将即将执行的命令首先输出到 GUI 日志。
        #
        # 前缀 "$ " 是类似 shell prompt 的视觉效果。
        #
        # invocation.as_list() 应返回完整命令，然后用空格拼接用于展示。
        self.logReceived.emit(
            "$ " + " ".join(invocation.as_list()) + "\n"
        )

        # 设置 build 子进程的工作目录。
        self.build_process.setWorkingDirectory(
            str(invocation.working_directory)
        )

        # 设置 build 子进程使用的环境变量。
        self.build_process.setProcessEnvironment(
            self._environment(invocation.environment)
        )

        # 异步启动 build。
        #
        # QProcess.start() 调用后立即返回，
        # 不会等待外部进程结束。
        #
        # 后续输出通过 readyReadStandardOutput signal 获取，
        # 结束通过 finished signal 获取。
        self.build_process.start(
            invocation.program, list(invocation.arguments)
        )

    def start_run(
        self, invocation: RunInvocation, log_path: Path
    ) -> None:
        """
        启动一次正式运行任务，并同步将输出写入日志文件。

        与 build 相比，run 多了持久化日志功能：

        1. 命令本身写入 log_path；
        2. 运行过程中 stdout/stderr 同时：
           - 发送给 GUI；
           - append 到日志文件。

        参数
        ----------
        invocation:
            描述 run 命令的 RunInvocation。

        log_path:
            本次 run 对应的日志文件路径。
        """

        # 与 start_build() 相同，不允许同时存在 build/run。
        if self.build_running or self.run_running:
            raise RuntimeError("another process is already active")

        log_path.parent.mkdir(parents=True, exist_ok=True)

        # 以 append 模式打开日志文件。
        #
        # 使用 "a" 而不是 "w"：
        # 如果日志文件已经存在，不会覆盖旧内容，
        # 新内容会追加到文件尾部。
        self._run_log = log_path.open("a", encoding="utf-8")

        # 新 run 刚开始时，重置“用户主动停止”状态。
        self._stopping = False

        # 通知 UI 当前进入 running 状态。
        self.stateChanged.emit("running")

        # 生成一条 shell 风格的命令日志。
        command = "$ " + " ".join(invocation.as_list()) + "\n"

        # 将命令展示到 UI 日志。
        self.logReceived.emit(command)

        # 同时将命令写入持久化日志文件。
        self._run_log.write(command)

        # 立即 flush，确保命令尽快真正写入文件。
        # 如果程序后续异常退出，至少已经记录了执行的具体命令。
        self._run_log.flush()

        # 设置 run 子进程工作目录。
        self.run_process.setWorkingDirectory(
            str(invocation.working_directory)
        )

        # 设置 run 子进程环境变量。
        self.run_process.setProcessEnvironment(
            self._environment(invocation.environment)
        )

        # 异步启动 run。
        self.run_process.start(
            invocation.program, list(invocation.arguments)
        )

    def _read_build_output(self) -> None:
        """
        处理 build 进程新产生的输出。

        该方法由：

            build_process.readyReadStandardOutput

        signal 自动触发。

        这里仅负责：
        读取 -> 解码 -> 转发到 logReceived。

        build 日志当前没有像 run 那样持久化到独立日志文件。
        """

        # 读取当前已经缓存的全部 build 输出。
        text = self._decode(self.build_process)

        # 避免发送空字符串 signal。
        if text:
            self.logReceived.emit(text)

    def _read_run_output(self) -> None:
        """
        处理 run 进程新产生的输出。

        与 build 不同，run 输出有两个目的地：

        1. 通过 logReceived signal 发给 UI；
        2. 同时写入本次 run 对应的日志文件。
        """

        # 读取 run 当前所有可用输出。
        text = self._decode(self.run_process)

        # 没有内容时直接结束。
        if not text:
            return

        # 将内容发送给 GUI / 其他监听器。
        self.logReceived.emit(text)

        # 正常情况下，run 启动时已经打开日志文件。
        if self._run_log is not None:
            # 将进程输出追加写入文件。
            self._run_log.write(text)

            # 立即刷新，方便用户在运行过程中实时查看日志，
            # 也减少异常退出时尚未落盘的数据。
            self._run_log.flush()

    def _build_finished(
        self, exit_code: int, _status: QProcess.ExitStatus
    ) -> None:
        """
        build 进程结束后的回调。

        参数由 QProcess.finished signal 提供：

        exit_code:
            外部程序返回的退出码。

        _status:
            QProcess.ExitStatus，例如 NormalExit / CrashExit。

        当前实现没有使用 _status，因此变量名前加下划线表示
        “有这个参数，但有意不使用”。

        当前成功判定规则非常直接：

            exit_code == 0 -> 成功
            exit_code != 0 -> 失败
        """

        # Unix/CLI 程序通常约定退出码 0 表示正常成功。
        success = exit_code == 0

        # 无论成功失败，build 已结束，所以 Controller 恢复 idle。
        self.stateChanged.emit("idle")

        # 将结果通知上层。
        self.buildFinished.emit(
            success,
            "" if success else f"Build exited with code {exit_code}",
        )

    def _run_finished(
        self, exit_code: int, _status: QProcess.ExitStatus
    ) -> None:
        """
        run 进程结束后的回调。

        主要负责三个工作：

        1. 关闭日志文件；
        2. 根据退出码以及 _stopping 判断运行是否成功；
        3. 恢复 Controller 状态。

        注意 run 的成功规则与 build 略有区别：

            success = exit_code == 0 and not self._stopping

        也就是说：

        - exit_code == 0 且自然结束 -> 成功；
        - exit_code != 0 -> 失败；
        - 即使 exit_code == 0，但这是用户主动 stop 的 -> 仍视为失败/未正常完成。
        """

        # 如果本次 run 打开了日志文件，则在进程结束后关闭。
        if self._run_log is not None:
            self._run_log.close()
            self._run_log = None

        # 只有“退出码为 0 + 非主动停止”才算成功完成。
        success = exit_code == 0 and not self._stopping

        # run 已经结束，Controller 恢复空闲状态。
        self.stateChanged.emit("idle")

        # 将运行结果和真实 exit code 通知上层。
        self.runFinished.emit(success, exit_code)

        # 清理停止标记，为下一次 run 做准备。
        self._stopping = False

    def stop_run(self) -> None:
        """
        尝试优雅停止当前 run 进程。
        """

        # 当前没有 run 时，stop 是 no-op。
        if not self.run_running:
            return

        # 标记这是用户主动停止。
        # 后续 _run_finished() 会利用该标记，
        # 避免把主动停止的任务判断为成功完成。
        self._stopping = True

        # 通知 UI 当前进入 stopping 状态。
        self.stateChanged.emit("stopping")

        # 获取 QProcess 对应操作系统层面的 PID。
        self_pid = self.run_process.processId()

        # 转成普通 int，方便传给 os.killpg。
        pid = int(self_pid)

        # 用于记录是否已经成功向进程组发送 SIGTERM。
        terminated_group = False

        # killpg 是 Unix/POSIX 平台相关能力。
        #
        # 先确认：
        #   1. PID 有效；
        #   2. 当前 Python/os 模块提供 killpg。
        if pid > 0 and hasattr(os, "killpg"):
            try:
                # 向 PID 对应的进程组发送 SIGTERM。
                os.killpg(pid, signal.SIGTERM)

                # 如果没有异常，则认为已经成功尝试终止整个进程组。
                terminated_group = True

            except ProcessLookupError:
                # PID / 进程组已经不存在。
                #
                # 说明目标可能恰好已经结束，因此直接返回。
                return

            except OSError:
                # killpg 因其他系统错误失败。
                #
                # 此时不直接放弃，而是走下面的 QProcess.terminate()
                # fallback。
                terminated_group = False

        # 如果没有成功通过 killpg 停止进程组，
        # 则只针对 QProcess 所管理的进程调用 terminate()。
        if not terminated_group:
            self.run_process.terminate()

        # 设置一个 3 秒后的兜底强杀。
        #
        # 如果进程收到 SIGTERM / terminate 后在 3 秒内退出，
        # _force_stop() 中会发现 run_running == False，从而什么都不做。
        #
        # 如果 3 秒后仍未退出，则升级到强制终止。
        QTimer.singleShot(3000, self._force_stop)

    def _force_stop(self) -> None:
        """
        stop_run() 的强制终止阶段。

        该方法通常在发送 SIGTERM / terminate 3 秒后执行。

        如果进程此时仍然存活：

        1. 优先尝试向进程组发送 SIGKILL；
        2. 如果失败，则调用 QProcess.kill()。

        SIGKILL 与 SIGTERM 的关键区别：

        - SIGTERM:
            可以被程序捕获、处理或忽略，允许优雅清理。

        - SIGKILL:
            不能被目标程序捕获或忽略，由操作系统直接终止。
        """

        # 如果进程在之前 3 秒内已经退出，就无需继续处理。
        if not self.run_running:
            return

        # 再次读取当前运行进程 PID。
        pid = int(self.run_process.processId())

        # 在支持进程组信号的平台上，
        # 优先尝试杀掉整个进程组。
        if pid > 0 and hasattr(os, "killpg"):
            try:
                # SIGKILL 强制终止整个进程组。
                os.killpg(pid, signal.SIGKILL)

                # 信号成功发送后直接返回，
                # 无需再调用 QProcess.kill()。
                return
            except OSError:
                # 进程组强杀失败时，继续 fallback。
                pass

        # 最终兜底：
        # 强制结束 QProcess 直接管理的进程。
        self.run_process.kill()
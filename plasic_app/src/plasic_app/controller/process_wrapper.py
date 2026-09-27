from __future__ import annotations

# - os.setsid(): 创建新的 session
# - os.execvpe(): 用指定命令替换当前 Python 进程。
#   exec + vector arguments + search path + custom environment
# - os.environ: 获取当前进程的环境变量
import os

import sys


def _start_session() -> None:
    """Create a session unless the launcher already made us group leader.

    POSIX ``setsid`` rejects a process-group leader with EPERM. QProcess and
    terminal launchers may already put the wrapper in exactly that state; it
    is still safe for ProcessController to signal the group by this PID.
    """
    if not hasattr(os, "setsid"):
        return
    try:
        os.setsid()
    except PermissionError:
        if not (
            hasattr(os, "getpgrp")
            and hasattr(os, "getpid")
            and os.getpgrp() == os.getpid()
        ):
            raise


def main() -> int:
    """
    作为 process wrapper 的主入口。

    预期调用方式类似：

        python process_wrapper.py <command> [arg1 arg2 ...]

    该 wrapper 的主要职责是：

    1. 检查是否提供了需要执行的目标命令；
    2. 如果当前操作系统支持 setsid()，则创建一个新的 session；
    3. 使用 execvpe() 将当前 Python 进程直接替换成目标命令。

    返回值：
        int:
            2   - 没有提供需要执行的命令；
            127 - execvpe() 之后的理论兜底返回值。

    注意：
        正常情况下，如果 os.execvpe() 执行成功，
        当前 Python 进程会被新程序完全替换，
        因此不会继续执行到 `return 127`。
    """

    # sys.argv 包含当前 Python 程序的全部命令行参数。
    #
    # 例如：
    #
    #     python process_wrapper.py python worker.py --foo bar
    #
    # 此时 sys.argv 大致为：
    #
    #     [
    #         "process_wrapper.py",
    #         "python",
    #         "worker.py",
    #         "--foo",
    #         "bar",
    #     ]
    #
    # 因此：
    # - sys.argv[0] 是当前 wrapper 脚本本身；
    # - sys.argv[1] 才是需要实际执行的目标命令。
    #
    # 如果参数数量小于 2，说明用户没有提供目标命令。
    if len(sys.argv) < 2:
        print("process_wrapper requires a command", file=sys.stderr)
        return 2

    # 创建独立 session；如果启动器已经建好了同 PID 的进程组，
    # 则保留该进程组，以便 ProcessController 后续整组停止。
    _start_session()

    # 使用 execvpe() 执行目标命令。
    os.execvpe(sys.argv[1], sys.argv[1:], os.environ)

    # 理论上的兜底退出码。
    #
    # 按照 exec 系列函数的语义：
    #
    # - 成功：不会返回到这里；
    # - 失败：通常会抛出异常。
    #
    # 因此在正常 Python 语义下，这一行基本不可达（unreachable）。
    #
    # 这里仍然保留原代码中的 return 127，
    # 没有对实现做任何修改。
    return 127


if __name__ == "__main__":
    raise SystemExit(main())

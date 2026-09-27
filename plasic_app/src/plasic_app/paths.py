
from __future__ import annotations

import os
from pathlib import Path


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def repository_root() -> Path:
    return project_root().parent


def discover_c_root() -> Path:
    override = os.environ.get("PLASIC_C_ROOT")
    candidate = (
        Path(override).expanduser()
        if override
        else repository_root() / "src"
    )
    candidate = candidate.resolve()
    if not (candidate / "Makefile").is_file():
        raise FileNotFoundError(
            f"Cannot find standalone C PlaSiC at {candidate}. "
            "Set PLASIC_C_ROOT to the directory that contains the Makefile."
        )
    return candidate


def repository_build_root() -> Path:
    return repository_root() / "build"


def app_build_root() -> Path:
    return repository_build_root() / "app"


def default_runs_root() -> Path:
    return project_root() / "runs"



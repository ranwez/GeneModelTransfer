#!/usr/bin/env python3
"""Record provenance information for LRRtransfer workflow runs."""

import importlib.metadata
import platform
import socket
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

import yaml


def create_run_id() -> str:
    return datetime.now().astimezone().strftime("%Y%m%dT%H%M%S%f%z")


def _get_git_info(repo_root: Path) -> tuple[Optional[str], Optional[bool]]:
    """Return the current Git commit and dirty status."""
    try:
        commit = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()

        status = subprocess.run(
            ["git", "-C", str(repo_root), "status", "--porcelain"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return None, None

    return commit or None, bool(status.strip())


def _get_snakemake_version() -> str:
    try:
        return importlib.metadata.version("snakemake")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _write_yaml(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(path.suffix + ".tmp")

    with temporary_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(
            data,
            handle,
            sort_keys=False,
            allow_unicode=True,
        )

    temporary_path.replace(path)


def _append_workflow_log(
    log_path: Path,
    run_id: str,
    status: str,
    timestamp: datetime,
) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(
            f"{timestamp.isoformat(timespec='seconds')}\t"
            f"{run_id}\t{status}\n"
        )


def start_run_record(
    *,
    run_id: str,
    history_path: str,
    latest_path: str,
    log_path: str,
    config: dict[str, Any],
    repo_root: str,
    argv: list[str],
    working_directory: str,
) -> None:
    started_at = datetime.now().astimezone()
    git_commit, git_dirty = _get_git_info(Path(repo_root))

    data = {
        "run": {
            "id": run_id,
            "status": "running",
            "started_at": started_at.isoformat(timespec="seconds"),
            "finished_at": None,
            "duration_seconds": None,
        },
        "software": {
            "git_commit": git_commit,
            "git_dirty": git_dirty,
            "snakemake_version": _get_snakemake_version(),
            "python_version": platform.python_version(),
        },
        "execution": {
            "command": argv,
            "working_directory": working_directory,
            "hostname": socket.gethostname(),
        },
        "resolved_config": config,
    }

    _write_yaml(Path(history_path), data)
    _write_yaml(Path(latest_path), data)

    _append_workflow_log(
        Path(log_path),
        run_id,
        "STARTED",
        started_at,
    )


def finish_run_record(
    *,
    history_path: str,
    latest_path: str,
    log_path: str,
    status: str,
) -> None:
    history = Path(history_path)

    if not history.is_file():
        return

    with history.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle)

    finished_at = datetime.now().astimezone()
    started_at = datetime.fromisoformat(data["run"]["started_at"])

    data["run"].update(
        {
            "status": status,
            "finished_at": finished_at.isoformat(timespec="seconds"),
            "duration_seconds": round(
                (finished_at - started_at).total_seconds(),
                3,
            ),
        }
    )

    _write_yaml(history, data)
    _write_yaml(Path(latest_path), data)

    _append_workflow_log(
        Path(log_path),
        data["run"]["id"],
        status.upper(),
        finished_at,
    )

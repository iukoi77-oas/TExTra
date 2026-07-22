"""Run prep external commands with optional debug log capture."""

import os
import shlex
import subprocess


def _decode_stream(value):
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _format_command(command):
    return " ".join(shlex.quote(str(part)) for part in command)


def _append_extend_log(args, section, command, cwd, stdout_text="", stderr_text=""):
    log_dir = os.path.join(args.out_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, "prep_extend.log")
    lines = [f"========={section}==========\n"]
    if cwd:
        lines.append(f"[cwd] {cwd}\n")
    lines.append(f"[command] {_format_command(command)}\n")
    if stdout_text:
        lines.append("[stdout]\n")
        lines.append(stdout_text)
        if not stdout_text.endswith("\n"):
            lines.append("\n")
    if stderr_text:
        lines.append("[stderr]\n")
        lines.append(stderr_text)
        if not stderr_text.endswith("\n"):
            lines.append("\n")
    lines.append("\n")
    payload = "".join(lines).encode("utf-8")
    fd = os.open(log_path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
    try:
        os.write(fd, payload)
    finally:
        os.close(fd)


def run_prep_command(command, args, section, cwd=None, stdout=None):
    """Run a prep external command and append stdout/stderr to prep_extend.log in debug mode."""
    debug = bool(getattr(args, "debug", False))
    stdout_target = subprocess.PIPE if debug and stdout is None else stdout
    if stdout_target is None:
        stdout_target = subprocess.DEVNULL

    try:
        result = subprocess.run(
            command,
            check=True,
            cwd=cwd,
            stdout=stdout_target,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as exc:
        if debug:
            _append_extend_log(
                args,
                section,
                command,
                cwd,
                stdout_text=_decode_stream(exc.stdout if stdout is None else ""),
                stderr_text=_decode_stream(exc.stderr),
            )
        raise

    if debug:
        _append_extend_log(
            args,
            section,
            command,
            cwd,
            stdout_text=_decode_stream(result.stdout if stdout is None else ""),
            stderr_text=_decode_stream(result.stderr),
        )
    return result

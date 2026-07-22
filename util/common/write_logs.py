"""Write shared console and file log messages."""

import time
import os

_LOG_FILE = None

def get_current_time():
    """Return a human-readable local timestamp."""
    formatted_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    return formatted_time


def set_log_file(path):
    """Enable optional append logging to a plain-text file."""
    global _LOG_FILE
    os.makedirs(os.path.dirname(path), exist_ok=True)
    _LOG_FILE = path


def log_message(level, message, bold=False, color="default", console=True):
    """Write a standardized status message to stdout and the active log file."""
    colors = {
        "default": "\033[0m",
        "success": "\033[1;32m",   # Green
        "step": "\033[1;36m",   # Cyan
        "warning": "\033[1;33m",# Yellow
        "error": "\033[1;31m" # Red
    }
    start_color = colors.get(color, colors["default"])
    end_color = "\033[0m"
    bold_code = "\033[1m" if bold else ""
    
    if console:
        print(f"{start_color}{level}{end_color} {bold_code}{message}{end_color}")
    if _LOG_FILE:
        with open(_LOG_FILE, "a", encoding="utf-8") as handle:
            handle.write(f"{get_current_time()}\t{level}\t{message}\n")

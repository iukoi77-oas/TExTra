"""Small helpers for custom command help screens."""

import os
import sys
import textwrap


def supports_color():
    """Return True when ANSI color is appropriate for help output."""
    return sys.stdout.isatty() and not os.environ.get("NO_COLOR")


def color_text(text, code, enabled):
    """Wrap text in an ANSI color code when enabled."""
    if not enabled:
        return text
    return f"\033[{code}m{text}\033[0m"


def help_title(text, *, color_enabled=None, title_color="95"):
    """Format a module help title with the TExTra title style."""
    if color_enabled is None:
        color_enabled = supports_color()
    return "\n" + color_text(text, f"1;{title_color}", color_enabled)


def help_box(title, rows, *, width=116, left_width=None, color_enabled=None, key_color="96", title_color="2"):
    """Format key/value rows in an ASCII help box."""
    if color_enabled is None:
        color_enabled = supports_color()
    rows = list(rows)
    title_plain = f" {title} "
    content_width = max(width, len(title_plain) + 4)
    left_width = left_width or max((len(left) for left, _right in rows), default=0)
    right_width = max(20, content_width - left_width - 6)
    top = "+-" + title_plain + "-" * max(0, content_width - len(title_plain)) + "-+"
    bottom = "+" + "-" * (content_width + 2) + "+"
    body = []
    for left, right in rows:
        wrapped = textwrap.wrap(str(right), width=right_width) or [""]
        for idx, part in enumerate(wrapped):
            label = left if idx == 0 else ""
            label_text = color_text(label, f"1;{key_color}", color_enabled) if label else ""
            plain = f"  {label:<{left_width}}  {part}"
            padding = " " * max(0, content_width - len(plain))
            rendered = f"  {label_text}{' ' * (left_width - len(label))}  {part}{padding}"
            body.append(f"|{rendered}|")
    return "\n".join([color_text(top, title_color, color_enabled), *body, color_text(bottom, title_color, color_enabled)])

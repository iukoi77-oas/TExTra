#!/usr/bin/env python
"""Top-level TExTra command dispatcher."""

import argparse
import os
import sys

package_dir = os.path.dirname(os.path.realpath(__file__))
project_root = os.path.dirname(package_dir)
sys.path.insert(0, os.path.abspath(project_root))

VERSION = "1.1.0"


def _supports_color():
    """Return True when ANSI color is appropriate for top-level help."""
    return sys.stdout.isatty() and not os.environ.get("NO_COLOR")


def _color(text, code, enabled):
    if not enabled:
        return text
    return f"\033[{code}m{text}\033[0m"


def _box(title, rows, color_enabled, *, key_color="96", title_color="2", width=96):
    """Format key/value rows in an ASCII help box."""
    rows = list(rows)
    content_width = max(
        width,
        len(title) + 4,
        *(len(left) + 2 + len(right) for left, right in rows),
    )
    top_label = f" {title} "
    top = "+-" + top_label + "-" * max(0, content_width - len(top_label)) + "-+"
    body = []
    left_width = max((len(left) for left, _right in rows), default=0)
    for left, right in rows:
        plain = f"  {left:<{left_width}}  {right}"
        padding = " " * max(0, content_width - len(plain))
        rendered_left = _color(left, f"1;{key_color}", color_enabled)
        rendered = f"  {rendered_left}{' ' * (left_width - len(left))}  {right}{padding}"
        body.append(f"|{rendered}|")
    bottom = "+" + "-" * (content_width + 2) + "+"
    return "\n".join(
        [
            _color(top, title_color, color_enabled),
            *body,
            _color(bottom, title_color, color_enabled),
        ]
    )


def _format_banner_line(line, enabled):
    """Color the leading TE block in the TExTra banner."""
    if not enabled:
        return line
    te_width = 16
    return _color(line[:te_width], "95", enabled) + _color(line[te_width:], "97", enabled)


def _format_top_level_help():
    """Build the custom top-level help screen."""
    color = _supports_color()
    cyan = "96"
    bold = "1"

    banner = [
        " _______ ______        _______        ",
        "|__   __|  ____|      |__   __|       ",
        "   | |  | |__   __  __   | |_ __ __ _ ",
        "   | |  |  __|  \\ \\/ /   | | '__/ _` |",
        "   | |  | |____  >  <    | | | | (_| |",
        "   |_|  |______|/_/\\_\\   |_|_|  \\__,_|",
    ]
    banner_text = "\n".join(_format_banner_line(line, color) for line in banner)
    lines = [
        banner_text,
        f"Version {VERSION}",
        "",
        "TExTra: TE-overlap exon discovery, quantification, and downstream analysis.",
        "",
        _box(
            "Usage",
            [
                ("TExTra <command> [options]", "Run a TExTra module"),
                ("TExTra <command> --help", "Show module-specific help"),
            ],
            color,
        ),
        "",
        _box(
            "Global options",
            [
                ("-h, --help", "Show this message and exit"),
                ("-v, --version", "Show version and exit"),
            ],
            color,
        ),
        "",
        _box(
            "Available commands",
            [
                ("prep", "Map reads, assemble transcripts, and convert annotations"),
                ("qual", "Identify TE-overlap exons and HITindex evidence"),
                ("quant", "Quantify transcripts and TE-overlap exon usage"),
                ("diff", "Test differential TE-overlap exon usage"),
                ("upstream", "Run prep + qual + quant in one command"),
                ("test", "Run prep + qual + quant + diff on a test dataset"),
            ],
            color,
        ),
        "",
        _box(
            "Examples",
            [
                ("TExTra prep --help", "Show prep module options"),
                ("TExTra upstream --help", "Show prep + qual + quant workflow options"),
                ("TExTra test --help", "Show bundled test workflow options"),
            ],
            color,
        ),
    ]
    return "\n".join(lines) + "\n"


def _dispatch(command):
    """Lazy-load command entrypoints to keep top-level CLI help lightweight."""
    if command == "prep":
        from TExTra.bin.mode0_pipeline import main as prep_main

        prep_main()
        return
    if command == "qual":
        from TExTra.bin.mode1_pipeline import main as qual_main

        qual_main()
        return
    if command == "quant":
        from TExTra.bin.mode2_pipeline import main as quant_main

        quant_main()
        return
    if command == "diff":
        from TExTra.bin.mode3_pipeline import main as diff_main

        diff_main()
        return
    if command == "upstream":
        from TExTra.bin.upstream_pipeline import main as upstream_main

        upstream_main()
        return
    if command == "test":
        from TExTra.bin.test_pipeline import main as test_main

        test_main()
        return


def main():
    """Route CLI subcommands to the corresponding pipeline entrypoint."""
    if len(sys.argv) == 1 or (len(sys.argv) > 1 and sys.argv[1] in {"-h", "--help"}):
        print(_format_top_level_help(), end="")
        return

    main_parser = argparse.ArgumentParser(
        prog="TExTra",
        description="TExTra: TE-overlap exon discovery, quantification, and downstream analysis.",
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    main_parser.add_argument("-v", "--version", action="version", version=f"TExTra {VERSION}")

    subparsers = main_parser.add_subparsers(
        title="Available Commands",
        dest="command",
        metavar="<command>",
        required=True,
    )

    subparsers.add_parser(
        "prep",
        add_help=False,
        help="Map reads, assemble transcripts, and convert annotations",
    )
    subparsers.add_parser(
        "qual",
        add_help=False,
        help="Identify TE-overlap exons and HITindex evidence",
    )
    subparsers.add_parser(
        "quant",
        add_help=False,
        help="Quantify transcripts and TE-overlap exon usage",
    )
    subparsers.add_parser(
        "diff",
        add_help=False,
        help="Test differential TE-overlap exon usage",
    )
    subparsers.add_parser(
        "upstream",
        add_help=False,
        help="Run prep + qual + quant in one command",
    )
    subparsers.add_parser(
        "test",
        add_help=False,
        help="Run prep + qual + quant + diff on a test dataset",
    )

    main_args, _ = main_parser.parse_known_args()

    _dispatch(main_args.command)

if __name__ == "__main__":
    main()

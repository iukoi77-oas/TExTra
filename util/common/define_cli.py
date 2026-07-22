"""Define shared CLI arguments and validation helpers."""


READTYPE_STRAND_CHOICES = {
    "paired": {"none", "rf", "fr"},
    "single": {"none", "r", "f"},
}


def _with_default(help_text, default_text):
    """Append a default sentence without duplicating punctuation."""
    return f"{str(help_text).rstrip('.')}. Default: {default_text}."


def add_threading_arguments(
    parser,
    threads_help="Number of threads to use. Default: 4.",
    njobs_help="Maximum number of parallel jobs; if omitted, use --threads.",
):
    """Add common thread and worker-count options."""
    parser.add_argument("-t", "--threads", type=int, default=4, help=threads_help)
    parser.add_argument(
        "--njobs",
        type=int,
        default=None,
        help=njobs_help,
    )


def add_read_layout_arguments(parser, strand_default="none", strand_help=None):
    """Add common readtype/strand options."""
    if strand_help is None:
        strand_help = (
            "Strand specificity: none, rf (reverse-forward), fr (forward-reverse), "
            "r (reverse, single-end), f (forward, single-end)"
        )
    parser.add_argument(
        "--strand",
        choices=["none", "rf", "fr", "r", "f"],
        default=strand_default,
        help=_with_default(strand_help, strand_default),
    )
    parser.add_argument(
        "--readtype",
        choices=["paired", "single"],
        default="paired",
        help="Read type: paired or single. Default: paired.",
    )


def add_project_argument(parser, default="project", help_text="Project prefix."):
    """Add common project prefix option."""
    parser.add_argument("--project", default=default, help=_with_default(help_text, default))


def add_run_mode_arguments(parser, debug_help=None):
    """Add common run-mode options.

    ``--debug`` keeps intermediate files and enables developer-oriented outputs.
    ``--detail`` is reserved for result-checking detail tables and summaries.
    """
    debug_help = debug_help or "Enable debug mode and keep intermediate files."
    parser.add_argument("--debug", action="store_true", help=f"{debug_help} Default: off.")
    parser.add_argument(
        "--detail",
        action="store_true",
        help="Enable detail mode for additional result-checking tables and summaries. Default: off.",
    )


def validate_threading_args(args):
    """Validate common thread and worker-count options."""
    if int(args.threads) < 1:
        raise RuntimeError("--threads must be a positive integer.")
    if getattr(args, "njobs", None) is not None and int(args.njobs) < 1:
        raise RuntimeError("--njobs must be a positive integer when provided.")


def validate_read_layout(args, module_name="module"):
    """Validate readtype/strand combinations shared by prep, qual, and quant."""
    allowed = READTYPE_STRAND_CHOICES.get(args.readtype, set())
    if args.strand not in allowed:
        allowed_text = ", ".join(sorted(allowed))
        raise RuntimeError(
            f"Invalid {module_name} read layout: --readtype {args.readtype} cannot be used with "
            f"--strand {args.strand}. Allowed --strand values: {allowed_text}."
        )

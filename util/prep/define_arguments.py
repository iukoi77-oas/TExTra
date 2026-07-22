"""Define prep CLI arguments and pass-through helpers."""


def add_prep_specific_arguments(parser, *, assembly_group=None, assignment_group=None, star_group=None):
    """Add arguments used by prep and by upstream when it invokes prep."""
    assembly = assembly_group or parser
    assignment = assignment_group or parser
    star = star_group or parser

    assembly.add_argument("--index", default=None, help="Existing STAR genome index directory. Default: unset; prep builds an index under the output directory for FASTQ input.")
    assembly.add_argument("--extend", default=None, help="Optional extra TE annotation file merged with --te before TE BED conversion. Default: unset.")
    assembly.add_argument(
        "--merge-method",
        choices=["taco", "stringtie"],
        default="taco",
        help="Transcript merge backend after per-replicate StringTie assembly: TACO or StringTie merge. Default: taco.",
    )
    assembly.add_argument("--taco-path", default=None, help="Path to taco_run or a TACO installation directory; used only with --merge-method taco. Default: auto-discover bundled/external TACO when available.")
    assembly.add_argument(
        "--assembly-mode",
        choices=["de-novo", "reference-guided"],
        default="de-novo",
        help="Per-replicate StringTie assembly mode: de-novo assembly or reference-guided assembly with --gene. Default: de-novo.",
    )
    assembly.add_argument("--optimized", action="store_true", help="Use optimized merge filtering defaults; explicit --min-expr/--min-length/--min-frac override the corresponding optimized values. Default: off.")
    assignment.add_argument("--gene-assignment", choices=["strict"], default="strict", help="Transcript-to-gene assignment policy after gffcompare. Default: strict.")
    assignment.add_argument("--include-overlap", action="store_true", help="Also retain lower-confidence reference-overlap transcript candidates during gene assignment. Default: off.")
    assembly.add_argument(
        "-j",
        "--junction",
        dest="junction",
        type=float,
        default=2.0,
        help="StringTie first-pass minimum junction coverage. Default: 2.0.",
    )
    assembly.add_argument(
        "--min-expr",
        type=float,
        default=None,
        help="Merge minimum expression (StringTie -F or TACO --filter-min-expr). If omitted, no expression filter is added unless --optimized is set.",
    )
    assembly.add_argument(
        "--min-length",
        type=int,
        default=None,
        help="Merge minimum transcript length (StringTie -m or TACO --filter-min-length). If omitted, no length filter is added unless --optimized is set.",
    )
    assembly.add_argument(
        "--min-frac",
        type=float,
        default=None,
        help="Merge minimum isoform fraction (StringTie -f or TACO --isoform-frac). If omitted, no isoform-fraction filter is added unless --optimized is set.",
    )
    star.add_argument("--sjdbOverhang", type=int, default=100, help="STAR sjdbOverhang used when prep builds a STAR index. Default: 100.")
    star.add_argument("--seedMultimapNmax", type=int, default=50000, help="STAR seedMultimapNmax passed to STAR alignment. Default: 50000.")
    star.add_argument("--winAnchorMultimapNmax", type=int, default=100, help="STAR winAnchorMultimapNmax passed to STAR alignment. Default: 100.")
    star.add_argument("--outFilterMultimapNmax", type=int, default=100, help="STAR outFilterMultimapNmax passed to STAR alignment. Default: 100.")


def build_prep_passthrough_args(args):
    """Build prep CLI arguments shared by upstream and standalone prep defaults."""
    prep_args = [
        "--assembly-mode", args.assembly_mode,
        "--merge-method", args.merge_method,
        "--gene-assignment", args.gene_assignment,
        "--junction", str(args.junction),
        "--sjdbOverhang", str(args.sjdbOverhang),
        "--seedMultimapNmax", str(args.seedMultimapNmax),
        "--winAnchorMultimapNmax", str(args.winAnchorMultimapNmax),
        "--outFilterMultimapNmax", str(args.outFilterMultimapNmax),
    ]

    optional_args = {
        "--min-expr": args.min_expr,
        "--min-length": args.min_length,
        "--min-frac": args.min_frac,
    }
    for flag, value in optional_args.items():
        if value is not None:
            prep_args.extend([flag, str(value)])

    optional_paths = {
        "--index": args.index,
        "--extend": args.extend,
        "--taco-path": args.taco_path,
    }
    for flag, value in optional_paths.items():
        if value:
            prep_args.extend([flag, value])

    optional_flags = [
        ("--debug", args.debug),
        ("--detail", getattr(args, "detail", False)),
        ("--optimized", args.optimized),
        ("--include-overlap", args.include_overlap),
    ]
    for flag, enabled in optional_flags:
        if enabled:
            prep_args.append(flag)

    return prep_args

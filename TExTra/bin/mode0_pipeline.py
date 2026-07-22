"""CLI entrypoint for TExTra prep: input parsing, alignment, assembly, and conversion."""

import sys
import os
import argparse
import json
import textwrap
import pandas as pd

# Allow direct execution from the source tree without requiring installation.
script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.insert(0, os.path.abspath(parent_dir))

from util.common.write_logs import log_message, set_log_file
from util.common.run_config import update_global_config
from util.prep.define_arguments import add_prep_specific_arguments
from util.common.define_cli import (
    add_read_layout_arguments,
    add_run_mode_arguments,
    add_threading_arguments,
    validate_read_layout,
    validate_threading_args,
)
from util.common.help_format import help_title


class RawDefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Preserve epilog formatting while still showing default values."""


def _supports_color():
    return sys.stdout.isatty() and not os.environ.get("NO_COLOR")


def _color(text, code, enabled):
    if not enabled:
        return text
    return f"\033[{code}m{text}\033[0m"


def _help_box(title, rows, *, width=116, left_width=None):
    color = _supports_color()
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
            label_text = _color(label, "1;96", color) if label else ""
            plain = f"  {label:<{left_width}}  {part}"
            padding = " " * max(0, content_width - len(plain))
            rendered = f"  {label_text}{' ' * (left_width - len(label))}  {part}{padding}"
            body.append(f"|{rendered}|")
    return "\n".join([_color(top, "2", color), *body, _color(bottom, "2", color)])


class PrepArgumentParser(argparse.ArgumentParser):
    """argparse parser with custom top-level prep help layout."""

    def format_help(self):
        lines = [
            help_title("TExTra prep"),
            "Map reads, assemble consensus transcripts, assign transcripts to genes, and convert gene/TE annotations.",
            "",
            _help_box(
                "Usage",
                [
                    ("TExTra prep -i INPUT -g GENOME -G GENE -o OUT_DIR -r TE [options]", "Run prep"),
                    ("TExTra prep --help", "Show prep help"),
                ],
            ),
            "",
            _help_box(
                "Required arguments",
                [
                    ("-i, --input", "Input TSV. First column is sample/condition; following columns are replicate inputs."),
                    ("-g, --genome", "Genome FASTA used for STAR index generation and annotation conversion."),
                    ("-G, --gene", "Reference gene annotation GTF used by StringTie/gffcompare and converted to gene BED."),
                    ("-o, --out_dir", "Prep output root directory."),
                    ("-r, --te", "TE annotation file converted to TE BED; GTF/BED/RepeatMasker-style inputs are supported."),
                ],
            ),
            "",
            _help_box(
                "Optional arguments",
                [
                    ("-h, --help", "Show this message and exit."),
                    ("-t, --threads", "Number of threads to use. Default: 4."),
                    ("--njobs", "Maximum number of parallel jobs. Default: omitted, use --threads."),
                    ("--strand", "Strand specificity: none, rf, fr, r, f. Default: none."),
                    ("--readtype", "Read type: paired or single. Default: paired."),
                    ("--debug", "Enable debug mode and keep intermediate files. Default: off."),
                    ("--detail", "Enable detail mode for additional result-checking tables and summaries. Default: off."),
                    ("--assembly-mode", "Per-replicate StringTie assembly mode: de-novo or reference-guided. Default: de-novo."),
                    ("--merge-method", "Transcript merge backend: taco or stringtie. Default: taco."),
                    ("--optimized", "Use optimized merge filtering defaults; explicit --min-expr/--min-length/--min-frac override corresponding values. Default: off."),
                ],
            ),
            "",
            _help_box(
                "Advanced options",
                [
                    ("--index", "Existing STAR genome index directory. Default: unset; prep builds one for FASTQ input."),
                    ("--extend", "Optional extra TE annotation file merged with --te before TE BED conversion. Default: unset."),
                    ("--taco-path", "Path to taco_run or a TACO installation directory; used only with --merge-method taco. Default: auto-discover bundled/external TACO."),
                    ("--gene-assignment", "Transcript-to-gene assignment policy after gffcompare. Default: strict."),
                    ("--include-overlap", "Also retain lower-confidence reference-overlap transcript candidates during gene assignment. Default: off."),
                    ("-j, --junction", "StringTie first-pass minimum junction coverage. Default: 2.0."),
                    ("--min-expr", "Merge minimum expression; if omitted, no expression filter is added unless --optimized is set."),
                    ("--min-length", "Merge minimum transcript length; if omitted, no length filter is added unless --optimized is set."),
                    ("--min-frac", "Merge minimum isoform fraction; if omitted, no isoform-fraction filter is added unless --optimized is set."),
                    ("--sjdbOverhang", "STAR sjdbOverhang used when prep builds a STAR index. Default: 100."),
                    ("--seedMultimapNmax", "STAR seedMultimapNmax passed to STAR alignment. Default: 50000."),
                    ("--winAnchorMultimapNmax", "STAR winAnchorMultimapNmax passed to STAR alignment. Default: 100."),
                    ("--outFilterMultimapNmax", "STAR outFilterMultimapNmax passed to STAR alignment. Default: 100."),
                ],
            ),
            "",
            _help_box(
                "Examples",
                [
                    ("TExTra prep --input input.tsv --out_dir result --genome genome.fa --gene gene.gtf --te TE_annotation", "Minimal run"),
                    ("TExTra prep -i input.tsv -o result -g genome.fa -G gene.gtf -r TE_annotation --strand rf", "Stranded run"),
                ],
            ),
        ]
        return "\n".join(lines) + "\n"


def parse_arguments(args_list):
    parser = PrepArgumentParser(
        prog="TExTra prep",
        description="Map reads, assemble consensus transcripts, assign transcripts to genes, and convert gene/TE annotations for downstream TExTra modules.",
        add_help=True,
        formatter_class=RawDefaultsHelpFormatter,
        epilog=(
            "Examples:\n"
            "  TExTra prep --input input.tsv --out_dir result --genome genome.fa --gene gene.gtf --te TE_annotation\n"
            "  TExTra prep -i input.tsv -o result -g genome.fa -G gene.gtf -r TE_annotation --strand rf --readtype paired"
        ),
    )

    required = parser.add_argument_group("Required arguments")
    execution = parser.add_argument_group("Execution options")
    read_layout = parser.add_argument_group("Read layout")
    run_modes = parser.add_argument_group("Run modes")
    assembly = parser.add_argument_group("Assembly and merge options")
    assignment = parser.add_argument_group("Gene assignment options")
    star = parser.add_argument_group("STAR native options")

    required.add_argument(
        "-i",
        "--input",
        required=True,
        help=(
            "Input TSV. First column is sample/condition name; following columns are replicate inputs. "
            "Use read1,read2 within one cell for paired FASTQ; use one BAM or one FASTQ per cell otherwise."
        ),
    )
    required.add_argument("-g", "--genome", required=True, help="Genome FASTA used for STAR index generation and annotation conversion.")
    required.add_argument("-G", "--gene", required=True, help="Reference gene annotation GTF used by StringTie/gffcompare and converted to gene BED.")
    required.add_argument("-o", "--out_dir", required=True, help="Prep output root directory.")
    required.add_argument("-r", "--te", required=True, help="TE annotation file converted to TE BED; GTF/BED/RepeatMasker-style inputs are supported.")
    add_threading_arguments(execution)
    add_read_layout_arguments(read_layout)
    add_run_mode_arguments(run_modes)
    add_prep_specific_arguments(
        parser,
        assembly_group=assembly,
        assignment_group=assignment,
        star_group=star,
    )

    return parser.parse_args(args_list)


def _build_input_info(input_path):
    """Convert input table into canonical mate1/mate2/sample records."""
    input_df = pd.read_csv(input_path, sep="\t", header=None)
    input_df.columns = ["samples"] + [f"rep{i+1}" for i in range(input_df.shape[1] - 1)]
    input_df = input_df.replace(["NA", "NaN", "nan", ""], pd.NA)

    mate1_list = []
    mate2_list = []
    group_list = []
    seen_input_samples = set()
    for row in input_df.itertuples(index=False, name=None):
        raw_sample = row[0]
        if pd.isna(raw_sample):
            raise RuntimeError("Input TSV contains an empty sample name.")
        sample = str(raw_sample).strip()
        if not sample or sample in {"NA", "NaN", "nan"}:
            raise RuntimeError("Input TSV contains an empty sample name.")
        if sample in seen_input_samples:
            raise RuntimeError(
                f"Duplicate sample '{sample}' found in input TSV. Put all replicates for one sample on the same row."
            )
        seen_input_samples.add(sample)
        rep_idx = 0
        for file_path in row[1:]:
            if pd.isna(file_path):
                continue
            rep_idx += 1
            if "," in str(file_path):
                mate1, mate2 = [token.strip() for token in str(file_path).split(",", 1)]
            else:
                mate1, mate2 = str(file_path).strip(), "-"
            run_id = f"{sample}_rep{rep_idx}"
            mate1_list.append(mate1)
            mate2_list.append(mate2)
            group_list.append(run_id)

    return pd.DataFrame({"mate1": mate1_list, "mate2": mate2_list, "sample": group_list})


def _validate_input_table(input_info, args):
    """Validate resolved replicate rows from the user input table."""
    if input_info.empty:
        raise RuntimeError("Input table produced no valid replicate rows.")

    seen_samples = set()
    for row in input_info.itertuples(index=False):
        sample = str(row.sample).strip()
        mate1 = str(row.mate1).strip()
        mate2 = str(row.mate2).strip()
        if not sample:
            raise RuntimeError("Encountered an empty sample/replicate identifier in input table.")
        if sample in seen_samples:
            raise RuntimeError(f"Duplicate resolved replicate name detected: {sample}")
        seen_samples.add(sample)

        mate1_is_bam = mate1.lower().endswith(".bam")
        mate2_is_bam = mate2.lower().endswith(".bam")
        if not os.path.isfile(mate1):
            raise FileNotFoundError(f"Input file not found for {sample}: {mate1}")
        if mate2 not in {"", "-"} and not os.path.isfile(mate2):
            raise FileNotFoundError(f"Input mate2 file not found for {sample}: {mate2}")

        if mate1_is_bam:
            if mate2 not in {"", "-"}:
                raise RuntimeError(
                    f"BAM input must be provided as a single file per replicate ({sample}); got: {mate1}, {mate2}"
                )
            continue

        if str(args.readtype).lower() == "paired" and mate2 in {"", "-"}:
            raise RuntimeError(
                f"Paired-end mode requires mate2 for FASTQ replicate {sample}. Use mate1,mate2 in the input TSV."
            )
        if str(args.readtype).lower() == "single" and mate2 not in {"", "-"}:
            raise RuntimeError(
                f"Single-end mode requires exactly one FASTQ per replicate ({sample}); got mate2: {mate2}"
            )
        if mate2_is_bam:
            raise RuntimeError(
                f"Mixed FASTQ/BAM input is not allowed within one replicate ({sample}): {mate1}, {mate2}"
            )
        if mate2 not in {"", "-"}:
            mate1_gz = mate1.lower().endswith((".fq.gz", ".fastq.gz"))
            mate2_gz = mate2.lower().endswith((".fq.gz", ".fastq.gz"))
            if mate1_gz != mate2_gz:
                raise RuntimeError(
                    f"Mixed compressed/uncompressed FASTQ input is not supported within one replicate ({sample}): {mate1}, {mate2}"
                )


def _validate_inputs(args):
    """Validate required user-provided input paths."""
    required_files = {
        "input table": args.input,
        "genome FASTA": args.genome,
        "gene annotation": args.gene,
        "TE annotation": args.te,
    }
    for label, path in required_files.items():
        if not os.path.isfile(path):
            raise FileNotFoundError(f"{label} not found: {path}")
    if args.extend and not os.path.isfile(args.extend):
        raise FileNotFoundError(f"extended TE annotation not found: {args.extend}")


def _validate_parameters(args):
    """Validate prep CLI parameters before running external tools."""
    validate_threading_args(args)
    validate_read_layout(args, module_name="prep")
    if int(args.sjdbOverhang) < 1:
        raise RuntimeError("--sjdbOverhang must be a positive integer.")
    for name in ["seedMultimapNmax", "winAnchorMultimapNmax", "outFilterMultimapNmax"]:
        if int(getattr(args, name)) < 1:
            raise RuntimeError(f"--{name} must be a positive integer.")
    if args.min_expr is not None and float(args.min_expr) < 0:
        raise RuntimeError("--min-expr must be non-negative when provided.")
    if args.min_length is not None and int(args.min_length) < 1:
        raise RuntimeError("--min-length must be a positive integer when provided.")
    if args.min_frac is not None and not (0.0 <= float(args.min_frac) <= 1.0):
        raise RuntimeError("--min-frac must be between 0 and 1 when provided.")


def _on_off(value):
    return "on" if bool(value) else "off"


def _merge_filter_value(args, attr_name):
    value = getattr(args, attr_name, None)
    if value is not None:
        return str(value)
    return "optimized" if getattr(args, "optimized", False) else "not_set"


def _log_prep_parameters(args, input_info):
    """Print the prep mode choices that change outputs or downstream context."""
    mate1_values = input_info["mate1"].fillna("").astype(str).str.strip().str.lower()
    bam_n = int(mate1_values.str.endswith(".bam").sum())
    fastq_n = int(len(input_info) - bam_n)
    if fastq_n == 0:
        star_index = "not_used"
    else:
        star_index = "provided" if args.index else "auto"
    merge_method = "StringTie" if args.merge_method == "stringtie" else "TACO"
    log_message(
        "[INFO]",
        (
            "Prep parameters: "
            f"replicates={len(input_info)}, FASTQ={fastq_n}, BAM={bam_n}, "
            f"readtype={args.readtype}, strand={args.strand}, "
            f"assembly_mode={args.assembly_mode}, "
            f"merge={merge_method}, optimized={_on_off(args.optimized)}, "
            f"junction={args.junction}, "
            f"min_expr={_merge_filter_value(args, 'min_expr')}, "
            f"min_length={_merge_filter_value(args, 'min_length')}, "
            f"min_frac={_merge_filter_value(args, 'min_frac')}, "
            f"STAR index={star_index}"
        ),
        color="info",
    )
    if getattr(args, "detail", False) or getattr(args, "debug", False):
        log_message(
            "[INFO]",
            (
                "Prep STAR parameters: "
                f"sjdbOverhang={args.sjdbOverhang}, "
                f"seedMultimapNmax={args.seedMultimapNmax}, "
                f"winAnchorMultimapNmax={args.winAnchorMultimapNmax}, "
                f"outFilterMultimapNmax={args.outFilterMultimapNmax}"
            ),
            color="info",
        )


def _prepare_debug_logs(args):
    """Prepare prep debug log files and remove legacy debug artifacts."""
    if not args.debug:
        return
    log_dir = os.path.join(args.out_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    for name in ["prep_extend.log", "extend.log", "prep.parameters.json", "prep.run_config.json"]:
        path = os.path.join(log_dir, name)
        if os.path.exists(path):
            os.remove(path)


def _write_run_config(args, input_info):
    """Write minimal prep context that downstream modules may reuse."""
    replicates = {}
    for replicate in input_info["sample"].astype(str).tolist():
        sample = replicate
        if "_rep" in replicate:
            prefix, suffix = replicate.rsplit("_rep", 1)
            if prefix and suffix.isdigit():
                sample = prefix
        replicates.setdefault(sample, []).append(replicate)
    module_payload = {
        "out_dir": os.path.abspath(args.out_dir),
        "strand": args.strand,
        "readtype": args.readtype,
        "assembly_mode": args.assembly_mode,
        "merge_method": args.merge_method,
        "optimized": bool(args.optimized),
        "junction": float(args.junction),
        "min_expr": args.min_expr if args.min_expr is None else float(args.min_expr),
        "min_length": args.min_length if args.min_length is None else int(args.min_length),
        "min_frac": args.min_frac if args.min_frac is None else float(args.min_frac),
        "samples": sorted(replicates),
        "replicates": {sample: sorted(values) for sample, values in sorted(replicates.items())},
        "detail": bool(getattr(args, "detail", False)),
        "debug": bool(getattr(args, "debug", False)),
    }
    update_global_config(args.out_dir, "prep", module_payload)
    if args.debug:
        log_dir = os.path.join(args.out_dir, "logs")
        os.makedirs(log_dir, exist_ok=True)
        debug_payload = {"module": "prep", "schema_version": 1, "run_mode": "debug", **module_payload}
        with open(os.path.join(log_dir, "prep_config.json"), "w", encoding="utf-8") as handle:
            json.dump(debug_payload, handle, indent=2, sort_keys=True)


def main(args_list=None):
    if args_list is None:
        args_list = sys.argv[2:]
    args = parse_arguments(args_list)
    if args.debug:
        args.detail = True
    from TExTra.src.mode0.step1_alignment import align_func
    from TExTra.src.mode0.step2_assembly import assemble_func
    from TExTra.src.mode0.step3_teconvert import convert_func

    _validate_parameters(args)
    _validate_inputs(args)
    set_log_file(os.path.join(args.out_dir, "logs", "prep.log"))
    _prepare_debug_logs(args)

    input_info = _build_input_info(args.input)
    _validate_input_table(input_info, args)
    _write_run_config(args, input_info)
    _log_prep_parameters(args, input_info)

    bamfiles, samples = align_func(input_info, args)
    consensus_gtf = assemble_func(bamfiles, samples, args)
    convert_func(consensus_gtf, args)


if __name__ == "__main__":
    main()

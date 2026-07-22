"""CLI entrypoint for TExTra quant: transcript quantification and TE-exon usage."""

import json
import os
import sys
import argparse

# Allow direct execution from the source tree without requiring installation.
script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
project_root = os.path.dirname(parent_dir)
sys.path.insert(0, os.path.abspath(project_root))
sys.path.insert(0, os.path.abspath(parent_dir))

from util.common.define_cli import (
    add_project_argument,
    add_read_layout_arguments,
    add_run_mode_arguments,
    add_threading_arguments,
    validate_read_layout,
    validate_threading_args,
)
from util.common.define_layout import CLASSIFY_DIR, ALIGNMENT_DIR, ASSEMBLY_DIR, TRANSCRIPT_GENE_ASSIGNMENT_TSV, resolve_consensus_gtf
from util.common.run_config import extract_module_config, global_config_path, read_global_config, update_global_config
from util.common.help_format import help_box, help_title
from util.common.write_logs import log_message, set_log_file


def _option_was_provided(args_list, *names):
    tokens = list(args_list or [])
    for token in tokens:
        for name in names:
            if token == name or token.startswith(name + "="):
                return True
    return False


class QuantArgumentParser(argparse.ArgumentParser):
    """argparse parser with custom quant help layout."""

    def format_help(self):
        lines = [
            help_title("TExTra quant"),
            "Quantify consensus transcripts and project abundance to TE-overlap exon usage.",
            "",
            help_box(
                "Usage",
                [
                    ("TExTra quant --prep PREP --qual QUAL -o OUT_DIR [options]", "Run quant from prep/qual outputs"),
                    ("TExTra quant --reuse -o OUT_DIR [options]", "Infer prep/qual context when possible"),
                ],
            ),
            "",
            help_box(
                "Required arguments",
                [
                    ("-o, --out_dir", "Quant output root directory."),
                    ("--prep", "Path to TExTra prep output; required unless inferred by --reuse."),
                    ("--qual", "Path to TExTra qual output; required unless inferred by --reuse."),
                    ("-s, --samples", "Sample/condition names, comma-separated; required unless inferred by --reuse."),
                ],
            ),
            "",
            help_box(
                "Optional arguments",
                [
                    ("-h, --help", "Show this message and exit."),
                    ("-t, --threads", "Number of threads to use. Default: 4."),
                    ("--njobs", "Maximum number of parallel jobs. Default: omitted, use --threads."),
                    ("-i, --input", "Original prep input TSV used to resolve quant reads. Default: unset; BAM fallback may be used."),
                    ("-g, --genome", "Genome FASTA used to build RSEM reference or Salmon transcript FASTA when needed. Default: infer from config when possible."),
                    ("--quantifier", "Quantification backend: rsem or salmon. Default: rsem."),
                    ("--project", "Project name. Default: project."),
                    ("--strand", "Strand-specific RNA-seq library type. Default: none."),
                    ("--readtype", "Read type: paired or single. Default: paired."),
                    ("--compute-gene-abundance", "Compute and export gene-level abundance tables. Default: off."),
                    ("--reuse", "Infer omitted prep/qual context; does not reuse quantification results. Default: off."),
                    ("--debug", "Enable debug mode and keep intermediate files. Default: off."),
                    ("--detail", "Enable detail mode for additional result-checking tables and summaries. Default: off."),
                ],
            ),
            "",
            help_box(
                "Advanced options",
                [
                    ("--quant-result-dir", "Directory containing reusable RSEM/Salmon backend outputs. Default: unset; compute backend quantification as needed."),
                ],
            ),
            "",
            help_box(
                "Examples",
                [
                    ("TExTra quant --prep result --qual result -o result -s sample1,sample2 -g genome.fa", "Run quant"),
                    ("TExTra quant --reuse -o result -g genome.fa", "Infer context from config"),
                ],
            ),
        ]
        return "\n".join(lines) + "\n"


def parse_arguments(args_list):
    """Parse quant CLI arguments."""
    parser = QuantArgumentParser(
        prog="TExTra quant",
        description="TExTra quantitative analysis on consensus transcripts using RSEM or Salmon.",
        add_help=True,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("-o", "--out_dir", type=str, required=True, help="Output root directory for quant results.")

    add_threading_arguments(parser)
    parser.add_argument(
        "-s",
        "--samples",
        default=None,
        help="Sample/condition names, comma-separated; inferred with --reuse when omitted if possible.",
    )

    add_read_layout_arguments(parser, strand_help="Strand-specific RNA-seq library type.")

    parser.add_argument("--prep", default=None, help="Path to TExTra prep output; required unless inferred by --reuse.")
    parser.add_argument("--qual", default=None, help="Path to TExTra qual output; required unless inferred by --reuse.")
    parser.add_argument("-i", "--input", default=None, help="Original prep input TSV used to resolve quant reads.")
    parser.add_argument("-g", "--genome", type=str, default=None, help="Genome FASTA used to build RSEM reference or Salmon transcript FASTA.")
    parser.add_argument("--quantifier", choices=["rsem", "salmon"], default="rsem", help="Quantification backend.")
    parser.add_argument(
        "--quant-result-dir",
        default=None,
        help=(
            "Directory containing reusable RSEM/Salmon quantification outputs "
            "(typically 04_quantification from a debug quant run). "
            "If omitted, backend quantification is computed as needed; final usage tables are not reused."
        ),
    )
    add_project_argument(parser, default="project", help_text="Project name")
    add_run_mode_arguments(parser)
    parser.add_argument(
        "--compute-gene-abundance",
        action="store_true",
        default=False,
        help="Compute and export gene-level abundance tables.",
    )
    parser.add_argument(
        "--reuse",
        action="store_true",
        help=(
            "Infer omitted --prep, --qual, --samples, --strand, and --readtype from run context; "
            "does not reuse quantification results."
        ),
    )

    args = parser.parse_args(args_list)
    internal_defaults = {
        "show_tool_output": False,
        "force": False,
        "rsem_reference": None,
        "rsem_build_reference": False,
        "rsem_read1": None,
        "rsem_read2": None,
        "rsem_single_reads": None,
        "rsem_no_qualities": False,
        "rsem_aligner": "auto",
        "rsem_star_path": None,
        "rsem_strandedness": "auto",
        "salmon_index": None,
        "salmon_build_index": False,
        "salmon_libtype": "auto",
        "rsem_auto_bam2fastq": True,
        "rsem_bam2fastq_dir": None,
    }
    for name, value in internal_defaults.items():
        setattr(args, name, value)
    args._explicit_prep = _option_was_provided(args_list, "--prep")
    args._explicit_qual = _option_was_provided(args_list, "--qual")
    args._explicit_samples = _option_was_provided(args_list, "-s", "--samples")
    args._explicit_strand = _option_was_provided(args_list, "--strand")
    args._explicit_readtype = _option_was_provided(args_list, "--readtype")
    return args


def _validate_parameters(args):
    """Validate shared quant CLI parameters."""
    validate_threading_args(args)
    validate_read_layout(args, module_name="quant")


def _is_prep_output_dir(path):
    if not path or not os.path.isdir(path):
        return False
    return (
        os.path.isfile(resolve_consensus_gtf(path))
        and os.path.isfile(os.path.join(path, ASSEMBLY_DIR, TRANSCRIPT_GENE_ASSIGNMENT_TSV))
    )


def _is_qual_output_dir(path):
    if not path or not os.path.isdir(path):
        return False
    return os.path.isfile(os.path.join(path, CLASSIFY_DIR, "transcript_exon_te_annotation.tsv"))


def _infer_samples_from_alignment(prep_dir):
    alignment_dir = os.path.join(prep_dir, ALIGNMENT_DIR)
    if not os.path.isdir(alignment_dir):
        return []
    samples = set()
    for name in os.listdir(alignment_dir):
        path = os.path.join(alignment_dir, name)
        if not os.path.isdir(path):
            continue
        if "_rep" in name:
            prefix, suffix = name.rsplit("_rep", 1)
            if prefix and suffix.isdigit():
                samples.add(prefix)
    return sorted(samples)


def _config_candidates(out_dir):
    candidates = []
    for path in [out_dir, os.getcwd()]:
        if not path:
            continue
        abs_path = os.path.abspath(path)
        if abs_path in candidates:
            continue
        candidates.append(abs_path)
    return candidates


def _discover_reuse_context(out_dir):
    prep_candidates = []
    qual_candidates = []
    prep_config = None
    qual_config = None
    used_config_path = None

    for root in _config_candidates(out_dir):
        config_path = global_config_path(root)
        if os.path.isfile(config_path):
            config, path = read_global_config(root)
            used_config_path = path
            prep_cfg = extract_module_config(config, "prep")
            qual_cfg = extract_module_config(config, "qual")
            if prep_cfg and prep_config is None:
                prep_config = prep_cfg
            if qual_cfg and qual_config is None:
                qual_config = qual_cfg
            if prep_cfg:
                prep_dir = os.path.abspath(prep_cfg.get("out_dir", root))
                if _is_prep_output_dir(prep_dir) and prep_dir not in prep_candidates:
                    prep_candidates.append(prep_dir)
            if qual_cfg:
                qual_dir = os.path.abspath(qual_cfg.get("out_dir", root))
                if _is_qual_output_dir(qual_dir) and qual_dir not in qual_candidates:
                    qual_candidates.append(qual_dir)
        if _is_prep_output_dir(root) and root not in prep_candidates:
            prep_candidates.append(root)
        if _is_qual_output_dir(root) and root not in qual_candidates:
            qual_candidates.append(root)

    return prep_candidates, qual_candidates, prep_config, qual_config, used_config_path


def _single_or_error(candidates, label):
    if len(candidates) > 1:
        raise RuntimeError(
            f"--reuse found multiple possible {label} output directories. "
            f"Please pass --{label} explicitly: " + ", ".join(candidates)
        )
    return candidates[0] if candidates else None


def _apply_reuse_context(args):
    if not args.reuse:
        missing = []
        if not args.prep:
            missing.append("--prep")
        if not args.qual:
            missing.append("--qual")
        if not args.samples:
            missing.append("--samples")
        if missing:
            raise RuntimeError(
                "Missing required quant argument(s): "
                + ", ".join(missing)
                + ". Provide them explicitly, or rerun with --reuse to infer them from run config."
            )
        return args

    prep_candidates, qual_candidates, prep_config, qual_config, config_path = _discover_reuse_context(args.out_dir)
    reused = []

    if not args.prep:
        args.prep = _single_or_error(prep_candidates, "prep")
        if args.prep:
            reused.append(f"prep={args.prep}")
    if not args.qual:
        args.qual = _single_or_error(qual_candidates, "qual")
        if args.qual:
            reused.append(f"qual={args.qual}")
    if not args.prep:
        raise RuntimeError("--reuse could not locate a prep output directory. Pass --prep explicitly.")
    if not args.qual:
        raise RuntimeError("--reuse could not locate a qual output directory. Pass --qual explicitly.")

    context = qual_config or prep_config or {}
    if not args.samples:
        samples = []
        if qual_config:
            samples = [str(sample) for sample in qual_config.get("samples", []) if str(sample).strip()]
        if not samples and prep_config:
            samples = [str(sample) for sample in prep_config.get("samples", []) if str(sample).strip()]
        if not samples:
            samples = _infer_samples_from_alignment(args.prep)
        if samples:
            args.samples = ",".join(samples)
            reused.append(f"samples={len(samples)}")
    if not args.samples:
        raise RuntimeError("--reuse could not infer --samples from run context. Pass --samples explicitly.")

    if context and not args._explicit_strand and context.get("strand"):
        args.strand = str(context["strand"])
        reused.append(f"strand={args.strand}")
    if context and not args._explicit_readtype and context.get("readtype"):
        args.readtype = str(context["readtype"])
        reused.append(f"readtype={args.readtype}")

    if not (prep_config or qual_config):
        log_message(
            "[WARNING]",
            f"--reuse did not find module run config in {config_path or global_config_path(args.out_dir)}; file-based inference was used.",
            color="warning",
        )
    if reused:
        log_message("[INFO]", "Reused quant context: " + ", ".join(reused), color="info")
    return args


def _abs_or_none(path):
    return os.path.abspath(path) if path else None


def _split_samples(samples):
    return [item.strip() for item in str(samples or "").split(",") if item.strip()]


def _write_quant_run_config(args):
    payload = {
        "out_dir": os.path.abspath(args.out_dir),
        "prep": _abs_or_none(args.prep),
        "qual": _abs_or_none(args.qual),
        "input": _abs_or_none(getattr(args, "input", None)),
        "genome": _abs_or_none(getattr(args, "genome", None)),
        "project": args.project,
        "samples": _split_samples(args.samples),
        "strand": args.strand,
        "readtype": args.readtype,
        "threads": int(args.threads),
        "njobs": None if args.njobs is None else int(args.njobs),
        "quantifier": str(getattr(args, "quantifier", "rsem")).lower(),
        "compute_gene_abundance": bool(getattr(args, "compute_gene_abundance", False)),
        "reuse": bool(getattr(args, "reuse", False)),
        "quant_result_dir": _abs_or_none(getattr(args, "quant_result_dir", None)),
        "detail": bool(getattr(args, "detail", False)),
        "debug": bool(getattr(args, "debug", False)),
    }
    update_global_config(args.out_dir, "quant", payload)
    if getattr(args, "debug", False):
        logs_dir = os.path.join(args.out_dir, "logs")
        os.makedirs(logs_dir, exist_ok=True)
        debug_payload = {"module": "quant", "schema_version": 1, "run_mode": "debug", **payload}
        with open(os.path.join(logs_dir, "quant_config.json"), "w", encoding="utf-8") as handle:
            json.dump(debug_payload, handle, indent=2, sort_keys=True)
    return payload


def main(args_list=None):
    """Run the quant CLI entrypoint."""
    if args_list is None:
        args_list = sys.argv[2:] if len(sys.argv) > 1 and sys.argv[1] == "quant" else sys.argv[1:]
    args = parse_arguments(args_list)
    if args.debug:
        args.detail = True
    set_log_file(os.path.join(args.out_dir, "logs", "quant.log"))
    args = _apply_reuse_context(args)
    _validate_parameters(args)
    _write_quant_run_config(args)

    from TExTra.src.mode2.step1_count import count_func

    count_func(args)


if __name__ == "__main__":
    main()

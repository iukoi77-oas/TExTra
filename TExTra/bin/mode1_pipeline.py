"""CLI entrypoint for TExTra qual: TE-exon annotation and HITindex classification."""

import argparse
import json
import os
import sys

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
)
from util.common.define_layout import ALIGNMENT_DIR, ASSEMBLY_DIR, CONSENSUS_GTF, CONVERT_DIR
from util.common.run_config import (
    extract_module_config,
    global_config_path,
    read_global_config,
)
from util.common.help_format import help_box, help_title
from util.common.write_logs import log_message, set_log_file
from TExTra.src.mode1.step1_identify import identify_func


def _option_was_provided(args_list, *names):
    tokens = list(args_list or [])
    for token in tokens:
        for name in names:
            if token == name or token.startswith(name + "="):
                return True
    return False


class QualArgumentParser(argparse.ArgumentParser):
    """argparse parser with custom qual help layout."""

    def format_help(self):
        lines = [
            help_title("TExTra qual"),
            "Identify TE-overlap exons and optional HITindex positional evidence from prep outputs.",
            "",
            help_box(
                "Usage",
                [
                    ("TExTra qual --prep PREP -o OUT_DIR [options]", "Run qual from prep output"),
                    ("TExTra qual --reuse -o OUT_DIR [options]", "Infer prep context when possible"),
                ],
            ),
            "",
            help_box(
                "Required arguments",
                [
                    ("-o, --out_dir", "Qual output root directory."),
                    ("--prep", "Path to TExTra prep output; required unless inferred by --reuse."),
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
                    ("--project", "Project prefix for classify combined outputs. Default: project."),
                    ("--strand", "Strand specificity: none, rf, fr, r, f. Default: none."),
                    ("--readtype", "Read type: paired or single. Default: paired."),
                    ("--calculate-afe-ale", "Calculate AFE/ALE exon usage outputs. Requires HITindex. Default: off."),
                    ("--ignore-junction", "Ignore TE-overlap junction support/degradation checks. Default: off; cannot be combined with explicit --junction."),
                    ("-j, --junction", "Minimum TE-side junction reads for TE-overlap degradation. Default: 2.0; inherits modified prep --junction when omitted."),
                    ("--reuse", "Infer omitted prep context; does not reuse qual result files. Default: off."),
                    ("--debug", "Enable debug mode and keep intermediate files. Default: off."),
                    ("--detail", "Enable detail mode for additional result-checking tables and summaries. Default: off."),
                ],
            ),
            "",
            help_box(
                "Advanced options",
                [
                    ("--ss3-buffer", "3' splice-site buffer size for HITindex splice-site roles. Default: 20."),
                    ("--ss5-buffer", "5' splice-site buffer size for HITindex splice-site roles. Default: 50."),
                    ("--genmodel-iters", "ADVI iterations in HITindex generative model. Default: 100000."),
                    ("--bootstrap-n", "Bootstrap iterations for HITindex significance. Default: 1000."),
                    ("--seed", "Optional random seed for HITindex model fitting and bootstrap. Default: unset."),
                    ("--te-overlap-min-bp", "Minimum overlap bp for TE-overlap metrics. Default: 10."),
                    ("--te-overlap-min-frac", "Minimum overlap fraction for TE-overlap metrics. Default: 0.1."),
                    ("--splice-site-flank-bp", "Flank window around splice sites for boundary/anchor/var-site hit checks. Default: 10."),
                    ("--skip-hitindex", "Skip HITindex positional classification. Default: off; cannot be combined with --calculate-afe-ale or --hitindex-dir."),
                    ("--hitindex-dir", "Directory containing reusable per-replicate HITindex outputs. Default: unset; recompute HITindex."),
                ],
            ),
            "",
            help_box(
                "Examples",
                [
                    ("TExTra qual --prep result -o result -s sample1,sample2", "Run qual after prep"),
                    ("TExTra qual --reuse -o result", "Infer context from TExTra.config.json"),
                ],
            ),
        ]
        return "\n".join(lines) + "\n"


def parse_arguments(args_list):
    """Parse qual CLI arguments."""
    parser = QualArgumentParser(
        prog="TExTra qual",
        description="Identify TE-overlap exons and optional HITindex positional evidence from prep outputs.",
        add_help=True,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-s",
        "--samples",
        default=None,
        help="Sample/condition names, comma-separated; inferred with --reuse when omitted if possible.",
    )
    add_threading_arguments(parser)
    parser.add_argument("--prep", default=None, help="Path to TExTra prep output; required unless inferred by --reuse.")
    parser.add_argument("-o", "--out_dir", required=True, help="Qual output root directory.")
    add_project_argument(parser, default="project", help_text="Project prefix for classify combined outputs")
    add_read_layout_arguments(parser)
    parser.add_argument(
        "--ss3-buffer",
        dest="ss3buffer",
        type=int,
        default=20,
        help="3' splice-site buffer size for HITindex splice-site roles.",
    )
    parser.add_argument(
        "--ss5-buffer",
        dest="ss5buffer",
        type=int,
        default=50,
        help="5' splice-site buffer size for HITindex splice-site roles.",
    )
    parser.add_argument(
        "--calculate-afe-ale",
        dest="calculate_afe_ale",
        action="store_true",
        default=False,
        help="Calculate AFE/ALE exon usage outputs. Requires HITindex; cannot be combined with --skip-hitindex.",
    )
    parser.set_defaults(te_overlap_junction_evidence=True)
    parser.add_argument(
        "--ignore-junction",
        dest="te_overlap_junction_evidence",
        action="store_false",
        default=argparse.SUPPRESS,
        help="Ignore TE-overlap junction support/degradation checks; TE overlap annotation is still generated.",
    )
    parser.add_argument(
        "-j",
        "--junction",
        type=float,
        default=2.0,
        help=(
            "Minimum TE-side junction reads required to retain a TE-overlap call during junction evidence "
            "degradation. Default: 2.0. If omitted and prep used a non-default --junction, inherit the prep value "
            "from TExTra.config.json and print an INFO message. "
            "Cannot be combined with --ignore-junction or --skip-hitindex."
        ),
    )
    parser.add_argument("--genmodel-iters", type=int, default=100000, help="ADVI iterations in HITindex generative model.")
    parser.add_argument("--bootstrap-n", type=int, default=1000, help="Bootstrap iterations for HITindex significance.")
    parser.add_argument("--seed", type=int, default=None, help="Optional random seed for HITindex model fitting and bootstrap.")
    parser.add_argument("--te-overlap-min-bp", type=int, default=10, help="Minimum overlap bp for TE-overlap metrics.")
    parser.add_argument("--te-overlap-min-frac", type=float, default=0.1, help="Minimum overlap fraction for TE-overlap metrics.")
    parser.add_argument(
        "--splice-site-flank-bp",
        type=int,
        default=10,
        help="Flank window (bp) around splice sites for boundary/anchor/var-site hit checks.",
    )
    parser.add_argument(
        "--skip-hitindex",
        action="store_true",
        help=(
            "Skip HITindex positional classification and produce annotation/TE-overlap tables only. "
            "Cannot be combined with --calculate-afe-ale or --hitindex-dir; junction evidence is disabled."
        ),
    )
    parser.add_argument(
        "--hitindex-dir",
        default=None,
        help=(
            "Directory containing reusable per-replicate HITindex outputs. "
            "If omitted, HITindex is recomputed. Cannot be combined with --skip-hitindex. "
            "If junction support is missing/invalid, qual warns and continues as --ignore-junction; "
            "complete AFE/ALE files auto-enable AFE/ALE outputs."
        ),
    )
    parser.add_argument(
        "--reuse",
        action="store_true",
        help=(
            "Infer omitted --prep, --samples, --strand, and --readtype from the previous prep run context; "
            "does not reuse qual result files."
        ),
    )
    add_run_mode_arguments(parser, debug_help="Enable debug mode and keep intermediate files.")

    args = parser.parse_args(args_list)
    args._explicit_prep = _option_was_provided(args_list, "--prep")
    args._explicit_samples = _option_was_provided(args_list, "-s", "--samples")
    args._explicit_strand = _option_was_provided(args_list, "--strand")
    args._explicit_readtype = _option_was_provided(args_list, "--readtype")
    args._explicit_junction = _option_was_provided(args_list, "-j", "--junction")
    if getattr(args, "_explicit_junction", False) and float(args.junction) <= 0:
        raise RuntimeError("--junction must be a positive number when provided.")
    if getattr(args, "_explicit_junction", False) and not bool(args.te_overlap_junction_evidence):
        raise RuntimeError("--junction cannot be combined with --ignore-junction.")
    if getattr(args, "_explicit_junction", False) and bool(args.skip_hitindex):
        raise RuntimeError("--junction cannot be combined with --skip-hitindex because junction evidence is disabled.")
    return args


def _is_prep_output_dir(path):
    if not path:
        return False
    return (
        os.path.isfile(os.path.join(path, CONVERT_DIR, "TE_anno.bed"))
        and os.path.isfile(os.path.join(path, ASSEMBLY_DIR, CONSENSUS_GTF))
    )


def _discover_prep_dir(out_dir):
    candidates = []
    for path in [out_dir, os.getcwd()]:
        if not path:
            continue
        abs_path = os.path.abspath(path)
        config_path = global_config_path(abs_path)
        if not os.path.isfile(config_path):
            continue
        config, _path = read_global_config(abs_path)
        prep_config = extract_module_config(config, "prep")
        prep_dir = os.path.abspath(prep_config.get("out_dir", abs_path)) if prep_config else abs_path
        if prep_dir in candidates:
            continue
        if _is_prep_output_dir(prep_dir):
            candidates.append(prep_dir)
    for path in [out_dir, os.getcwd()]:
        if not path:
            continue
        abs_path = os.path.abspath(path)
        if abs_path in candidates:
            continue
        if _is_prep_output_dir(abs_path):
            candidates.append(abs_path)
    if len(candidates) > 1:
        raise RuntimeError(
            "--reuse found multiple possible prep output directories. "
            "Please pass --prep explicitly: " + ", ".join(candidates)
        )
    return candidates[0] if candidates else None


def _load_prep_config(prep_dir):
    config, config_path = read_global_config(prep_dir)
    prep_config = extract_module_config(config, "prep")
    if prep_config:
        return prep_config, config_path

    debug_config_path = os.path.join(prep_dir, "logs", "prep_config.json")
    if not os.path.isfile(debug_config_path):
        return None, config_path
    with open(debug_config_path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise RuntimeError(f"Invalid prep run config: {debug_config_path} must contain a JSON object.")
    return payload, debug_config_path


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


def _set_junction_threshold_from_prep(args, config):
    default_threshold = 2.0
    value = getattr(args, "junction", default_threshold)
    if not getattr(args, "_explicit_junction", False) and isinstance(config, dict):
        prep_value = config.get("junction")
        try:
            prep_threshold = float(prep_value)
        except (TypeError, ValueError):
            prep_threshold = default_threshold
        if prep_threshold > 0 and prep_threshold != default_threshold:
            value = prep_threshold
            log_message(
                "[INFO]",
                (
                    "Qual junction threshold inherited from prep: "
                    f"--junction {prep_threshold:g} (qual default: {default_threshold:g})."
                ),
                color="info",
            )
    try:
        threshold = float(value)
    except (TypeError, ValueError):
        threshold = default_threshold
    if threshold <= 0:
        threshold = default_threshold
    args.junction = threshold
    args.te_boundary_min_side_reads = threshold


def _apply_reuse_context(args):
    if not args.reuse:
        missing = []
        if not args.prep:
            missing.append("--prep")
        if not args.samples:
            missing.append("--samples")
        if missing:
            raise RuntimeError(
                "Missing required qual argument(s): "
                + ", ".join(missing)
                + ". Provide them explicitly, or rerun with --reuse to infer them from prep output."
            )
        config, _config_path = _load_prep_config(args.prep)
        _set_junction_threshold_from_prep(args, config)
        return args

    reused = []
    if not args.prep:
        discovered = _discover_prep_dir(args.out_dir)
        if discovered:
            args.prep = discovered
            reused.append(f"prep={discovered}")
    if not args.prep:
        raise RuntimeError("--reuse could not locate a prep output directory. Pass --prep explicitly.")

    config, config_path = _load_prep_config(args.prep)
    if config is None:
        log_message(
            "[WARNING]",
            f"--reuse requested, but prep config was not found: {config_path}. "
            "Only file-based inference will be used.",
            color="warning",
        )

    if not args.samples:
        samples = []
        if config:
            samples = [str(sample) for sample in config.get("samples", []) if str(sample).strip()]
        if not samples:
            samples = _infer_samples_from_alignment(args.prep)
        if samples:
            args.samples = ",".join(samples)
            reused.append(f"samples={len(samples)}")
    if not args.samples:
        raise RuntimeError("--reuse could not infer --samples from prep output. Pass --samples explicitly.")

    if config and not args._explicit_strand and config.get("strand"):
        args.strand = str(config["strand"])
        reused.append(f"strand={args.strand}")
    if config and not args._explicit_readtype and config.get("readtype"):
        args.readtype = str(config["readtype"])
        reused.append(f"readtype={args.readtype}")
    _set_junction_threshold_from_prep(args, config)

    if reused:
        log_message("[INFO]", "Reused prep context: " + ", ".join(reused), color="info")
    return args


def main(args_list=None):
    """Run the qual CLI entrypoint."""
    if args_list is None:
        args_list = sys.argv[2:] if len(sys.argv) > 1 and sys.argv[1] == "qual" else sys.argv[1:]
    args = parse_arguments(args_list)
    if args.debug:
        args.detail = True
    set_log_file(os.path.join(args.out_dir, "logs", "qual.log"))
    args = _apply_reuse_context(args)
    identify_func(args)


if __name__ == "__main__":
    main()

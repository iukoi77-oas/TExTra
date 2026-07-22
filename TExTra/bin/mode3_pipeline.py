#!/usr/bin/env python
"""CLI entrypoint for TExTra diff and optional ncPred analysis."""

import sys
import os
import argparse
import json
import glob

# Allow direct execution from the source tree without requiring installation.
script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
project_root = os.path.dirname(parent_dir)
sys.path.insert(0, os.path.abspath(project_root))
sys.path.insert(0, os.path.abspath(parent_dir))

from util.common.write_logs import log_message, set_log_file
from util.common.define_layout import DOWNSTREAM_DIR, QUANT_DIR, resolve_consensus_gtf
from util.common.external_tools import resolve_external_file
from util.common.run_config import extract_module_config, global_config_path, read_global_config, update_global_config
from util.common.help_format import help_box, help_title
from util.common.define_cli import (
    add_project_argument,
    add_run_mode_arguments,
)


def _option_was_provided(args_list, *names):
    tokens = list(args_list or [])
    for token in tokens:
        for name in names:
            if token == name or token.startswith(name + "="):
                return True
    return False


class DiffArgumentParser(argparse.ArgumentParser):
    """argparse parser with custom diff help layout."""

    def format_help(self):
        lines = [
            help_title("TExTra diff"),
            "Test differential TE-overlap exon usage and optionally run coding-potential prediction.",
            "",
            help_box(
                "Usage",
                [
                    ("TExTra diff --prep PREP --quant QUANT -o OUT_DIR --groups A,B [options]", "Run differential usage"),
                    ("TExTra diff --reuse -o OUT_DIR [options]", "Infer prep/quant context when possible"),
                    ("TExTra diff ... --ncpred --genome GENOME", "Also run ncPred coding-potential prediction"),
                ],
            ),
            "",
            help_box(
                "Required arguments",
                [
                    ("-o, --out_dir", "Output root directory; diff writes the 05_downstream subdirectory."),
                    ("--prep", "Path to TExTra prep output; required unless inferred by --reuse."),
                    ("--quant", "Path to TExTra quant output; required unless inferred by --reuse."),
                    ("--groups", "Two condition names to compare, comma-separated; required unless --reuse can infer exactly two groups."),
                ],
            ),
            "",
            help_box(
                "Optional arguments",
                [
                    ("-h, --help", "Show this message and exit."),
                    ("--project", "Project name. Default: sample1_vs_sample2."),
                    ("--test-method", "Differential usage test method: classical or empirical. Default: classical."),
                    ("--padj", "Adjusted p-value threshold for significant differential usage. Default: 0.05."),
                    ("--pvalue", "Optional raw p-value threshold; disabled by default."),
                    ("--delta-exon-usage", "Minimum absolute exon usage difference for significant differential usage. Default: 0.1."),
                    ("--paired", "Use paired replicate testing for classical mode. Default: off."),
                    ("--ncpred", "Run ncPred and write ncPred/plek_result.csv plus ncPred/selected_transcripts.gtf. Default: off."),
                    ("-g, --genome", "Genome FASTA file. Default: unset; required only when --ncpred is enabled."),
                    ("--plek-path", "Path to PLEK2.py or a PLEK installation directory. Default: auto-discover bundled/external PLEK."),
                    ("--ncpred-model", "PLEK2 coding-potential model: ve or pl. Default: ve."),
                    ("--min-length", "Minimum transcript length retained for ncPred candidates; used only with --ncpred. Default: 200."),
                    ("--debug", "Enable debug mode and keep intermediate files. Default: off."),
                    ("--detail", "Enable detail mode for additional result-checking tables and summaries. Default: off."),
                    ("--reuse", "Infer omitted prep/quant/project/groups context; does not reuse diff result files. Default: off."),
                ],
            ),
            "",
            help_box(
                "Advanced options",
                [
                    ("--empirical-background-n", "Local background size for empirical mode; nearest replicate deltas used per event. Default: 1000."),
                ],
            ),
            "",
            help_box(
                "Examples",
                [
                    ("TExTra diff --prep result --quant result -o result --groups condA,condB", "Run classical diff"),
                    ("TExTra diff --reuse -o result", "Infer context from TExTra.config.json"),
                    ("TExTra diff --prep result --quant result -o result --groups condA,condB --ncpred -g genome.fa", "Run diff + ncPred"),
                ],
            ),
        ]
        return "\n".join(lines) + "\n"


def parse_arguments(args_list):
    parser = DiffArgumentParser(
        prog="TExTra diff",
        description="Differential TE-overlap exon usage analysis with optional coding-potential prediction.",
        add_help=True,  
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-g", "--genome", type=str, required=False,
                        help="Genome FASTA file; required only when --ncpred is enabled.")
    parser.add_argument("-o", "--out_dir", type=str, required=True,
                        help=f"Output root directory; diff writes the {DOWNSTREAM_DIR} subdirectory.")
    parser.add_argument("--prep", default=None, help="Path to TExTra prep output; required unless inferred by --reuse.")
    parser.add_argument("--quant", default=None, help="Path to TExTra quant output containing 04_quantification/<project>.TE_overlap.exon_usage.tsv; required unless inferred by --reuse.")
    
    add_project_argument(parser, default="sample1_vs_sample2", help_text="Project name.")
    parser.add_argument("--ncpred", action="store_true",
                        help="Run ncPred coding-potential analysis and write ncPred/plek_result.csv plus ncPred/selected_transcripts.gtf.")
    parser.add_argument(
        "--plek-path",
        default=None,
        help=(
            "Path to PLEK2.py or a PLEK installation directory; used only with --ncpred. "
            "If omitted, TExTra searches TEXTRA_EXTERNAL_DIR/PLEK*/PLEK2.py and util/external/PLEK*/PLEK2.py."
        ),
    )
    parser.add_argument("--test-method", choices=["classical", "empirical"], default="classical",
                        help="Differential usage test method. 'classical' compares replicate-level usage directly; 'empirical' uses transcript-abundance-aware local background.")
    parser.add_argument("--min-length", type=int, default=200,
                        help="Minimum transcript length retained for ncPred candidate transcripts; used only with --ncpred.")
    parser.add_argument("--padj", type=float, default=0.05,
                        help="Adjusted p-value threshold for significant differential usage.")
    parser.add_argument("--pvalue", type=float, default=None,
                        help="Optional raw p-value threshold. By default only --padj and --delta-exon-usage define significance.")
    parser.add_argument("--delta-exon-usage", dest="delta_exon_usage", type=float, default=0.1,
                        help="Minimum absolute exon usage difference for significant differential usage.")
    parser.add_argument("--empirical-background-n", type=int, default=1000,
                        help="Local background size for empirical mode; number of nearest replicate deltas used per event.")
    parser.add_argument("--paired", action="store_true",
                        help="Use paired replicate testing for classical mode.")
    
    parser.add_argument("--groups", default=None,
                        help="Two condition names to compare, comma-separated; names must match condition prefixes inferred from quant sample columns.")
    parser.add_argument("--ncpred-model", choices=["ve", "pl"], default="ve",
                        help="PLEK2 coding-potential model: ve for vertebrate, pl for plant.")
    add_run_mode_arguments(parser)
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="Infer omitted --prep, --quant, --project, --groups, and --genome from run context; does not reuse diff result files.",
    )

    args = parser.parse_args(args_list)
    args._explicit_prep = _option_was_provided(args_list, "--prep")
    args._explicit_quant = _option_was_provided(args_list, "--quant")
    args._explicit_groups = _option_was_provided(args_list, "--groups")
    args._explicit_project = _option_was_provided(args_list, "--project")
    args._explicit_genome = _option_was_provided(args_list, "-g", "--genome")
    return args


def _is_prep_output_dir(path):
    if not path:
        return False
    return os.path.isfile(resolve_consensus_gtf(os.path.abspath(path)))


def _is_quant_output_dir(path):
    if not path:
        return False
    quant_dir = os.path.join(os.path.abspath(path), QUANT_DIR)
    return os.path.isdir(quant_dir) and bool(glob.glob(os.path.join(quant_dir, "*.TE_overlap.exon_usage.tsv")))


def _single_or_error(candidates, label):
    if len(candidates) > 1:
        raise RuntimeError(
            f"--reuse found multiple possible {label} output directories. "
            f"Please pass --{label} explicitly: " + ", ".join(candidates)
        )
    return candidates[0] if candidates else None


def _discover_reuse_context(out_dir):
    prep_candidates = []
    quant_candidates = []
    prep_config = None
    quant_config = None
    used_config_path = None

    for root in [out_dir, os.getcwd()]:
        if not root:
            continue
        root = os.path.abspath(root)
        config_path = global_config_path(root)
        if os.path.isfile(config_path):
            config, path = read_global_config(root)
            used_config_path = path
            prep_cfg = extract_module_config(config, "prep")
            quant_cfg = extract_module_config(config, "quant")
            if prep_cfg and prep_config is None:
                prep_config = prep_cfg
            if quant_cfg and quant_config is None:
                quant_config = quant_cfg
            if prep_cfg:
                prep_dir = os.path.abspath(prep_cfg.get("out_dir", root))
                if _is_prep_output_dir(prep_dir) and prep_dir not in prep_candidates:
                    prep_candidates.append(prep_dir)
            if quant_cfg:
                quant_dir = os.path.abspath(quant_cfg.get("out_dir", root))
                if _is_quant_output_dir(quant_dir) and quant_dir not in quant_candidates:
                    quant_candidates.append(quant_dir)
                prep_dir = os.path.abspath(quant_cfg.get("prep", root))
                if _is_prep_output_dir(prep_dir) and prep_dir not in prep_candidates:
                    prep_candidates.append(prep_dir)
        if _is_prep_output_dir(root) and root not in prep_candidates:
            prep_candidates.append(root)
        if _is_quant_output_dir(root) and root not in quant_candidates:
            quant_candidates.append(root)

    return prep_candidates, quant_candidates, prep_config, quant_config, used_config_path


def _infer_groups_from_quant_config(quant_config):
    samples = [str(sample).strip() for sample in (quant_config or {}).get("samples", []) if str(sample).strip()]
    if len(samples) == 2:
        return ",".join(samples)
    if len(samples) > 2:
        raise RuntimeError(
            "--reuse found more than two quant sample groups. "
            "Pass --groups explicitly to choose the comparison: " + ", ".join(samples)
        )
    return None


def _apply_reuse_context(args):
    if not args.reuse:
        missing = []
        if not args.prep:
            missing.append("--prep")
        if not args.quant:
            missing.append("--quant")
        if not args.groups:
            missing.append("--groups")
        if missing:
            raise RuntimeError(
                "Missing required diff argument(s): "
                + ", ".join(missing)
                + ". Provide them explicitly, or rerun with --reuse to infer them from run config."
            )
        return args

    prep_candidates, quant_candidates, prep_config, quant_config, config_path = _discover_reuse_context(args.out_dir)
    reused = []

    if not args.prep:
        args.prep = _single_or_error(prep_candidates, "prep")
        if args.prep:
            reused.append(f"prep={args.prep}")
    if not args.quant:
        args.quant = _single_or_error(quant_candidates, "quant")
        if args.quant:
            reused.append(f"quant={args.quant}")
    if not args.prep:
        raise RuntimeError("--reuse could not locate a prep output directory. Pass --prep explicitly.")
    if not args.quant:
        raise RuntimeError("--reuse could not locate a quant output directory. Pass --quant explicitly.")

    if quant_config and not args._explicit_project and quant_config.get("project"):
        args.project = str(quant_config["project"])
        reused.append(f"project={args.project}")
    if not args.groups:
        inferred_groups = _infer_groups_from_quant_config(quant_config)
        if inferred_groups:
            args.groups = inferred_groups
            reused.append(f"groups={args.groups}")
    if not args.groups:
        raise RuntimeError("--reuse could not infer --groups from quant context. Pass --groups explicitly.")

    if args.ncpred and not args.genome:
        genome = None
        if quant_config and quant_config.get("genome"):
            genome = quant_config.get("genome")
        elif prep_config and prep_config.get("genome"):
            genome = prep_config.get("genome")
        if genome:
            args.genome = str(genome)
            reused.append(f"genome={args.genome}")

    if not (prep_config or quant_config):
        log_message(
            "[WARNING]",
            f"--reuse did not find module run config in {config_path or global_config_path(args.out_dir)}; file-based inference was used.",
            color="warning",
        )
    if reused:
        log_message("[INFO]", "Reused diff context: " + ", ".join(reused), color="info")
    return args


def _validate_inputs(args):
    """Validate mode3 upstream outputs and ncPred-specific inputs."""
    if not os.path.isdir(args.prep):
        raise FileNotFoundError(f"prep directory not found: {args.prep}")
    if not os.path.isdir(args.quant):
        raise FileNotFoundError(f"quant directory not found: {args.quant}")
    quant_table = os.path.join(args.quant, QUANT_DIR)
    if not os.path.isdir(quant_table):
        raise FileNotFoundError(f"{QUANT_DIR} directory not found: {quant_table}")
    if args.ncpred:
        if not args.genome:
            raise ValueError("--ncpred requires --genome")
        if not os.path.isfile(args.genome):
            raise FileNotFoundError(f"genome FASTA not found: {args.genome}")
        args.plek2_path = _resolve_explicit_or_auto_plek2_path(args.plek_path)


def _resolve_explicit_or_auto_plek2_path(plek_path):
    if plek_path:
        path = os.path.abspath(plek_path)
        if os.path.isfile(path):
            return path
        if os.path.isdir(path):
            direct = os.path.join(path, "PLEK2.py")
            if os.path.isfile(direct):
                return os.path.abspath(direct)
            matches = sorted(glob.glob(os.path.join(path, "PLEK*", "PLEK2.py")))
            for match in matches:
                if os.path.isfile(match):
                    return os.path.abspath(match)
        raise FileNotFoundError(
            "--plek-path must point to PLEK2.py or a directory containing PLEK2.py "
            f"(directly or under PLEK*/). Got: {plek_path}"
        )
    plek2_path, candidates = resolve_external_file("PLEK*", "PLEK2.py", anchor_file=__file__)
    if plek2_path:
        return plek2_path
    searched = ", ".join(candidates) if candidates else os.path.abspath(os.path.join("util", "external"))
    raise FileNotFoundError(
        "PLEK2.py not found for --ncpred. "
        "Pass --plek-path, or place PLEK2 under TEXTRA_EXTERNAL_DIR/PLEK*/PLEK2.py or util/external/PLEK*/PLEK2.py. "
        f"Searched external roots: {searched}"
    )


def _validate_parameters(args):
    if int(args.min_length) < 1:
        raise RuntimeError("--min-length must be a positive integer.")
    if not (0.0 <= float(args.padj) <= 1.0):
        raise RuntimeError("--padj must be between 0 and 1.")
    if args.pvalue is not None and not (0.0 <= float(args.pvalue) <= 1.0):
        raise RuntimeError("--pvalue must be between 0 and 1 when provided.")
    if float(args.delta_exon_usage) < 0:
        raise RuntimeError("--delta-exon-usage must be non-negative.")
    if int(args.empirical_background_n) < 1:
        raise RuntimeError("--empirical-background-n must be a positive integer.")


def _write_diff_run_config(args):
    payload = {
        "out_dir": os.path.abspath(args.out_dir),
        "output_root": os.path.abspath(args.output_root),
        "prep": os.path.abspath(args.prep),
        "quant": os.path.abspath(args.quant),
        "genome": os.path.abspath(args.genome) if args.genome else None,
        "project": args.project,
        "groups": [group.strip() for group in str(args.groups).split(",") if group.strip()],
        "test_method": args.test_method,
        "padj": float(args.padj),
        "pvalue": None if args.pvalue is None else float(args.pvalue),
        "delta_exon_usage": float(args.delta_exon_usage),
        "empirical_background_n": int(args.empirical_background_n),
        "paired": bool(args.paired),
        "ncpred": bool(args.ncpred),
        "ncpred_model": args.ncpred_model,
        "plek_path": os.path.abspath(args.plek_path) if args.plek_path else None,
        "plek2_path": os.path.abspath(getattr(args, "plek2_path", "")) if getattr(args, "plek2_path", None) else None,
        "min_length": int(args.min_length),
        "reuse": bool(getattr(args, "reuse", False)),
        "detail": bool(args.detail),
        "debug": bool(args.debug),
    }
    update_global_config(args.output_root, "diff", payload)
    if args.debug:
        logs_dir = os.path.join(args.output_root, "logs")
        os.makedirs(logs_dir, exist_ok=True)
        debug_payload = {"module": "diff", "schema_version": 1, "run_mode": "debug", **payload}
        with open(os.path.join(logs_dir, "diff_config.json"), "w", encoding="utf-8") as handle:
            json.dump(debug_payload, handle, indent=2, sort_keys=True)
    return payload


def main(args_list=None):
    """Run differential usage analysis and optional ncPred from parsed CLI arguments."""
    if args_list is None:
        args_list = sys.argv[2:] if len(sys.argv) > 1 and sys.argv[1] == "diff" else sys.argv[1:]
    args = parse_arguments(args_list)
    if args.debug:
        args.detail = True
    set_log_file(os.path.join(args.out_dir, "logs", "diff.log"))
    args = _apply_reuse_context(args)
    _validate_parameters(args)
    _validate_inputs(args)
    args.output_root = os.path.abspath(args.out_dir)
    args.out_dir = os.path.join(args.output_root, DOWNSTREAM_DIR)
    if args.debug:
        args.diff_extend_log = os.path.join(args.output_root, "logs", "diff_extend.log")
        os.makedirs(os.path.dirname(args.diff_extend_log), exist_ok=True)
        with open(args.diff_extend_log, "w", encoding="utf-8") as handle:
            handle.write("diff debug command/output log\n")
    else:
        args.diff_extend_log = None
    args.diff_total_steps = 2 if args.ncpred else 1
    os.makedirs(args.out_dir, exist_ok=True)
    _write_diff_run_config(args)

    from TExTra.src.mode3.DE import DE_func

    DE_func(args)
    
    if args.ncpred:
        from TExTra.src.mode3.ncPred import ncPred_func

        ncPred_func(args)

if __name__ == "__main__":
    main()

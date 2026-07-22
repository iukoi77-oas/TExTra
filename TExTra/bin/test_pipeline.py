"""CLI entrypoint for the repository test workflow."""

import argparse
import os
import sys

# Allow direct execution from the source tree and locate repository test assets.
script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
repo_dir = os.path.dirname(parent_dir)
sys.path.insert(0, os.path.abspath(repo_dir))
sys.path.insert(0, os.path.abspath(parent_dir))

from util.common.define_cli import (
    add_project_argument,
    add_read_layout_arguments,
    add_run_mode_arguments,
    add_threading_arguments,
    validate_read_layout,
    validate_threading_args,
)
from util.common.help_format import help_box, help_title
from util.common.write_logs import log_message

TEST_GENMODEL_ITERS = 1000
TEST_BOOTSTRAP_N = 100


def _print_module_separator(index, module_name):
    line = f"==================== TExTra {index}: {module_name} ===================="
    print(f"\033[1;95m{line}\033[0m")


class TestArgumentParser(argparse.ArgumentParser):
    """argparse parser with custom test help layout."""

    def format_help(self):
        lines = [
            help_title("TExTra test"),
            "Run prep + qual + quant + diff on a small test dataset.",
            "",
            help_box(
                "Usage",
                [
                    ("TExTra test [options]", "Run the bundled test workflow"),
                    ("TExTra test --test-data-dir TEST_DATA_DIR [options]", "Run with an external test-data directory"),
                    ("TExTra test ... --skip-diff", "Run prep + qual + quant only"),
                ],
            ),
            "",
            help_box(
                "Optional arguments",
                [
                    ("-h, --help", "Show this message and exit."),
                    ("-o, --out_dir", "Output root directory. Default: test/result under the detected test root."),
                    ("--test-data-dir", "Directory containing test reads and default references. Default: auto-detect test/example_data."),
                    ("-g, --genome", "Genome FASTA path; if omitted, use <test-data-dir>/reference.fa."),
                    ("-G, --gene", "Gene annotation GTF path; if omitted, use <test-data-dir>/gencode.vM21.gtf."),
                    ("-r, --te", "TE annotation path; if omitted, use <test-data-dir>/TE_rmsk.gtf."),
                    ("--input-tsv", "Input TSV template resolved against --test-data-dir. Default: auto-detect test/input.tsv."),
                    ("-t, --threads", "Threads per module. Default: 4."),
                    ("--njobs", "Maximum number of parallel jobs. Default: omitted, use --threads."),
                    ("--seed", "Optional random seed passed to qual HITindex steps. Default: unset."),
                    ("--project", "Project prefix. Default: project."),
                    ("--strand", "Library strand type: none, rf, fr, r, f. Default: rf."),
                    ("--readtype", "Read type: paired or single. Default: paired."),
                    ("--debug", "Enable debug mode and keep intermediate files. Default: off."),
                    ("--detail", "Enable detail mode for additional result-checking tables and summaries. Default: off."),
                ],
            ),
            "",
            help_box(
                "Test workflow",
                [
                    ("--assembly-mode", "StringTie assembly mode passed to prep: de-novo or reference-guided. Default: reference-guided."),
                    ("--te-overlap-min-bp", "Minimum overlap bp for TE-overlap metrics. Default: 10."),
                    ("--te-overlap-min-frac", "Minimum overlap fraction for TE-overlap metrics. Default: 0.1."),
                    ("--splice-site-flank-bp", "Flank window around splice sites for boundary/anchor/var-site checks. Default: 10."),
                    ("--ignore-junction", "Ignore TE-overlap junction support/degradation checks. Default: off; TE overlap annotation is still generated."),
                    ("--quantifier", "Quantification backend passed to quant: rsem or salmon. Default: rsem."),
                    ("--skip-diff", "Skip the diff step and stop after quant. Default: off."),
                    ("--diff-samples", "Two condition names for diff, comma-separated; if omitted, infer the first two input groups."),
                    ("--diff-padj", "Adjusted p-value threshold for diff. Default: 0.05."),
                    ("HITindex speed", f"Test workflow uses --genmodel-iters {TEST_GENMODEL_ITERS} and --bootstrap-n {TEST_BOOTSTRAP_N}."),
                ],
            ),
            "",
            help_box(
                "Examples",
                [
                    ("TExTra test", "Run the bundled full test with default paths"),
                    ("TExTra test --test-data-dir test/example_data", "Run with an explicit test-data directory"),
                    ("TExTra test --detail", "Run with detail result tables"),
                ],
            ),
        ]
        return "\n".join(lines) + "\n"


def parse_arguments(args_list):
    parser = TestArgumentParser(
        prog="TExTra test",
        description="Run prep + qual + quant + diff on a small test dataset.",
        add_help=True,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--test-data-dir",
        default=None,
        help="Directory containing test BAM/FASTQ files and default references.",
    )
    parser.add_argument("-o", "--out_dir", default=None, help="Output root directory.")
    parser.add_argument("-g", "--genome", default=None, help="Genome FASTA path; if omitted, use <test-data-dir>/reference.fa.")
    parser.add_argument("-G", "--gene", default=None, help="Gene annotation GTF path; if omitted, use <test-data-dir>/gencode.vM21.gtf.")
    parser.add_argument("-r", "--te", default=None, help="TE annotation path; if omitted, use <test-data-dir>/TE_rmsk.gtf.")
    parser.add_argument("--input-tsv", default=None, help="Input TSV template path resolved against --test-data-dir.")
    add_threading_arguments(parser, threads_help="Threads per module.")
    parser.add_argument("--seed", type=int, default=None, help="Optional random seed passed to qual HITindex steps.")
    add_project_argument(parser, default="project", help_text="Project prefix.")
    add_read_layout_arguments(parser, strand_default="rf", strand_help="Library strand type.")
    parser.add_argument(
        "--assembly-mode",
        choices=["de-novo", "reference-guided"],
        default="reference-guided",
        help="Per-replicate StringTie assembly mode passed to prep.",
    )
    add_run_mode_arguments(parser, debug_help="Enable debug mode and keep intermediate files.")
    parser.add_argument("--te-overlap-min-bp", type=int, default=10, help="Minimum overlap bp for TE-overlap metrics.")
    parser.add_argument("--te-overlap-min-frac", type=float, default=0.1, help="Minimum overlap fraction for TE-overlap metrics.")
    parser.add_argument("--splice-site-flank-bp", type=int, default=10, help="Flank window (bp) around splice sites for boundary/anchor/var-site hit checks.")
    parser.set_defaults(te_overlap_junction_evidence=True)
    parser.add_argument(
        "--ignore-junction",
        dest="te_overlap_junction_evidence",
        action="store_false",
        default=argparse.SUPPRESS,
        help="Ignore TE-overlap junction support/degradation checks in qual; TE overlap annotation is still generated.",
    )
    parser.add_argument(
        "--skip-diff",
        action="store_true",
        help="Skip diff step in test workflow.",
    )
    parser.add_argument(
        "--quantifier",
        choices=["rsem", "salmon"],
        default="rsem",
        help="Quantification backend passed to quant.",
    )
    parser.add_argument(
        "--diff-samples",
        default=None,
        help="Two condition names for diff, comma-separated; if omitted, infer the first two conditions from input TSV first column.",
    )
    parser.add_argument("--diff-padj", type=float, default=0.05, help="Adjusted p-value threshold for diff step.")
    return parser.parse_args(args_list)


def _validate_parameters(args):
    validate_threading_args(args)
    validate_read_layout(args, module_name="test")
    if int(args.te_overlap_min_bp) < 1:
        raise RuntimeError("--te-overlap-min-bp must be a positive integer.")
    if not (0.0 <= float(args.te_overlap_min_frac) <= 1.0):
        raise RuntimeError("--te-overlap-min-frac must be between 0 and 1.")
    if int(args.splice_site_flank_bp) < 0:
        raise RuntimeError("--splice-site-flank-bp must be a non-negative integer.")
    if not (0.0 <= float(args.diff_padj) <= 1.0):
        raise RuntimeError("--diff-padj must be between 0 and 1.")


def _resolve_cell_path(test_data_dir, token):
    value = str(token).strip()
    if not value or value in {"NA", "NaN", "nan", "-"}:
        return value
    if os.path.isabs(value):
        return value
    return os.path.abspath(os.path.join(test_data_dir, value))


def _build_resolved_input(input_tsv, test_data_dir, out_dir):
    if not os.path.isfile(input_tsv):
        raise FileNotFoundError(f"test input TSV not found: {input_tsv}")
    if not os.path.isdir(test_data_dir):
        raise FileNotFoundError(f"test data directory not found: {test_data_dir}")

    rows = []
    with open(input_tsv, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            fixed = [parts[0].strip()]
            for cell in parts[1:]:
                chunks = [_resolve_cell_path(test_data_dir, x) for x in str(cell).split(",")]
                fixed.append(",".join(chunks))
            rows.append(fixed)

    if not rows:
        raise RuntimeError(f"No valid rows found in {input_tsv}")

    tmp_dir = os.path.join(out_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    resolved_path = os.path.join(tmp_dir, "test_input.resolved.tsv")
    with open(resolved_path, "w", encoding="utf-8") as out:
        for row in rows:
            out.write("\t".join(row) + "\n")
    return resolved_path


def _resolve_reference_path(test_data_dir, user_value, default_name, label):
    if user_value:
        path = user_value if os.path.isabs(user_value) else os.path.abspath(user_value)
    else:
        path = os.path.abspath(os.path.join(test_data_dir, default_name))
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{label} not found: {path}")
    return path


def _existing_default_path(user_value, candidates, label, expect_dir=False):
    if user_value:
        path = user_value if os.path.isabs(user_value) else os.path.abspath(user_value)
    else:
        path = None
        for candidate in candidates:
            if expect_dir and os.path.isdir(candidate):
                path = candidate
                break
            if not expect_dir and os.path.isfile(candidate):
                path = candidate
                break
        if path is None:
            searched = ", ".join(os.path.abspath(x) for x in candidates)
            raise FileNotFoundError(f"{label} not found. Searched: {searched}")
    path = os.path.abspath(path)
    if expect_dir and not os.path.isdir(path):
        raise FileNotFoundError(f"{label} not found: {path}")
    if not expect_dir and not os.path.isfile(path):
        raise FileNotFoundError(f"{label} not found: {path}")
    return path


def _default_roots():
    roots = []
    for root in [os.environ.get("TEXTRA_HOME"), os.getcwd(), repo_dir]:
        if not root:
            continue
        root = os.path.abspath(root)
        if root not in roots:
            roots.append(root)
    return roots


def _resolve_test_defaults(args):
    roots = _default_roots()
    args.test_data_dir = _existing_default_path(
        args.test_data_dir,
        [os.path.join(root, "test", "example_data") for root in roots],
        "test data directory",
        expect_dir=True,
    )
    args.input_tsv = _existing_default_path(
        args.input_tsv,
        [os.path.join(root, "test", "input.tsv") for root in roots],
        "test input TSV",
        expect_dir=False,
    )
    if args.out_dir is None:
        for root in roots:
            input_candidate = os.path.join(root, "test", "input.tsv")
            data_candidate = os.path.join(root, "test", "example_data")
            if os.path.isfile(input_candidate) or os.path.isdir(data_candidate):
                args.out_dir = os.path.join(root, "test", "result")
                break
        if args.out_dir is None:
            args.out_dir = os.path.join(os.getcwd(), "test", "result")
    elif not os.path.isabs(args.out_dir):
        args.out_dir = os.path.abspath(args.out_dir)
    return args


def _infer_diff_groups(input_tsv):
    groups = []
    seen = set()
    with open(input_tsv, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            group = line.rstrip("\n").split("\t", 1)[0].strip()
            if not group or group in {"NA", "NaN", "nan"}:
                continue
            if group in seen:
                continue
            seen.add(group)
            groups.append(group)
    if len(groups) < 2:
        raise RuntimeError(
            f"Cannot infer diff groups from {input_tsv}: need at least two condition names in first column."
        )
    return f"{groups[0]},{groups[1]}"


def main(args_list=None):
    if args_list is None:
        args_list = sys.argv[2:] if len(sys.argv) > 1 and sys.argv[1] == "test" else sys.argv[1:]
    args = parse_arguments(args_list)
    args = _resolve_test_defaults(args)
    _validate_parameters(args)
    from TExTra.bin.upstream_pipeline import main as upstream_main
    from TExTra.bin.mode3_pipeline import main as diff_main

    os.makedirs(args.out_dir, exist_ok=True)
    genome = _resolve_reference_path(args.test_data_dir, args.genome, "reference.fa", "genome FASTA")
    gene = _resolve_reference_path(args.test_data_dir, args.gene, "gencode.vM21.gtf", "gene annotation GTF")
    te = _resolve_reference_path(args.test_data_dir, args.te, "TE_rmsk.gtf", "TE annotation")
    resolved_input = _build_resolved_input(args.input_tsv, args.test_data_dir, args.out_dir)

    pipeline_args = [
        "-i", resolved_input,
        "-o", args.out_dir,
        "-g", genome,
        "-G", gene,
        "-r", te,
        "-t", str(args.threads),
        "--project", args.project,
        "--strand", args.strand,
        "--readtype", args.readtype,
        "--quantifier", args.quantifier,
        "--calculate-afe-ale",
        "--te-overlap-min-bp", str(args.te_overlap_min_bp),
        "--te-overlap-min-frac", str(args.te_overlap_min_frac),
        "--splice-site-flank-bp", str(args.splice_site_flank_bp),
        "--genmodel-iters", str(TEST_GENMODEL_ITERS),
        "--bootstrap-n", str(TEST_BOOTSTRAP_N),
    ]
    log_message(
        "[INFO]",
        "TExTra test uses fast HITindex settings: "
        f"--genmodel-iters {TEST_GENMODEL_ITERS}, --bootstrap-n {TEST_BOOTSTRAP_N}",
        color="info",
    )
    if args.njobs is not None:
        pipeline_args.extend(["--njobs", str(args.njobs)])
    if args.seed is not None:
        pipeline_args.extend(["--seed", str(args.seed)])
    if args.debug:
        pipeline_args.append("--debug")
    if args.detail:
        pipeline_args.append("--detail")
    pipeline_args.extend(["--assembly-mode", args.assembly_mode])
    if not args.te_overlap_junction_evidence:
        pipeline_args.append("--ignore-junction")

    upstream_main(pipeline_args)
    if args.skip_diff:
        return

    diff_samples = args.diff_samples.strip() if args.diff_samples else _infer_diff_groups(args.input_tsv)
    diff_args = [
        "--prep", args.out_dir,
        "--quant", args.out_dir,
        "-o", args.out_dir,
        "--project", args.project,
        "--groups", diff_samples,
        "--padj", str(args.diff_padj),
    ]
    if args.debug:
        diff_args.append("--debug")
    if args.detail:
        diff_args.append("--detail")
    _print_module_separator(4, "diff")
    diff_main(diff_args)


if __name__ == "__main__":
    main()

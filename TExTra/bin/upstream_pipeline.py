"""CLI entrypoint for the combined prep + qual + quant workflow."""

import argparse
import os
import sys

# Allow direct execution from the source tree without requiring installation.
script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
project_root = os.path.dirname(parent_dir)
sys.path.insert(0, os.path.abspath(project_root))
sys.path.insert(0, os.path.abspath(parent_dir))

from util.prep.define_arguments import add_prep_specific_arguments, build_prep_passthrough_args
from util.common.define_cli import (
    add_project_argument,
    add_read_layout_arguments,
    add_run_mode_arguments,
    add_threading_arguments,
    validate_read_layout,
    validate_threading_args,
)
from util.common.help_format import help_box, help_title


def _print_module_separator(index, module_name):
    line = f"==================== TExTra {index}: {module_name} ===================="
    print(f"\033[1;95m{line}\033[0m")


class UpstreamArgumentParser(argparse.ArgumentParser):
    """argparse parser with custom upstream help layout."""

    def format_help(self):
        lines = [
            help_title("TExTra upstream"),
            "Run prep + qual + quant in one command.",
            "",
            help_box(
                "Usage",
                [
                    ("TExTra upstream -i INPUT -o OUT_DIR -g GENOME -G GENE -r TE [options]", "Run the upstream workflow"),
                    ("TExTra upstream ... --calculate-afe-ale", "Also calculate AFE/ALE usage during qual"),
                ],
            ),
            "",
            help_box(
                "Required arguments",
                [
                    ("-i, --input", "Input TSV used by prep and quant read resolution."),
                    ("-o, --out_dir", "Output root directory shared by prep, qual, and quant."),
                    ("-g, --genome", "Genome FASTA used by prep and quant reference/index construction."),
                    ("-G, --gene", "Reference gene annotation GTF used by prep."),
                    ("-r, --te", "TE annotation file used by prep; supports GTF/GFF/BED and RepeatMasker OUT/TXT."),
                ],
            ),
            "",
            help_box(
                "Workflow options",
                [
                    ("-h, --help", "Show this message and exit."),
                    ("-s, --samples", "Sample/condition names, comma-separated; if omitted, infer from input first column."),
                    ("-t, --threads", "Threads per module. Default: 4."),
                    ("--njobs", "Maximum number of parallel jobs. Default: omitted, use --threads."),
                    ("--project", "Project prefix. Default: project."),
                    ("--strand", "Library strand type: none, rf, fr, r, f. Default: none."),
                    ("--readtype", "Read type: paired or single. Default: paired."),
                    ("--debug", "Enable debug mode and keep intermediate files. Default: off."),
                    ("--detail", "Enable detail mode for additional result-checking tables and summaries. Default: off."),
                ],
            ),
            "",
            help_box(
                "Prep options",
                [
                    ("--index", "Existing STAR genome index directory. Default: unset; prep builds one for FASTQ input."),
                    ("--extend", "Optional extra TE annotation merged with --te before TE BED conversion. Default: unset."),
                    ("--merge-method", "Transcript merge backend after StringTie assembly: taco or stringtie. Default: taco."),
                    ("--taco-path", "Path to taco_run or a TACO installation directory; used only with --merge-method taco. Default: auto-discover bundled/external TACO."),
                    ("--assembly-mode", "StringTie assembly mode: de-novo or reference-guided. Default: de-novo."),
                    ("--optimized", "Use optimized merge filtering defaults; explicit --min-* options override corresponding defaults. Default: off."),
                    ("--include-overlap", "Also retain lower-confidence reference-overlap transcript candidates during gene assignment. Default: off."),
                    ("-j, --junction", "StringTie first-pass minimum junction coverage and qual TE-side junction degradation threshold. Default: 2.0."),
                    ("--min-expr", "Merge minimum expression filter; unset unless provided or --optimized is used."),
                    ("--min-length", "Merge minimum transcript length filter; unset unless provided or --optimized is used."),
                    ("--min-frac", "Merge minimum isoform-fraction filter; unset unless provided or --optimized is used."),
                ],
            ),
            "",
            help_box(
                "Qual options",
                [
                    ("--calculate-afe-ale", "Enable AFE/ALE usage outputs. Requires HITindex. Default: off."),
                    ("--skip-hitindex", "Skip HITindex positional classification. Default: off; cannot be combined with --calculate-afe-ale."),
                    ("--ignore-junction", "Ignore TE-overlap junction support/degradation checks. Default: off; TE overlap annotation is still generated."),
                    ("--seed", "Optional random seed passed to HITindex steps. Default: unset."),
                    ("--te-overlap-min-bp", "Minimum overlap bp for TE-overlap metrics. Default: 10."),
                    ("--te-overlap-min-frac", "Minimum overlap fraction for TE-overlap metrics. Default: 0.1."),
                    ("--splice-site-flank-bp", "Flank window around splice sites for boundary/anchor/var-site checks. Default: 10."),
                ],
            ),
            "",
            help_box(
                "Quant options",
                [
                    ("--quantifier", "Quantification backend for quant: rsem or salmon. Default: rsem."),
                    ("--quant-result-dir", "Directory containing reusable RSEM/Salmon backend outputs. Default: unset; compute backend quantification as needed."),
                    ("--compute-gene-abundance", "Compute and export gene-level abundance tables in quant. Default: off."),
                ],
            ),
            "",
            help_box(
                "Advanced options",
                [
                    ("--gene-assignment", "Transcript-to-gene assignment policy after gffcompare. Default: strict."),
                    ("--sjdbOverhang", "STAR sjdbOverhang used when prep builds a STAR index. Default: 100."),
                    ("--seedMultimapNmax", "STAR seedMultimapNmax passed to STAR alignment. Default: 50000."),
                    ("--winAnchorMultimapNmax", "STAR winAnchorMultimapNmax passed to STAR alignment. Default: 100."),
                    ("--outFilterMultimapNmax", "STAR outFilterMultimapNmax passed to STAR alignment. Default: 100."),
                    ("--ss3-buffer", "3' splice-site buffer size for HITindex splice-site roles. Default: 20."),
                    ("--ss5-buffer", "5' splice-site buffer size for HITindex splice-site roles. Default: 50."),
                    ("--genmodel-iters", "HITindex ADVI iterations. Default: 100000."),
                    ("--bootstrap-n", "HITindex bootstrap iterations. Default: 1000."),
                ],
            ),
            "",
            help_box(
                "Examples",
                [
                    ("TExTra upstream -i input.tsv -o result -g genome.fa -G gene.gtf -r TE_annotation", "Run prep + qual + quant"),
                    ("TExTra upstream ... --calculate-afe-ale --detail", "Run with AFE/ALE and detail result tables"),
                ],
            ),
        ]
        return "\n".join(lines) + "\n"


def parse_arguments(args_list):
    parser = UpstreamArgumentParser(
        prog="TExTra upstream",
        description="Run prep + qual + quant in one command.",
        add_help=True,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True, help="Input TSV used by prep and quant read resolution.")
    parser.add_argument("-o", "--out_dir", required=True, help="Output root directory shared by prep, qual, and quant.")
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA used by prep and quant reference/index construction.")
    parser.add_argument("-G", "--gene", required=True, help="Reference gene annotation GTF used by prep.")
    parser.add_argument("-r", "--te", required=True, help="TE annotation file used by prep; supports GTF/GFF/BED and RepeatMasker OUT/TXT.")
    parser.add_argument("-s", "--samples", default=None, help="Sample/condition names, comma-separated; if omitted, infer from input first column.")
    add_threading_arguments(parser, threads_help="Threads per module.")
    add_project_argument(parser, default="project", help_text="Project prefix.")
    add_read_layout_arguments(parser, strand_help="Library strand type.")
    add_run_mode_arguments(parser, debug_help="Enable debug mode and keep intermediate files.")

    add_prep_specific_arguments(parser)

    parser.add_argument(
        "--calculate-afe-ale",
        action="store_true",
        default=False,
        help="Enable AFE/ALE usage outputs in qual. Requires HITindex; cannot be combined with --skip-hitindex.",
    )
    parser.add_argument(
        "--skip-hitindex",
        action="store_true",
        help=(
            "Skip HITindex positional classification in qual and produce annotation/TE-overlap tables only. "
            "Cannot be combined with --calculate-afe-ale; junction evidence is disabled."
        ),
    )
    parser.add_argument("--genmodel-iters", type=int, default=100000, help="HITindex ADVI iterations.")
    parser.add_argument("--bootstrap-n", type=int, default=1000, help="Bootstrap iterations.")
    parser.add_argument("--seed", type=int, default=None, help="Optional random seed passed to qual HITindex steps.")
    parser.add_argument("--ss3-buffer", dest="ss3buffer", type=int, default=20, help="3' splice-site buffer size for HITindex splice-site roles.")
    parser.add_argument("--ss5-buffer", dest="ss5buffer", type=int, default=50, help="5' splice-site buffer size for HITindex splice-site roles.")
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

    parser.add_argument("--quantifier", choices=["rsem", "salmon"], default="rsem", help="Quantification backend for quant.")
    parser.add_argument(
        "--quant-result-dir",
        default=None,
        help="Directory containing reusable RSEM/Salmon backend outputs, typically 04_quantification from a debug quant run; final usage tables are not reused.",
    )
    parser.add_argument("--compute-gene-abundance", action="store_true", default=False, help="Compute and export gene-level abundance tables in quant.")
    return parser.parse_args(args_list)


def _validate_parameters(args):
    validate_threading_args(args)
    validate_read_layout(args, module_name="upstream")
    if bool(args.calculate_afe_ale) and bool(args.skip_hitindex):
        raise RuntimeError("--calculate-afe-ale requires HITindex and cannot be combined with --skip-hitindex.")
    if int(args.genmodel_iters) < 1:
        raise RuntimeError("--genmodel-iters must be a positive integer.")
    if int(args.bootstrap_n) < 1:
        raise RuntimeError("--bootstrap-n must be a positive integer.")
    if int(args.ss3buffer) < 0 or int(args.ss5buffer) < 0:
        raise RuntimeError("--ss3-buffer and --ss5-buffer must be non-negative integers.")
    if int(args.te_overlap_min_bp) < 1:
        raise RuntimeError("--te-overlap-min-bp must be a positive integer.")
    if not (0.0 <= float(args.te_overlap_min_frac) <= 1.0):
        raise RuntimeError("--te-overlap-min-frac must be between 0 and 1.")
    if int(args.splice_site_flank_bp) < 0:
        raise RuntimeError("--splice-site-flank-bp must be a non-negative integer.")
    if args.min_expr is not None and float(args.min_expr) < 0:
        raise RuntimeError("--min-expr must be non-negative when provided.")
    if args.min_length is not None and int(args.min_length) < 1:
        raise RuntimeError("--min-length must be a positive integer when provided.")
    if args.min_frac is not None and not (0.0 <= float(args.min_frac) <= 1.0):
        raise RuntimeError("--min-frac must be between 0 and 1 when provided.")


def _infer_samples(input_tsv):
    samples = []
    seen = set()
    with open(input_tsv, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            sample = line.rstrip("\n").split("\t", 1)[0].strip()
            if not sample:
                continue
            if sample in {"NA", "NaN", "nan"}:
                continue
            if sample not in seen:
                seen.add(sample)
                samples.append(sample)
    if not samples:
        raise RuntimeError(f"No sample names found in input TSV first column: {input_tsv}")
    return ",".join(samples)


def main(args_list=None):
    if args_list is None:
        args_list = sys.argv[2:] if len(sys.argv) > 1 and sys.argv[1] == "upstream" else sys.argv[1:]
    args = parse_arguments(args_list)
    _validate_parameters(args)
    from TExTra.bin.mode0_pipeline import main as prep_main
    from TExTra.bin.mode1_pipeline import main as qual_main
    from TExTra.bin.mode2_pipeline import main as quant_main

    samples = args.samples or _infer_samples(args.input)

    prep_args = [
        "-i", args.input,
        "-o", args.out_dir,
        "-g", args.genome,
        "-G", args.gene,
        "-r", args.te,
        "-t", str(args.threads),
        "--strand", args.strand,
        "--readtype", args.readtype,
    ]
    prep_args.extend(build_prep_passthrough_args(args))
    if args.njobs is not None:
        prep_args.extend(["--njobs", str(args.njobs)])

    qual_args = [
        "--prep", args.out_dir,
        "-o", args.out_dir,
        "-s", samples,
        "-t", str(args.threads),
        "--project", args.project,
        "--strand", args.strand,
        "--readtype", args.readtype,
        "--ss3-buffer", str(args.ss3buffer),
        "--ss5-buffer", str(args.ss5buffer),
        "--genmodel-iters", str(args.genmodel_iters),
        "--bootstrap-n", str(args.bootstrap_n),
        "--te-overlap-min-bp", str(args.te_overlap_min_bp),
        "--te-overlap-min-frac", str(args.te_overlap_min_frac),
        "--splice-site-flank-bp", str(args.splice_site_flank_bp),
    ]
    if args.njobs is not None:
        qual_args.extend(["--njobs", str(args.njobs)])
    if args.seed is not None:
        qual_args.extend(["--seed", str(args.seed)])
    if args.debug:
        qual_args.append("--debug")
    if args.detail:
        qual_args.append("--detail")
    if args.calculate_afe_ale:
        qual_args.append("--calculate-afe-ale")
    if args.skip_hitindex:
        qual_args.append("--skip-hitindex")
    if not args.te_overlap_junction_evidence:
        qual_args.append("--ignore-junction")

    quant_args = [
        "--prep", args.out_dir,
        "--qual", args.out_dir,
        "-o", args.out_dir,
        "-s", samples,
        "-i", args.input,
        "-t", str(args.threads),
        "--quantifier", args.quantifier,
        "--project", args.project,
        "--strand", args.strand,
        "--readtype", args.readtype,
    ]
    quant_args += ["-g", args.genome]
    if args.quant_result_dir:
        quant_args += ["--quant-result-dir", args.quant_result_dir]
    if args.compute_gene_abundance:
        quant_args.append("--compute-gene-abundance")
    if args.njobs is not None:
        quant_args.extend(["--njobs", str(args.njobs)])
    if args.debug:
        quant_args.append("--debug")
    if args.detail:
        quant_args.append("--detail")

    _print_module_separator(1, "prep")
    prep_main(prep_args)
    _print_module_separator(2, "qual")
    qual_main(qual_args)
    _print_module_separator(3, "quant")
    quant_main(quant_args)


if __name__ == "__main__":
    main()

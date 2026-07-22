"""Validate qual inputs, parameters, runtime dependencies, and consistency."""

import importlib.util
import os
import shutil
from collections import Counter

from util.common.define_cli import validate_read_layout, validate_threading_args
from util.common.define_layout import ALIGNMENT_DIR, CONVERT_DIR, resolve_consensus_gtf

def validate_inputs(args):
    """Validate mode1 dependencies from prep output."""
    if not os.path.isdir(args.prep):
        raise FileNotFoundError(f"prep directory not found: {args.prep}")
    if getattr(args, "hitindex_dir", None) and not os.path.isdir(args.hitindex_dir):
        raise FileNotFoundError(f"HITindex reuse directory not found: {args.hitindex_dir}")
    required_paths = [
        os.path.join(args.prep, CONVERT_DIR, "TE_anno.bed"),
        resolve_consensus_gtf(args.prep),
    ]
    if not args.skip_hitindex:
        required_paths.extend(
            [
                os.path.join(args.prep, ALIGNMENT_DIR),
                os.path.join(args.prep, CONVERT_DIR, "gene_anno.bed"),
            ]
        )
    for path in required_paths:
        if not os.path.exists(path):
            raise FileNotFoundError(f"required prep output missing: {path}")

def validate_parameters(args):
    """Validate qual CLI parameter ranges and shared read-layout options."""
    validate_threading_args(args)
    validate_read_layout(args, module_name="qual")
    if bool(getattr(args, "skip_hitindex", False)) and bool(getattr(args, "calculate_afe_ale", False)):
        raise RuntimeError("--calculate-afe-ale requires HITindex and cannot be combined with --skip-hitindex.")
    if bool(getattr(args, "skip_hitindex", False)) and getattr(args, "hitindex_dir", None):
        raise RuntimeError("--hitindex-dir reuses HITindex outputs and cannot be combined with --skip-hitindex.")
    if int(args.ss3buffer) < 0 or int(args.ss5buffer) < 0:
        raise RuntimeError("--ss3-buffer and --ss5-buffer must be non-negative integers.")
    if int(args.genmodel_iters) < 1:
        raise RuntimeError("--genmodel-iters must be a positive integer.")
    if int(args.bootstrap_n) < 1:
        raise RuntimeError("--bootstrap-n must be a positive integer.")
    if int(args.te_overlap_min_bp) < 1:
        raise RuntimeError("--te-overlap-min-bp must be a positive integer.")
    if not (0.0 <= float(args.te_overlap_min_frac) <= 1.0):
        raise RuntimeError("--te-overlap-min-frac must be between 0 and 1.")
    if int(args.splice_site_flank_bp) < 0:
        raise RuntimeError("--splice-site-flank-bp must be a non-negative integer.")

def validate_runtime(args):
    """Fail early with actionable messages for tools imported or executed by qual."""
    missing = []
    required_python_packages = ["filelock", "numpy", "pandas", "pybedtools", "scipy"]
    if not args.skip_hitindex:
        required_python_packages.append("pymc3")
    for package in required_python_packages:
        if importlib.util.find_spec(package) is None:
            missing.append(f"Python package '{package}'")
    if not (shutil.which("bedtools") or (shutil.which("sortBed") and shutil.which("mergeBed"))):
        missing.append("bedtools executables")
    if not args.skip_hitindex:
        for tool in ["samtools", "intersectBed"]:
            if shutil.which(tool) is None:
                missing.append(tool)
    if missing:
        raise RuntimeError(
            "Missing qual runtime dependency/dependencies: "
            + ", ".join(missing)
            + ". Check the TExTra environment and PATH before running qual."
        )

def validate_sample_list(sample_list):
    """Validate qual sample names parsed from --samples."""
    if not sample_list:
        raise RuntimeError("No samples provided to qual. Use --samples with a comma-separated sample list.")
    counts = Counter(sample_list)
    duplicates = sorted([sample for sample, count in counts.items() if count > 1])
    if duplicates:
        raise RuntimeError(
            "Duplicate sample name(s) in --samples are not allowed: "
            + ", ".join(duplicates)
        )

def validate_sample_bams(sample_list, bamfiles_dict, replicates_dict):
    """Validate prep BAM files and replicate labels required by HITindex."""
    missing = []
    empty_reps = []
    missing_files = []
    for sample in sample_list:
        bamfiles = bamfiles_dict.get(sample, [])
        replicates = replicates_dict.get(sample, [])
        if not bamfiles:
            missing.append(sample)
            continue
        if not replicates:
            empty_reps.append(sample)
        for bam in bamfiles:
            if not os.path.isfile(bam):
                missing_files.append(f"{sample}: {bam}")
    if missing:
        raise RuntimeError(
            "No accepted-hit BAM files found for sample(s): "
            + ", ".join(missing)
            + ". Check prep alignment outputs or --samples."
        )
    if empty_reps:
        raise RuntimeError(
            "No replicate labels found for sample(s): "
            + ", ".join(empty_reps)
            + ". Check prep alignment outputs."
        )
    if missing_files:
        raise FileNotFoundError(
            "BAM file(s) referenced by qual are missing: "
            + "; ".join(missing_files)
        )

def validate_gene_bed_matches_consensus(gene_bed_path, transcript_exon_rows):
    """Validate that prep gene_anno.bed matches consensus transcript exon rows."""
    consensus_keys = {
        (
            str(row.get("transcript_id", "")),
            str(row.get("exon_chrom", "")),
            int(row.get("exon_start", -1)),
            int(row.get("exon_end", -1)),
            str(row.get("exon_strand", "")),
        )
        for row in transcript_exon_rows
        if str(row.get("transcript_id", "")).strip()
    }
    if not consensus_keys:
        raise RuntimeError("No consensus transcript exon rows were loaded for qual.")

    bed_keys = set()
    with open(gene_bed_path, "r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise RuntimeError(f"Invalid gene BED at {gene_bed_path}:{line_no}: expected at least 6 columns.")
            chrom, start_s, end_s, name, _score, strand = fields[:6]
            parts = name.split(":")
            if len(parts) != 3:
                raise RuntimeError(
                    f"Invalid gene BED at {gene_bed_path}:{line_no}: name must use gene:transcript:exon format."
                )
            try:
                start0 = int(start_s)
                end0 = int(end_s)
            except ValueError as exc:
                raise RuntimeError(f"Invalid gene BED at {gene_bed_path}:{line_no}: start/end must be integers.") from exc
            bed_keys.add((parts[1], chrom, start0, end0, strand))

    missing = sorted(consensus_keys - bed_keys)
    if missing:
        preview = "; ".join(f"{tx}:{chrom}:{start}-{end}:{strand}" for tx, chrom, start, end, strand in missing[:5])
        raise RuntimeError(
            "gene_anno.bed is inconsistent with consensus_transcripts.gtf: "
            f"{len(missing)} consensus exon(s) are missing from {gene_bed_path}. "
            f"Examples: {preview}. Rerun prep so gene_anno.bed is regenerated from consensus_transcripts.gtf."
        )

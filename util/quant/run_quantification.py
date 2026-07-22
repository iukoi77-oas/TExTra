"""Run quant reference preparation and sample-level expression quantification."""

import os
import shutil

import pandas as pd

from util.common.define_layout import CLASSIFY_DIR, QUANT_DIR, resolve_consensus_gtf
from util.common.run_config import extract_module_config, read_global_config
from util.common.write_logs import log_message
from util.common.collect_files import parse_sample_csv
from util.quant.resolve_inputs import (
    infer_rsem_strandedness as _infer_rsem_strandedness,
    infer_star_path as _infer_star_path,
    is_valid_rsem_reference_prefix as _is_valid_rsem_reference_prefix,
    looks_like_fasta_path as _looks_like_fasta_path,
    looks_like_gzip_path as _looks_like_gzip_path,
    pick_existing_col as _pick_existing_col,
    resolve_rsem_aligner as _resolve_rsem_aligner,
    run_command as _run_command,
)


def _resolve_genome_fasta_arg(args):
    return getattr(args, "genome", None)


def _validate_inputs(args):
    """Validate mode2 dependencies for RSEM novel-transcript quantification."""
    if not os.path.isdir(args.prep):
        raise FileNotFoundError(f"prep directory not found: {args.prep}")
    if not os.path.isdir(args.qual):
        raise FileNotFoundError(f"qual directory not found: {args.qual}")
    if getattr(args, "quant_result_dir", None) and not os.path.isdir(args.quant_result_dir):
        raise FileNotFoundError(f"quant result reuse directory not found: {args.quant_result_dir}")

    annotation_path = os.path.join(args.qual, CLASSIFY_DIR, "transcript_exon_te_annotation.tsv")
    project_annotation_path = os.path.join(args.qual, CLASSIFY_DIR, f"{args.project}.transcript_exon_te_annotation.tsv")
    legacy_annotation_path = os.path.join(
        args.qual,
        CLASSIFY_DIR,
        "annotation",
        f"{args.project}.transcript_exon_te_annotation.tsv",
    )
    required_paths = [
        resolve_consensus_gtf(args.prep),
        annotation_path if os.path.isfile(annotation_path) else (
            project_annotation_path if os.path.isfile(project_annotation_path) else legacy_annotation_path
        ),
    ]
    for path in required_paths:
        if not os.path.exists(path):
            raise FileNotFoundError(f"required upstream output missing: {path}")

    sample_list = parse_sample_csv(args.samples)
    if not sample_list:
        raise RuntimeError("No valid sample names parsed from --samples.")
    if getattr(args, "quant_result_dir", None):
        return

    ref_prefix = args.rsem_reference
    if ref_prefix:
        ref_prefix = os.path.abspath(ref_prefix)
    else:
        ref_prefix = os.path.join(args.out_dir, QUANT_DIR, "rsem_novel_ref", args.project)
    need_build = bool(args.rsem_build_reference) or (not _is_valid_rsem_reference_prefix(ref_prefix))
    if str(getattr(args, "quantifier", "rsem")).lower() == "salmon":
        salmon_index = getattr(args, "salmon_index", None)
        if salmon_index:
            salmon_index = os.path.abspath(salmon_index)
        else:
            salmon_index = os.path.join(args.out_dir, QUANT_DIR, "salmon_novel_index", args.project)
        need_build = bool(getattr(args, "salmon_build_index", False)) or (not os.path.isdir(salmon_index))
    if need_build:
        genome_fa = _resolve_genome_fasta_arg(args)
        if not genome_fa:
            if str(getattr(args, "quantifier", "rsem")).lower() == "salmon":
                raise RuntimeError(
                    "Salmon index is missing and cannot be built because --genome was not provided."
                )
            raise RuntimeError(
                "RSEM reference is missing and cannot be built because --genome was not provided."
            )
        if not os.path.isfile(genome_fa):
            raise FileNotFoundError(f"genome FASTA not found: {genome_fa}")


def _load_qual_te_overlap_config(qual_dir):
    """Load qual TE-overlap provenance parameters for logging."""
    global_config, global_config_path = read_global_config(qual_dir)
    qual_config = extract_module_config(global_config, "qual")
    if qual_config:
        return qual_config

    log_message(
        "[WARNING]",
        f"qual run config not found in {global_config_path}; quant will continue using qual annotation tables.",
        color="warning",
    )
    return {}


def _on_off(value):
    return "on" if bool(value) else "off"


def _log_qual_te_overlap_config(config, detail=False):
    if not config:
        return
    log_message(
        "[INFO]",
        (
            "Qual context: "
            f"HITindex={_on_off(not bool(config.get('skip_hitindex', False)))}, "
            f"junction evidence={_on_off(config.get('te_overlap_junction_evidence_effective', False))}, "
            f"AFE/ALE={_on_off(config.get('calculate_afe_ale', False))}"
        ),
        color="info",
    )
    if not detail:
        return
    keys = [
        "te_overlap_min_bp",
        "te_overlap_min_frac",
        "splice_site_flank_bp",
        "calculate_afe_ale",
        "te_overlap_junction_evidence",
        "te_overlap_junction_evidence_effective",
        "skip_hitindex",
        "seed",
    ]
    summary = ", ".join(f"{key}={config.get(key)}" for key in keys if key in config)
    log_message("[INFO]", f"Detail qual parameters: {summary}", color="info")


def _prepare_rsem_reference(args, quant_dir, transcript_gtf_path):
    rsem_prepare = shutil.which("rsem-prepare-reference")
    if not rsem_prepare:
        raise RuntimeError("rsem-prepare-reference not found in PATH.")

    has_star = bool(shutil.which("STAR"))
    has_bowtie2 = bool(shutil.which("bowtie2"))
    resolved_aligner = _resolve_rsem_aligner(
        aligner=args.rsem_aligner,
        has_star=has_star,
        has_bowtie2=has_bowtie2,
    )
    star_path_dir = _infer_star_path(args.rsem_star_path) if resolved_aligner == "star" else None

    ref_dir = os.path.join(quant_dir, "rsem_novel_ref")
    os.makedirs(ref_dir, exist_ok=True)
    ref_prefix = os.path.abspath(args.rsem_reference) if args.rsem_reference else os.path.join(ref_dir, args.project)
    need_build = bool(args.rsem_build_reference) or (not _is_valid_rsem_reference_prefix(ref_prefix))

    if need_build:
        genome_arg = _resolve_genome_fasta_arg(args)
        if not genome_arg:
            raise RuntimeError("RSEM reference is missing and cannot be built because --genome was not provided.")
        genome_fasta = os.path.abspath(genome_arg)
        prep_log = os.path.join(quant_dir, f"{args.project}.rsem_prepare.log")
        cmd = [
            rsem_prepare,
            "--gtf",
            transcript_gtf_path,
            "--num-threads",
            str(max(int(args.threads), 1)),
        ]
        if resolved_aligner == "star":
            cmd.append("--star")
            if star_path_dir:
                cmd += ["--star-path", star_path_dir]
        elif resolved_aligner == "bowtie2":
            cmd.append("--bowtie2")
        cmd += [genome_fasta, ref_prefix]
        _run_command(cmd, show_tool_output=bool(args.show_tool_output), log_path=prep_log)

    return ref_prefix, resolved_aligner, star_path_dir


def _infer_salmon_libtype(strand, readtype, salmon_libtype="auto"):
    explicit = str(salmon_libtype).strip().upper()
    if explicit and explicit != "AUTO":
        return explicit
    rt = str(readtype).strip().lower()
    st = str(strand).strip().lower()
    if rt == "paired":
        if st in {"rf", "r"}:
            return "ISR"
        if st in {"fr", "f"}:
            return "ISF"
        return "IU"
    if st == "r":
        return "SR"
    if st == "f":
        return "SF"
    return "U"


def _prepare_salmon_index(args, quant_dir, transcript_gtf_path):
    salmon = shutil.which("salmon")
    if not salmon:
        raise RuntimeError("salmon not found in PATH.")
    gffread = shutil.which("gffread")
    if not gffread:
        raise RuntimeError("gffread not found in PATH; required to build Salmon transcript FASTA.")

    ref_dir = os.path.join(quant_dir, "salmon_novel_index")
    os.makedirs(ref_dir, exist_ok=True)
    index_dir = os.path.abspath(args.salmon_index) if getattr(args, "salmon_index", None) else os.path.join(ref_dir, args.project)
    transcript_fa = os.path.join(ref_dir, f"{args.project}.transcripts.fa")
    need_build = bool(getattr(args, "salmon_build_index", False)) or (not os.path.isdir(index_dir))

    if need_build:
        genome_arg = _resolve_genome_fasta_arg(args)
        if not genome_arg:
            raise RuntimeError("Salmon index is missing and cannot be built because --genome was not provided.")
        genome_fasta = os.path.abspath(genome_arg)
        fasta_log = os.path.join(quant_dir, f"{args.project}.salmon_transcripts.log")
        _run_command(
            [gffread, transcript_gtf_path, "-g", genome_fasta, "-w", transcript_fa],
            show_tool_output=bool(args.show_tool_output),
            log_path=fasta_log,
        )
        index_log = os.path.join(quant_dir, f"{args.project}.salmon_index.log")
        _run_command(
            [salmon, "index", "-t", transcript_fa, "-i", index_dir, "-p", str(max(int(args.threads), 1))],
            show_tool_output=bool(args.show_tool_output),
            log_path=index_log,
        )

    if not os.path.isdir(index_dir):
        raise RuntimeError(f"Salmon index directory not found after build: {index_dir}")
    if not os.path.isfile(transcript_fa):
        raise RuntimeError(f"Salmon transcript FASTA not found: {transcript_fa}")
    return index_dir, transcript_fa


def _same_path(path_a, path_b):
    try:
        return os.path.abspath(path_a) == os.path.abspath(path_b)
    except (TypeError, ValueError):
        return False


def _copy_file_if_needed(src, dst):
    if _same_path(src, dst):
        return False
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)
    return True


def _reusable_rsem_paths(args, sample):
    source_dir = getattr(args, "quant_result_dir", None)
    if not source_dir:
        return None
    candidates = [
        (
            os.path.join(source_dir, "RSEM", "result", f"{sample}.isoforms.result"),
            os.path.join(source_dir, "RSEM", "result", f"{sample}.genes.result"),
        ),
        (
            os.path.join(source_dir, f"{sample}.rsem_novel.isoforms.results"),
            os.path.join(source_dir, f"{sample}.rsem_novel.genes.results"),
        ),
    ]
    for src_isoforms, src_genes in candidates:
        if os.path.isfile(src_isoforms):
            return src_isoforms, src_genes if os.path.isfile(src_genes) else None
    return None


def _has_reusable_rsem_result(args, sample, require_gene=False):
    paths = _reusable_rsem_paths(args, sample)
    if paths is None:
        return False
    _src_isoforms, src_genes = paths
    return (not require_gene) or bool(src_genes and os.path.isfile(src_genes))


def _seed_reusable_rsem_result(args, sample, quant_dir):
    """Copy reusable RSEM result files for one sample into the active quant directory."""
    paths = _reusable_rsem_paths(args, sample)
    if paths is None:
        return False
    src_isoforms, src_genes = paths

    out_prefix = os.path.join(quant_dir, f"{sample}.rsem_novel")
    copied = False
    copied |= _copy_file_if_needed(src_isoforms, f"{out_prefix}.isoforms.results")
    if src_genes:
        copied |= _copy_file_if_needed(src_genes, f"{out_prefix}.genes.results")
    if copied:
        log_message("[INFO]", f"RSEM quant result reused: {sample}", color="info")
    return True


def _reusable_salmon_path(args, sample):
    source_dir = getattr(args, "quant_result_dir", None)
    if not source_dir:
        return None
    src_quant_sf = os.path.join(source_dir, f"{sample}.salmon", "quant.sf")
    if os.path.isfile(src_quant_sf):
        return src_quant_sf
    organized = os.path.join(source_dir, "salmon", "result", f"{sample}.quant.sf")
    if os.path.isfile(organized):
        return organized
    exported = os.path.join(source_dir, f"{sample}_quant.sf")
    if os.path.isfile(exported):
        return exported
    return None


def _has_reusable_salmon_result(args, sample):
    return _reusable_salmon_path(args, sample) is not None


def has_reusable_quant_result(args, sample, quantifier, require_gene=False):
    if str(quantifier).strip().lower() == "salmon":
        return _has_reusable_salmon_result(args, sample)
    return _has_reusable_rsem_result(args, sample, require_gene=require_gene)


def _seed_reusable_salmon_result(args, sample, quant_dir):
    """Copy reusable Salmon quant.sf for one sample into the active quant directory."""
    src_quant_sf = _reusable_salmon_path(args, sample)
    if src_quant_sf is None:
        return False

    out_dir = os.path.join(quant_dir, f"{sample}.salmon")
    dst_quant_sf = os.path.join(out_dir, "quant.sf")
    copied = _copy_file_if_needed(src_quant_sf, dst_quant_sf)
    if copied:
        log_message("[INFO]", f"Salmon quant result reused: {sample}", color="info")
    return True


def _run_rsem_for_sample(
    args,
    sample,
    reads_info,
    ref_prefix,
    resolved_aligner,
    star_path_dir,
    quant_dir,
    load_gene_abundance=False,
):
    out_prefix = os.path.join(quant_dir, f"{sample}.rsem_novel")
    isoforms_tsv = f"{out_prefix}.isoforms.results"
    genes_tsv = f"{out_prefix}.genes.results"
    _seed_reusable_rsem_result(args, sample, quant_dir)
    need_run = bool(args.force) or (not os.path.isfile(isoforms_tsv)) or (
        bool(load_gene_abundance) and not os.path.isfile(genes_tsv)
    )

    if need_run:
        rsem_calc = shutil.which("rsem-calculate-expression")
        if not rsem_calc:
            raise RuntimeError("rsem-calculate-expression not found in PATH.")
        calc_log = os.path.join(quant_dir, f"{sample}.rsem_calc.log")
        cmd = [rsem_calc, "--num-threads", str(max(int(args.threads), 1))]
        rs = _infer_rsem_strandedness(args.strand, args.rsem_strandedness)
        if rs in {"forward", "reverse"}:
            cmd += ["--strandedness", rs]

        if resolved_aligner == "star":
            cmd.append("--star")
            if star_path_dir:
                cmd += ["--star-path", star_path_dir]
        elif resolved_aligner == "bowtie2":
            cmd.append("--bowtie2")

        if str(args.readtype).lower() == "paired":
            r1 = os.path.abspath(reads_info["read1"])
            r2 = os.path.abspath(reads_info["read2"])
            if not os.path.isfile(r1):
                raise FileNotFoundError(f"read1 not found: {r1}")
            if not os.path.isfile(r2):
                raise FileNotFoundError(f"read2 not found: {r2}")
            cmd.append("--paired-end")
            auto_noq = _looks_like_fasta_path(r1) or _looks_like_fasta_path(r2)
            if bool(args.rsem_no_qualities) or auto_noq:
                cmd.append("--no-qualities")
            if resolved_aligner == "star" and (_looks_like_gzip_path(r1) or _looks_like_gzip_path(r2)):
                cmd.append("--star-gzipped-read-file")
            cmd += [r1, r2, ref_prefix, out_prefix]
        else:
            single_reads = [os.path.abspath(p) for p in list(reads_info["single"] or []) if str(p).strip()]
            if not single_reads:
                raise RuntimeError(f"No single-end reads resolved for sample: {sample}")
            for p in single_reads:
                if not os.path.isfile(p):
                    raise FileNotFoundError(f"single-end read not found: {p}")
            auto_noq = any(_looks_like_fasta_path(p) for p in single_reads)
            if bool(args.rsem_no_qualities) or auto_noq:
                cmd.append("--no-qualities")
            if resolved_aligner == "star" and any(_looks_like_gzip_path(p) for p in single_reads):
                cmd.append("--star-gzipped-read-file")
            cmd += single_reads + [ref_prefix, out_prefix]

        _run_command(cmd, show_tool_output=bool(args.show_tool_output), log_path=calc_log)

    if not os.path.isfile(isoforms_tsv):
        raise RuntimeError(f"RSEM isoforms output not found: {isoforms_tsv}")
    iso_df = pd.read_csv(isoforms_tsv, sep="\t")
    gene_df = None
    if load_gene_abundance:
        if not os.path.isfile(genes_tsv):
            raise RuntimeError(f"RSEM genes output not found: {genes_tsv}")
        gene_df = pd.read_csv(genes_tsv, sep="\t")
    return iso_df, gene_df, isoforms_tsv, genes_tsv


def _run_salmon_for_sample(
    args,
    sample,
    reads_info,
    salmon_index_dir,
    quant_dir,
):
    out_dir = os.path.join(quant_dir, f"{sample}.salmon")
    quant_sf = os.path.join(out_dir, "quant.sf")
    exported_quant_sf = os.path.join(quant_dir, f"{sample}_quant.sf")
    _seed_reusable_salmon_result(args, sample, quant_dir)
    need_run = bool(args.force) or (not os.path.isfile(quant_sf))
    if need_run:
        salmon = shutil.which("salmon")
        if not salmon:
            raise RuntimeError("salmon not found in PATH.")
        os.makedirs(out_dir, exist_ok=True)
        quant_log = os.path.join(quant_dir, f"{sample}.salmon_quant.log")
        cmd = [
            salmon,
            "quant",
            "-i",
            salmon_index_dir,
            "-l",
            _infer_salmon_libtype(args.strand, args.readtype, getattr(args, "salmon_libtype", "auto")),
            "-p",
            str(max(int(args.threads), 1)),
            "--validateMappings",
            "-o",
            out_dir,
        ]
        if str(args.readtype).lower() == "paired":
            r1 = os.path.abspath(reads_info["read1"])
            r2 = os.path.abspath(reads_info["read2"])
            if not os.path.isfile(r1):
                raise FileNotFoundError(f"read1 not found: {r1}")
            if not os.path.isfile(r2):
                raise FileNotFoundError(f"read2 not found: {r2}")
            cmd += ["-1", r1, "-2", r2]
        else:
            single_reads = [os.path.abspath(p) for p in list(reads_info["single"] or []) if str(p).strip()]
            if not single_reads:
                raise RuntimeError(f"No single-end reads resolved for sample: {sample}")
            for p in single_reads:
                if not os.path.isfile(p):
                    raise FileNotFoundError(f"single-end read not found: {p}")
            cmd += ["-r", " ".join(single_reads)] if len(single_reads) == 1 else ["-r"] + single_reads
        _run_command(cmd, show_tool_output=bool(args.show_tool_output), log_path=quant_log)

    if not os.path.isfile(quant_sf):
        raise RuntimeError(f"Salmon quant.sf not found: {quant_sf}")
    shutil.copyfile(quant_sf, exported_quant_sf)
    quant_df = pd.read_csv(quant_sf, sep="\t")
    return quant_df, exported_quant_sf


def _require_rsem_result_col(df, candidates, label, source_name):
    col = _pick_existing_col(df, candidates)
    if col is None:
        raise RuntimeError(
            f"Invalid {source_name}: required column for {label} not found. "
            f"Accepted names: {', '.join(candidates)}"
        )
    return col


def _normalize_transcript_quant_df(tool_name, sample, quant_df, tx_to_gene, te_overlap_exon_txs):
    source_name = f"{tool_name} transcript output"
    if str(tool_name).lower() == "salmon":
        tx_col = _require_rsem_result_col(quant_df, ["Name"], "transcript identifier", source_name)
        abundance_col = _require_rsem_result_col(quant_df, ["NumReads"], "NumReads", source_name)
        tpm_col = _require_rsem_result_col(quant_df, ["TPM"], "TPM", source_name)
        local = quant_df.copy()
        local["_txid"] = local[tx_col].astype(str)
        local["_abundance"] = pd.to_numeric(local[abundance_col], errors="coerce").fillna(0.0)
        local["_tpm"] = pd.to_numeric(local[tpm_col], errors="coerce").fillna(0.0)
        local["_gene"] = local["_txid"].map(lambda x: str(tx_to_gene.get(x, "")))
    else:
        tx_col = _require_rsem_result_col(
            quant_df,
            ["transcript_id", "transcript", "name"],
            "transcript identifier",
            "RSEM isoforms.results",
        )
        expected_col = _require_rsem_result_col(
            quant_df,
            ["expected_count", "expected_counts", "numreads"],
            "expected count",
            "RSEM isoforms.results",
        )
        tpm_col = _require_rsem_result_col(
            quant_df,
            ["tpm"],
            "TPM",
            "RSEM isoforms.results",
        )
        gene_col = _pick_existing_col(quant_df, ["gene_id", "gene"])
        local = quant_df.copy()
        local["_txid"] = local[tx_col].astype(str)
        local["_abundance"] = pd.to_numeric(local[expected_col], errors="coerce").fillna(0.0)
        local["_tpm"] = pd.to_numeric(local[tpm_col], errors="coerce").fillna(0.0)
        if gene_col:
            rsem_gene = local[gene_col].astype(str)
        else:
            rsem_gene = pd.Series([""] * len(local), index=local.index)
        local["_gene"] = local["_txid"].map(lambda x: str(tx_to_gene.get(x, "")))
        local["_gene"] = local["_gene"].where(local["_gene"].astype(str).str.strip() != "", rsem_gene)

    return pd.DataFrame(
        {
            "sample": sample,
            "transcript_id": local["_txid"].astype(str),
            "gene_id": local["_gene"].astype(str),
            "te_overlap_exon_transcript": local["_txid"].map(
                lambda x: "yes" if x in te_overlap_exon_txs else "no"
            ),
            "estimated_count": local["_abundance"].astype(float),
            "tpm": local["_tpm"].astype(float),
        }
    )

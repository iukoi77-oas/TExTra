"""Quant-stage orchestration for transcript abundance and TE-exon usage."""

import copy
import glob
import os
import shutil
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from util.common.define_layout import ALIGNMENT_DIR, QUANT_DIR, resolve_consensus_gtf
from util.common.write_logs import log_message, set_log_file
from util.common.collect_files import collect_bam_and_replicates, parse_sample_csv
from util.quant.resolve_inputs import (
    load_transcript_gene_map as _load_transcript_gene_map,
    load_prep_input_reads as _load_prep_input_reads,
    pick_existing_col as _pick_existing_col,
    resolve_sample_reads as _resolve_sample_reads,
)
from util.quant.adapt_qual_outputs import (
    _load_exon_event_annotation_table,
    _load_hitindex_position_tables,
    _load_te_exon_event_tx_map,
    _load_transcript_structure_position_table,
)
from util.quant.run_quantification import (
    _load_qual_te_overlap_config,
    _log_qual_te_overlap_config,
    _normalize_transcript_quant_df,
    _prepare_rsem_reference,
    _prepare_salmon_index,
    _require_rsem_result_col,
    _resolve_genome_fasta_arg,
    _run_rsem_for_sample,
    _run_salmon_for_sample,
    _validate_inputs,
    has_reusable_quant_result as _has_reusable_quant_result,
)
from util.quant.write_usage_tables import (
    _build_exon_event_usage_and_matrix,
    _cleanup_quant_outputs,
    _log_missing_supporting_transcripts,
    _write_project_exon_usage_tables,
)


def _on_off(value):
    return "on" if bool(value) else "off"


def _quant_read_source(args):
    if getattr(args, "input", None):
        return "input TSV"
    if getattr(args, "rsem_auto_bam2fastq", False):
        return "BAM fallback"
    return "unresolved"


def _log_quant_parameters(args, sample_n, quant_sample_n):
    genome = _resolve_genome_fasta_arg(args) or "not_set"
    log_message(
        "[INFO]",
        (
            "Quant parameters: "
            f"sample_groups={sample_n}, samples={quant_sample_n}, "
            f"quantifier={str(getattr(args, 'quantifier', 'rsem')).lower()}, "
            f"read source={_quant_read_source(args)}, "
            f"readtype={args.readtype}, strand={args.strand}, "
            f"genome={genome}, "
            f"gene abundance={_on_off(getattr(args, 'compute_gene_abundance', False))}, "
            f"quant result reuse={_on_off(bool(getattr(args, 'quant_result_dir', None)))}, "
            f"reuse context={_on_off(getattr(args, 'reuse', False))}"
        ),
        color="info",
    )


def _format_count_pct(count, total):
    count = int(count)
    total = int(total)
    if total <= 0:
        return f"{count} (NA)"
    return f"{count} ({count / total * 100:.1f}%)"


def _format_elapsed(seconds):
    seconds = float(seconds)
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes = seconds / 60.0
    if minutes < 60:
        return f"{minutes:.1f}min"
    return f"{minutes / 60.0:.1f}h"


def _format_value_counts(series, max_items=6):
    if series is None:
        return ""
    values = pd.Series(series).fillna("NA").astype(str).str.strip()
    values = values.mask(values.eq(""), "NA")
    counts = values.value_counts(dropna=False)
    if counts.empty:
        return ""
    items = [f"{label}={int(count)}" for label, count in counts.head(max_items).items()]
    if len(counts) > max_items:
        items.append(f"other={int(counts.iloc[max_items:].sum())}")
    return ", ".join(items)


def _append_log_section(handle, title, path):
    if not path or not os.path.isfile(path):
        return False
    abs_path = os.path.abspath(path)
    handle.write(f"========={title}=========\n")
    handle.write(f"source: {abs_path}\n\n")
    with open(path, "r", encoding="utf-8", errors="replace") as source:
        content = source.read()
    handle.write(content)
    if content and not content.endswith("\n"):
        handle.write("\n")
    handle.write("\n")
    return True


def _write_quant_extend_log(args, quant_dir, sample_list, quantifier):
    logs_dir = os.path.join(args.out_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    extend_log = os.path.join(logs_dir, "quant_extend.log")
    wrote_any = False

    with open(extend_log, "w", encoding="utf-8") as handle:
        project = args.project
        if quantifier == "salmon":
            wrote_any |= _append_log_section(
                handle,
                "SALMON-TRANSCRIPTS",
                os.path.join(quant_dir, f"{project}.salmon_transcripts.log"),
            )
            wrote_any |= _append_log_section(
                handle,
                "SALMON-INDEX",
                os.path.join(quant_dir, f"{project}.salmon_index.log"),
            )
            for sample in sample_list:
                wrote_any |= _append_log_section(
                    handle,
                    f"SALMON-{sample}",
                    os.path.join(quant_dir, f"{sample}.salmon_quant.log"),
                )
        else:
            wrote_any |= _append_log_section(
                handle,
                "RSEM-PREPARE",
                os.path.join(quant_dir, f"{project}.rsem_prepare.log"),
            )
            for sample in sample_list:
                wrote_any |= _append_log_section(
                    handle,
                    f"RSEM-{sample}",
                    os.path.join(quant_dir, f"{sample}.rsem_calc.log"),
                )
                wrote_any |= _append_log_section(
                    handle,
                    f"RSEM-{sample}-BACKEND",
                    os.path.join(quant_dir, f"{sample}.rsem_novel.log"),
                )

        bam2fastq_root = getattr(args, "rsem_bam2fastq_dir", None) or os.path.join(
            args.out_dir,
            "tmp",
            "rsem_bam2fastq",
        )
        for path in sorted(glob.glob(os.path.join(bam2fastq_root, "**", "*.bam2fastq.log"), recursive=True)):
            sample_name = os.path.basename(path).replace(".bam2fastq.log", "")
            wrote_any |= _append_log_section(handle, f"BAM2FASTQ-{sample_name}", path)

        if not wrote_any:
            handle.write("No backend command logs were found for this quant run.\n")

    return extend_log


def _copy_file_if_present(src, dst):
    if not os.path.isfile(src):
        return False
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)
    return True


def _copy_dir_contents(src_dir, dst_dir, skip_suffixes=()):
    if not os.path.isdir(src_dir):
        return False
    copied = False
    for root, dirs, files in os.walk(src_dir):
        rel_root = os.path.relpath(root, src_dir)
        out_root = dst_dir if rel_root == "." else os.path.join(dst_dir, rel_root)
        os.makedirs(out_root, exist_ok=True)
        for name in files:
            if any(name.endswith(suffix) for suffix in skip_suffixes):
                continue
            src = os.path.join(root, name)
            dst = os.path.join(out_root, name)
            shutil.copy2(src, dst)
            copied = True
        for name in dirs:
            os.makedirs(os.path.join(out_root, name), exist_ok=True)
    return copied


def _organize_quant_backend_debug_outputs(args, quant_dir, sample_list, quantifier, keep_gene_abundance=False):
    backend_dirname = "salmon" if quantifier == "salmon" else "RSEM"
    backend_dir = os.path.join(quant_dir, backend_dirname)
    index_dir = os.path.join(backend_dir, "index")
    result_dir = os.path.join(backend_dir, "result")
    os.makedirs(index_dir, exist_ok=True)
    os.makedirs(result_dir, exist_ok=True)

    if quantifier == "salmon":
        _copy_dir_contents(os.path.join(quant_dir, "salmon_novel_index"), index_dir, skip_suffixes=(".log", ".stat", "Log.out"))
        transcript_fa = os.path.join(quant_dir, "salmon_novel_index", f"{args.project}.transcripts.fa")
        _copy_file_if_present(transcript_fa, os.path.join(index_dir, f"{args.project}.transcripts.fa"))
        for sample in sample_list:
            quant_sf = os.path.join(quant_dir, f"{sample}.salmon", "quant.sf")
            exported_quant_sf = os.path.join(quant_dir, f"{sample}_quant.sf")
            _copy_file_if_present(quant_sf, os.path.join(result_dir, f"{sample}.quant.sf")) or _copy_file_if_present(
                exported_quant_sf,
                os.path.join(result_dir, f"{sample}.quant.sf"),
            )
    else:
        _copy_dir_contents(os.path.join(quant_dir, "rsem_novel_ref"), index_dir, skip_suffixes=(".log", ".stat", "Log.out"))
        for sample in sample_list:
            prefix = os.path.join(quant_dir, f"{sample}.rsem_novel")
            _copy_file_if_present(f"{prefix}.isoforms.results", os.path.join(result_dir, f"{sample}.isoforms.result"))
            if keep_gene_abundance:
                _copy_file_if_present(f"{prefix}.genes.results", os.path.join(result_dir, f"{sample}.genes.result"))
            for bam_path in sorted(glob.glob(f"{prefix}*.bam")):
                _copy_file_if_present(bam_path, os.path.join(result_dir, os.path.basename(bam_path)))

    bam2fastq_root = getattr(args, "rsem_bam2fastq_dir", None) or os.path.join(
        args.out_dir,
        "tmp",
        "rsem_bam2fastq",
    )
    _copy_dir_contents(
        bam2fastq_root,
        os.path.join(backend_dir, "bam2fastq"),
        skip_suffixes=(".log", ".stat", "Log.out"),
    )
    default_bam2fastq_root = os.path.join(args.out_dir, "tmp", "rsem_bam2fastq")
    if os.path.abspath(bam2fastq_root) == os.path.abspath(default_bam2fastq_root):
        shutil.rmtree(bam2fastq_root, ignore_errors=True)
    log_message("[SUCCESS]", f"Quant backend debug artifacts organized: {backend_dir}", color="success")
    return backend_dir


def count_func(args):
    """Run transcript quantification and TE-derived exon usage projection."""
    _validate_inputs(args)

    set_log_file(os.path.join(args.out_dir, "logs", "quant.log"))
    log_message("[INFO]", "Step 1/3: Quantification setup", bold=True, color="step")
    _log_qual_te_overlap_config(
        _load_qual_te_overlap_config(args.qual),
        detail=bool(getattr(args, "detail", False) or getattr(args, "debug", False)),
    )

    sample_list = parse_sample_csv(args.samples)
    candidate_quant_units = []
    if getattr(args, "quant_result_dir", None):
        if getattr(args, "input", None):
            candidate_quant_units = list(_load_prep_input_reads(args.input, sample_list, args.readtype).keys())
        else:
            candidate_quant_units = list(sample_list)
    quantifier = str(getattr(args, "quantifier", "rsem")).strip().lower()
    compute_gene_abundance = bool(getattr(args, "compute_gene_abundance", False))
    all_quant_results_reusable = bool(candidate_quant_units) and all(
        _has_reusable_quant_result(
            args,
            sample,
            quantifier,
            require_gene=compute_gene_abundance,
        )
        for sample in candidate_quant_units
    )
    if all_quant_results_reusable:
        sample_reads = {sample: {} for sample in candidate_quant_units}
    else:
        alignment_dir = os.path.join(args.prep, ALIGNMENT_DIR)
        bamfiles_dict, _ = collect_bam_and_replicates(alignment_dir, sample_list)
        sample_reads = _resolve_sample_reads(args, sample_list, bamfiles_dict=bamfiles_dict)
    quant_sample_list = list(sample_reads.keys())
    _log_quant_parameters(args, sample_n=len(sample_list), quant_sample_n=len(quant_sample_list))

    quant_dir = os.path.join(args.out_dir, QUANT_DIR)
    os.makedirs(quant_dir, exist_ok=True)

    transcript_gtf_path = resolve_consensus_gtf(args.prep)

    if quantifier == "salmon":
        if all_quant_results_reusable:
            salmon_index_dir = None
            log_message(
                "[INFO]",
                f"Salmon quant result reuse: samples={len(quant_sample_list)}, source={args.quant_result_dir}",
                color="info",
            )
        else:
            log_message("[INFO]", "Preparing Salmon transcript index.", color="info")
            index_started = time.monotonic()
            salmon_index_dir, _salmon_transcript_fa = _prepare_salmon_index(
                args=args,
                quant_dir=quant_dir,
                transcript_gtf_path=transcript_gtf_path,
            )
            log_message(
                "[SUCCESS]",
                f"Salmon index ready: {salmon_index_dir}, elapsed: {_format_elapsed(time.monotonic() - index_started)} finished",
                color="success",
            )
        ref_prefix = None
        resolved_aligner = None
        star_path_dir = None
    else:
        if all_quant_results_reusable:
            ref_prefix = None
            resolved_aligner = None
            star_path_dir = None
            log_message(
                "[INFO]",
                f"RSEM quant result reuse: samples={len(quant_sample_list)}, source={args.quant_result_dir}",
                color="info",
            )
        else:
            log_message("[INFO]", "Preparing RSEM consensus transcript reference.", color="info")
            reference_started = time.monotonic()
            ref_prefix, resolved_aligner, star_path_dir = _prepare_rsem_reference(
                args=args,
                quant_dir=quant_dir,
                transcript_gtf_path=transcript_gtf_path,
            )
            log_message(
                "[SUCCESS]",
                f"RSEM reference ready: {ref_prefix}, elapsed: {_format_elapsed(time.monotonic() - reference_started)} finished",
                color="success",
            )
        salmon_index_dir = None

    log_message("[INFO]", "Step 2/3: Transcript quantification", bold=True, color="step")
    tx_to_gene = _load_transcript_gene_map(args, transcript_gtf_path)
    exon_events, event_to_txs, event_to_gene = _load_te_exon_event_tx_map(args.qual, args.project)
    log_message(
        "[INFO]",
        f"Quantification universe: TE-overlap exon events={len(exon_events)}.",
        color="info",
    )
    te_overlap_exon_txs = set()
    for txs in event_to_txs.values():
        te_overlap_exon_txs.update(txs)

    all_tx_frames = []
    all_gene_frames = []
    max_jobs = int(args.njobs or args.threads or 1)
    worker_jobs = min(max_jobs, max(1, len(quant_sample_list)))
    per_worker_threads = max(1, int(args.threads) // worker_jobs)

    total_quant_samples = len(quant_sample_list)

    def _process_sample(sample_index, sample):
        local_args = copy.copy(args)
        local_args.threads = per_worker_threads
        sample_started = time.monotonic()
        status = (
            "reused"
            if getattr(args, "quant_result_dir", None)
            and _has_reusable_quant_result(
                args,
                sample,
                quantifier,
                require_gene=compute_gene_abundance,
            )
            else "computed"
        )
        log_message(
            "[INFO]",
            f"{quantifier.upper()} sample running [{sample_index}/{total_quant_samples}]: {sample}",
            color="info",
        )
        try:
            if quantifier == "salmon":
                quant_df, _exported_quant_sf = _run_salmon_for_sample(
                    args=local_args,
                    sample=sample,
                    reads_info=sample_reads[sample],
                    salmon_index_dir=salmon_index_dir,
                    quant_dir=quant_dir,
                )
                tx_df = _normalize_transcript_quant_df("salmon", sample, quant_df, tx_to_gene, te_overlap_exon_txs)
            else:
                iso_df, genes_df, _iso_path, _gene_path = _run_rsem_for_sample(
                    args=local_args,
                    sample=sample,
                    reads_info=sample_reads[sample],
                    ref_prefix=ref_prefix,
                    resolved_aligner=resolved_aligner,
                    star_path_dir=star_path_dir,
                    quant_dir=quant_dir,
                    load_gene_abundance=compute_gene_abundance,
                )
                tx_df = _normalize_transcript_quant_df("rsem", sample, iso_df, tx_to_gene, te_overlap_exon_txs)
            genes_df = None if quantifier == "salmon" else genes_df

            gene_df = None
            if compute_gene_abundance:
                if quantifier == "salmon":
                    gene_df = tx_df.groupby(["sample", "gene_id"], as_index=False)[["estimated_count", "tpm"]].sum()
                else:
                    g_col = _pick_existing_col(genes_df, ["gene_id", "gene", "name"])
                    g_expected = _require_rsem_result_col(
                        genes_df,
                        ["expected_count", "expected_counts", "numreads"],
                        "expected count",
                        "RSEM genes.results",
                    )
                    g_tpm = _require_rsem_result_col(
                        genes_df,
                        ["tpm"],
                        "TPM",
                        "RSEM genes.results",
                    )
                    if g_col:
                        g_local = genes_df.copy()
                        g_local["_gene"] = g_local[g_col].astype(str)
                        g_local["_abundance"] = pd.to_numeric(g_local[g_expected], errors="coerce").fillna(0.0)
                        g_local["_tpm"] = pd.to_numeric(g_local[g_tpm], errors="coerce").fillna(0.0)
                        gene_df = pd.DataFrame(
                            {
                                "sample": sample,
                                "gene_id": g_local["_gene"].astype(str),
                                "estimated_count": g_local["_abundance"].astype(float),
                                "tpm": g_local["_tpm"].astype(float),
                            }
                        )
                    else:
                        gene_df = tx_df.groupby(["sample", "gene_id"], as_index=False)[["estimated_count", "tpm"]].sum()
        except Exception as exc:
            elapsed = time.monotonic() - sample_started
            log_message(
                "[ERROR]",
                (
                    f"{quantifier.upper()} sample failed [{sample_index}/{total_quant_samples}]: "
                    f"{sample}, elapsed: {_format_elapsed(elapsed)}, error={exc}"
                ),
                color="error",
            )
            raise
        elapsed = time.monotonic() - sample_started
        return sample, status, tx_df, gene_df, elapsed

    log_message(
        "[INFO]",
        (
            f"Running {quantifier.upper()} quantification: "
            f"samples={len(quant_sample_list)}, jobs={worker_jobs}, threads/job={per_worker_threads}."
        ),
        color="info",
    )
    sample_results = {}
    quant_started = time.monotonic()
    status_counts = {"computed": 0, "reused": 0}

    def _record_sample_result(done_n, sample, status, tx_df, gene_df, elapsed):
        sample_results[sample] = (tx_df, gene_df)
        status_counts[status] = status_counts.get(status, 0) + 1
        gene_text = f", genes={len(gene_df)}" if gene_df is not None else ""
        action = "reused" if status == "reused" else "completed"
        log_message(
            "[INFO]",
            (
                f"{quantifier.upper()} sample {action} [{done_n}/{total_quant_samples}]: "
                f"{sample}, transcripts={len(tx_df)}{gene_text}, "
                f"elapsed: {_format_elapsed(elapsed)} finished"
            ),
            color="info",
        )

    if worker_jobs == 1:
        for sample_index, sample in enumerate(quant_sample_list, start=1):
            result = _process_sample(sample_index, sample)
            _record_sample_result(sample_index, *result)
    else:
        with ThreadPoolExecutor(max_workers=worker_jobs) as executor:
            future_map = {
                executor.submit(_process_sample, sample_index, sample): sample
                for sample_index, sample in enumerate(quant_sample_list, start=1)
            }
            for done_n, future in enumerate(as_completed(future_map), start=1):
                sample, status, tx_df, gene_df, elapsed = future.result()
                _record_sample_result(done_n, sample, status, tx_df, gene_df, elapsed)
    log_message(
        "[INFO]",
        (
            f"{quantifier.upper()} quantification completed: "
            f"computed={status_counts.get('computed', 0)}, reused={status_counts.get('reused', 0)}, "
            f"samples={total_quant_samples}, elapsed: {_format_elapsed(time.monotonic() - quant_started)} finished"
        ),
        color="info",
    )

    for sample in quant_sample_list:
        tx_df, gene_df = sample_results[sample]
        tx_out = os.path.join(quant_dir, f"{sample}_transcript_abundance.tsv")
        tx_df.to_csv(tx_out, sep="\t", index=False)
        all_tx_frames.append(tx_df)

        if compute_gene_abundance and gene_df is not None:
            gene_out = os.path.join(quant_dir, f"{sample}_gene_abundance.tsv")
            gene_df.to_csv(gene_out, sep="\t", index=False)
            all_gene_frames.append(gene_df)

    transcript_all = pd.concat(all_tx_frames, ignore_index=True) if all_tx_frames else pd.DataFrame()
    transcript_all_out = os.path.join(quant_dir, "transcript_abundance.tsv")
    transcript_all.to_csv(transcript_all_out, sep="\t", index=False)

    log_message("[SUCCESS]", f"Transcript abundance: {transcript_all_out}", color="success")
    if compute_gene_abundance:
        gene_all = pd.concat(all_gene_frames, ignore_index=True) if all_gene_frames else pd.DataFrame()
        gene_all_out = os.path.join(quant_dir, "gene_abundance.tsv")
        gene_all.to_csv(gene_all_out, sep="\t", index=False)
        log_message("[SUCCESS]", f"Gene abundance: {gene_all_out}", color="success")
    else:
        stale_gene_all = os.path.join(quant_dir, "gene_abundance.tsv")
        if os.path.isfile(stale_gene_all):
            os.remove(stale_gene_all)
        for stale_gene_sample in glob.glob(os.path.join(quant_dir, "*_gene_abundance.tsv")):
            try:
                os.remove(stale_gene_sample)
            except FileNotFoundError:
                pass

    log_message("[INFO]", "Step 3/3: TE-overlap exon usage", bold=True, color="step")
    structure_summary_df = _load_transcript_structure_position_table(args.qual, args.project)
    hitindex_summary_df, hitindex_detail_df = _load_hitindex_position_tables(
        args.qual,
        args.project,
        exon_events,
        structure_summary_df=structure_summary_df,
    )
    hitindex_detail_out = os.path.join(quant_dir, f"{args.project}.exon_hitindex_by_sample.tsv")
    hitindex_detail_df.to_csv(hitindex_detail_out, sep="\t", index=False)
    if hitindex_summary_df.empty:
        log_message(
            "[WARNING]",
            "No position summary was available for quant TE-exon events; ID_position columns will be NA.",
            color="warning",
        )
    else:
        source_counts = hitindex_summary_df["ID_position_source"].fillna("NA").astype(str).value_counts().to_dict()
        source_text = ", ".join(f"{key}={value}" for key, value in sorted(source_counts.items()))
        log_message(
            "[INFO]",
            f"Position evidence: events={len(hitindex_summary_df)}, {source_text}.",
            color="info",
        )
    event_anno_df = _load_exon_event_annotation_table(args.qual, args.project, hitindex_summary_df=hitindex_summary_df)
    usage_df, usage_wide = _build_exon_event_usage_and_matrix(
        transcript_df=transcript_all,
        exon_events=exon_events,
        event_to_txs=event_to_txs,
        event_to_gene=event_to_gene,
        tx_to_gene=tx_to_gene,
        event_anno_df=event_anno_df,
    )
    _log_missing_supporting_transcripts(transcript_all, event_to_txs)
    usage_out = os.path.join(quant_dir, "exon_usage_tx.tsv")
    usage_df.to_csv(usage_out, sep="\t", index=False)
    project_table, detail_table = _write_project_exon_usage_tables(
        usage_wide=usage_wide,
        sample_list=quant_sample_list,
        out_dir=quant_dir,
        project=args.project,
        event_anno_df=event_anno_df,
        write_detail=bool(getattr(args, "detail", False) or getattr(args, "debug", False)),
    )
    if getattr(args, "detail", False) or getattr(args, "debug", False):
        label_counts = _format_value_counts(event_anno_df.get("te_overlap_label"))
        position_counts = _format_value_counts(event_anno_df.get("ID_position_summary"))
        log_message(
            "[INFO]",
            (
                "Detail usage summary: "
                f"events={len(event_anno_df)}, "
                f"te_overlap_label={{{label_counts}}}, "
                f"ID_position={{{position_counts}}}"
            ),
            color="info",
        )
        if detail_table:
            log_message("[SUCCESS]", f"Detail TE-overlap exon usage: {detail_table}", color="success")
    log_message("[SUCCESS]", f"TE-overlap exon usage: {project_table}", color="success")

    if args.debug:
        extend_log = _write_quant_extend_log(args, quant_dir, quant_sample_list, quantifier)
        log_message("[SUCCESS]", f"Quant debug log: {extend_log}", color="success")
        backend_dir = _organize_quant_backend_debug_outputs(
            args,
            quant_dir,
            quant_sample_list,
            quantifier,
            keep_gene_abundance=compute_gene_abundance,
        )
        _cleanup_quant_outputs(
            quant_dir,
            args.project,
            keep_gene_abundance=compute_gene_abundance,
            keep_detail=True,
            keep_dirs={os.path.basename(backend_dir)},
        )

    if not args.debug:
        _cleanup_quant_outputs(
            quant_dir,
            args.project,
            keep_gene_abundance=compute_gene_abundance,
            keep_detail=bool(getattr(args, "detail", False)),
        )

"""Mode1 TE-exon identification step."""

import contextlib
import copy
import json
import os
import re
import shutil
import tempfile
import time

import pandas as pd
from joblib import Parallel, delayed

from util.common.define_layout import (
    ALIGNMENT_DIR,
    ASSEMBLY_DIR,
    CLASSIFY_DIR,
    CONVERT_DIR,
    TRANSCRIPT_GENE_ASSIGNMENT_TSV,
    resolve_consensus_gtf,
)
from util.common.run_config import update_global_config
from util.common.write_logs import log_message, set_log_file
from util.common.collect_files import collect_bam_and_replicates, parse_sample_csv
from util.qual.manage_outputs import (
    export_metaexon_exon_transcript_map,
    finalize_te_overlap_support_table,
    inspect_afe_ale_reuse,
    inspect_junction_support_reuse,
    is_complete_hitindex_replicate,
    load_te_overlap_support_records,
    prune_debug_hitindex_outputs,
    publish_classify_dir,
    rebuild_combined_afe_ale_outputs_from_hitindex,
    seed_hitindex_outputs,
)
from util.qual.validate_inputs import (
    validate_gene_bed_matches_consensus,
    validate_inputs,
    validate_parameters,
    validate_runtime,
    validate_sample_bams,
    validate_sample_list,
)


def _on_off(value):
    return "on" if bool(value) else "off"


def _format_elapsed(seconds):
    seconds = float(seconds or 0.0)
    if seconds >= 60:
        return f"{seconds / 60.0:.1f}min"
    return f"{seconds:.1f}s"


def _parse_average_loss(log_path):
    if not log_path or not os.path.isfile(log_path):
        return "Finished [100%]: Average Loss = NA"
    pattern = re.compile(r"Finished \[100%\]: Average Loss =\s*(.+)")
    last = None
    try:
        with open(log_path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                match = pattern.search(line)
                if match:
                    last = "Finished [100%]: Average Loss = " + match.group(1).strip()
    except OSError:
        return "Finished [100%]: Average Loss = NA"
    return last or "Finished [100%]: Average Loss = NA"


def _safe_log_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("._") or "replicate"


def _append_extend_log_section(extend_log_path, section, source_log_path):
    if not extend_log_path or not source_log_path or not os.path.isfile(source_log_path):
        return
    with open(extend_log_path, "a", encoding="utf-8") as out_handle:
        out_handle.write(f"\n========={section}=========\n")
        with open(source_log_path, "r", encoding="utf-8", errors="replace") as in_handle:
            shutil.copyfileobj(in_handle, out_handle)


def _count_te_overlap_rows(df):
    if df is None or df.empty:
        return 0
    for col in ["te_overlap_pass_final", "te_overlap_pass_any", "te_overlap_pass_raw"]:
        if col in df.columns:
            return int(pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int).sum())
    if "te_overlap_label" in df.columns:
        return int(df["te_overlap_label"].fillna("").astype(str).eq("TE_overlap").sum())
    return 0


def _format_count_pct(count, denominator):
    count = int(count or 0)
    denominator = int(denominator or 0)
    if denominator <= 0:
        return f"{count} (NA)"
    return f"{count} ({count / denominator * 100.0:.1f}%)"


def _nonempty_count(series):
    if series is None:
        return 0
    values = pd.Series(series).fillna("").astype(str).str.strip()
    return int((~values.str.lower().isin({"", "0", "false", "none", "na", "nan"})).sum())


def _positive_count(series):
    if series is None:
        return 0
    return int((pd.to_numeric(series, errors="coerce").fillna(0) > 0).sum())


def _format_value_counts(series, *, max_items=6):
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


def _count_tsv_rows(path):
    if not path or not os.path.isfile(path):
        return 0
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            return max(0, sum(1 for _line in handle) - 1)
    except OSError:
        return 0


def _read_tsv_safe(path):
    if not path or not os.path.isfile(path):
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except (OSError, ValueError, pd.errors.EmptyDataError, pd.errors.ParserError):
        return pd.DataFrame()


def _log_detail_qual_summary(args, annotation_df, teexon_path, classify_dir):
    if not (getattr(args, "detail", False) or getattr(args, "debug", False)):
        return

    total_tx_exons = int(len(annotation_df)) if annotation_df is not None else 0
    te_tx_exons = _count_te_overlap_rows(annotation_df)
    repeat_tx_exons = (
        _nonempty_count(annotation_df["te_splice_site_repeat_TE"])
        if annotation_df is not None and "te_splice_site_repeat_TE" in annotation_df.columns
        else 0
    )
    log_message(
        "[INFO]",
        (
            "Detail summary: "
            f"transcript_exons={total_tx_exons}, "
            f"TE_overlap={_format_count_pct(te_tx_exons, total_tx_exons)}, "
            f"splice_site_repeat_TE={_format_count_pct(repeat_tx_exons, total_tx_exons)}"
        ),
        color="info",
    )

    teexon_df = _read_tsv_safe(teexon_path)
    if not teexon_df.empty:
        side_counts = teexon_df.get("te_boundary_side", pd.Series(dtype=str)).fillna("none").astype(str).str.strip()
        side_counts = side_counts.mask(side_counts.eq(""), "none")
        side_5 = int(side_counts.eq("5'").sum())
        side_3 = int(side_counts.eq("3'").sum())
        side_both = int(side_counts.eq("both").sum())
        side_none = int(side_counts.eq("none").sum())
        log_message(
            "[INFO]",
            (
                "Detail TE boundary side: "
                f"5'={side_5}, 3'={side_3}, both={side_both}, none={side_none} "
                f"(denominator: TE-overlap exons={len(teexon_df)})"
            ),
            color="info",
        )

        candidate_counts = _format_value_counts(teexon_df.get("candidate_TE_event"))
        if candidate_counts:
            log_message("[INFO]", f"Detail HITindex candidate events: {candidate_counts}", color="info")
        confidence_counts = _format_value_counts(teexon_df.get("candidate_TE_confidence"))
        if confidence_counts:
            log_message("[INFO]", f"Detail HITindex candidate confidence: {confidence_counts}", color="info")
        position_counts = _format_value_counts(teexon_df.get("ID_position_summary"))
        if position_counts:
            log_message("[INFO]", f"Detail HITindex position summary: {position_counts}", color="info")

    if getattr(args, "te_overlap_junction_evidence", False) and annotation_df is not None:
        raw_te = (
            int(pd.to_numeric(annotation_df["te_overlap_pass_raw"], errors="coerce").fillna(0).astype(int).sum())
            if "te_overlap_pass_raw" in annotation_df.columns
            else te_tx_exons
        )
        supported = (
            int(pd.to_numeric(annotation_df["junction_supported"], errors="coerce").fillna(0).astype(int).sum())
            if "junction_supported" in annotation_df.columns
            else _positive_count(annotation_df.get("junction_support_sample_n"))
        )
        log_message(
            "[INFO]",
            (
                "Detail junction evidence: "
                f"raw_TE_overlap={raw_te}, retained={te_tx_exons}, "
                f"downgraded={max(0, raw_te - te_tx_exons)}, supported={supported}"
            ),
            color="info",
        )

    if getattr(args, "calculate_afe_ale", False):
        afe_rows = _count_tsv_rows(os.path.join(classify_dir, f"{args.project}.AFEPSI.detail.tsv"))
        ale_rows = _count_tsv_rows(os.path.join(classify_dir, f"{args.project}.ALEPSI.detail.tsv"))
        te_afe_rows = _count_tsv_rows(os.path.join(classify_dir, f"{args.project}.TE_overlap.AFEPSI.tsv"))
        te_ale_rows = _count_tsv_rows(os.path.join(classify_dir, f"{args.project}.TE_overlap.ALEPSI.tsv"))
        log_message(
            "[INFO]",
            (
                "Detail AFE/ALE outputs: "
                f"AFE_rows={afe_rows}, TE_overlap_AFE_rows={te_afe_rows}, "
                f"ALE_rows={ale_rows}, TE_overlap_ALE_rows={te_ale_rows}"
            ),
            color="info",
        )


@contextlib.contextmanager
def _redirect_worker_output(log_path, section):
    """Redirect noisy worker stdout/stderr away from the terminal."""
    target_path = log_path or os.devnull
    os.makedirs(os.path.dirname(target_path), exist_ok=True) if log_path else None
    with open(target_path, "a", encoding="utf-8") as handle:
        if log_path and section:
            handle.write(f"\n========={section}=========\n")
            handle.flush()
        saved_stdout = os.dup(1)
        saved_stderr = os.dup(2)
        try:
            os.dup2(handle.fileno(), 1)
            os.dup2(handle.fileno(), 2)
            yield
        finally:
            os.dup2(saved_stdout, 1)
            os.dup2(saved_stderr, 2)
            os.close(saved_stdout)
            os.close(saved_stderr)


def _log_qual_parameters(args, sample_n, replicate_n):
    junction_requested = bool(
        getattr(args, "_te_overlap_junction_evidence_requested", args.te_overlap_junction_evidence)
    )
    junction_effective = bool(args.te_overlap_junction_evidence and not args.skip_hitindex)
    if args.skip_hitindex and junction_requested:
        junction_note = "off (--skip-hitindex disables junction evidence)"
    else:
        junction_note = _on_off(junction_effective)
    log_message(
        "[INFO]",
        (
            "Qual parameters: "
            f"samples={sample_n}, replicates={replicate_n}, "
            f"reuse context={_on_off(getattr(args, 'reuse', False))}, "
            f"skip HITindex={_on_off(args.skip_hitindex)}, "
            f"HITindex reuse={_on_off(bool(args.hitindex_dir))}, "
            f"AFE/ALE={_on_off(args.calculate_afe_ale)}, "
            f"junction evidence={junction_note}, "
            f"junction_degrade_min_reads={float(getattr(args, 'te_boundary_min_side_reads', 2.0))}, "
            f"te_overlap_min_bp={int(args.te_overlap_min_bp)}, "
            f"te_overlap_min_frac={float(args.te_overlap_min_frac)}, "
            f"splice_site_flank_bp={int(args.splice_site_flank_bp)}"
        ),
        color="info",
    )


def _write_qual_run_config(args, sample_list, replicate_tasks, result_files):
    """Write qual run context to global and debug module configs."""
    replicates = {}
    for sample, _bam, replicate in replicate_tasks:
        replicates.setdefault(str(sample), []).append(str(replicate))
    module_payload = {
        "out_dir": os.path.abspath(args.out_dir),
        "prep": os.path.abspath(args.prep),
        "project": args.project,
        "strand": args.strand,
        "readtype": args.readtype,
        "samples": [str(sample) for sample in sample_list],
        "replicates": {sample: sorted(values) for sample, values in sorted(replicates.items())},
        "calculate_afe_ale": bool(args.calculate_afe_ale),
        "te_overlap_min_bp": int(args.te_overlap_min_bp),
        "te_overlap_min_frac": float(args.te_overlap_min_frac),
        "splice_site_flank_bp": int(args.splice_site_flank_bp),
        "junction": float(getattr(args, "junction", 2.0)),
        "junction_degrade_min_reads": float(getattr(args, "te_boundary_min_side_reads", 2.0)),
        "te_overlap_junction_evidence_requested": bool(
            getattr(args, "_te_overlap_junction_evidence_requested", args.te_overlap_junction_evidence)
        ),
        "te_overlap_junction_evidence": bool(args.te_overlap_junction_evidence),
        "te_overlap_junction_evidence_effective": bool(args.te_overlap_junction_evidence and not args.skip_hitindex),
        "skip_hitindex": bool(args.skip_hitindex),
        "seed": args.seed if args.seed is None else int(args.seed),
        "result_files": result_files,
        "detail": bool(getattr(args, "detail", False)),
        "debug": bool(getattr(args, "debug", False)),
    }
    update_global_config(args.out_dir, "qual", module_payload)
    if args.debug:
        log_dir = os.path.join(args.out_dir, "logs")
        os.makedirs(log_dir, exist_ok=True)
        debug_payload = {"module": "qual", "schema_version": 1, "run_mode": "debug", **module_payload}
        with open(os.path.join(log_dir, "qual_config.json"), "w", encoding="utf-8") as handle:
            json.dump(debug_payload, handle, indent=2, sort_keys=True)


def _adapt_hitindex_reuse_options(args, replicate_tasks):
    """Adjust effective qual options based on reusable HITindex directory contents."""
    args._te_overlap_junction_evidence_requested = bool(args.te_overlap_junction_evidence)
    if not args.hitindex_dir or args.skip_hitindex:
        return
    replicates = [replicate for _sample, _bam, replicate in replicate_tasks]
    if args.te_overlap_junction_evidence:
        junction_ok, junction_detail = inspect_junction_support_reuse(args.hitindex_dir, args.project)
        if junction_ok:
            log_message(
                "[INFO]",
                f"Junction evidence reuse: support table found in --hitindex-dir ({junction_detail}).",
                color="info",
            )
        else:
            log_message(
                "[WARNING]",
                (
                    "Junction evidence disabled for --hitindex-dir reuse because the support table is unavailable "
                    f"or invalid ({junction_detail}). Continuing as if --ignore-junction was used."
                ),
                color="warning",
            )
            args.te_overlap_junction_evidence = False

    afe_status, afe_detail = inspect_afe_ale_reuse(args.hitindex_dir, replicates)
    if afe_status == "complete":
        if not args.calculate_afe_ale:
            log_message(
                "[INFO]",
                "Reusable AFE/ALE HITindex outputs detected; enabling --calculate-afe-ale for this qual run.",
                color="info",
            )
        args.calculate_afe_ale = True
    elif args.calculate_afe_ale or afe_status == "invalid":
        log_message(
            "[WARNING]",
            (
                "AFE/ALE outputs from --hitindex-dir are unavailable or invalid "
                f"({afe_detail}). Continuing without --calculate-afe-ale."
            ),
            color="warning",
        )
        args.calculate_afe_ale = False


def identify_func(args):
    """Identify TE-associated exons and publish mode1 outputs."""
    validate_inputs(args)
    validate_parameters(args)
    validate_runtime(args)
    from util.qual.build_metaexons import metaexon_bed
    from util.qual.write_outputs import (
        _build_junction_summary,
        add_position_summary_to_teexon,
        enrich_afe_ale_outputs,
        export_project_exon_summary,
        write_project_exon_detail,
        publish_afe_ale_outputs,
        write_transcript_exon_te_outputs,
    )
    from util.quant.adapt_qual_outputs import (
        _load_hitindex_position_tables,
        _load_te_exon_event_tx_map,
        _load_transcript_structure_position_table,
    )
    from util.qual.annotate_te_events import _collect_te_overlap_event_annotations
    from util.qual.load_transcript_exons import load_consensus_transcript_exon_rows
    from TExTra.src.mode1.step2_classify import classify_func

    os.makedirs(args.out_dir, exist_ok=True)
    set_log_file(os.path.join(args.out_dir, "logs", "qual.log"))
    extend_log_path = os.path.join(args.out_dir, "logs", "qual_extend.log") if args.debug else None
    if extend_log_path:
        with open(extend_log_path, "w", encoding="utf-8") as handle:
            handle.write("qual debug command/output log\n")

    sample_list = parse_sample_csv(args.samples)
    validate_sample_list(sample_list)
    alignment_dir = os.path.join(args.prep, ALIGNMENT_DIR)
    bamfiles_dict, replicates_dict = collect_bam_and_replicates(alignment_dir, sample_list)
    if not args.skip_hitindex:
        validate_sample_bams(sample_list, bamfiles_dict, replicates_dict)

    exon_bed = os.path.join(args.prep, CONVERT_DIR, 'gene_anno.bed')
    log_message("[INFO]", "Step 1/3: Transcript-exon annotation", bold=True, color="step")

    final_classify_dir = os.path.join(args.out_dir, CLASSIFY_DIR)
    work_out_dir = tempfile.mkdtemp(prefix=f".{CLASSIFY_DIR}.", dir=args.out_dir)
    run_args = copy.copy(args)
    run_args.out_dir = work_out_dir
    classify_dir = os.path.join(run_args.out_dir, CLASSIFY_DIR)
    completed = False
    try:
        if os.path.exists(final_classify_dir):
            log_message(
                "[INFO]",
                f"Staging new classify directory; existing {final_classify_dir} will be replaced after success.",
                color="info",
                console=bool(args.debug),
            )
        else:
            log_message(
                "[INFO]",
                f"Staging classify directory for {final_classify_dir}",
                color="info",
                console=bool(args.debug),
            )
        annotation_dir = os.path.join(classify_dir, "annotation")
        os.makedirs(annotation_dir, exist_ok=True)
        final_annotation_dir = os.path.join(final_classify_dir, "annotation")

        transcript_gtf_path = resolve_consensus_gtf(args.prep)
        assignment_tsv_path = os.path.join(args.prep, ASSEMBLY_DIR, TRANSCRIPT_GENE_ASSIGNMENT_TSV)
        te_bed_path = os.path.join(args.prep, CONVERT_DIR, "TE_anno.bed")
        transcript_exon_rows = load_consensus_transcript_exon_rows(
            transcript_gtf_path,
            multiexon_only=True,
            assignment_tsv_path=assignment_tsv_path,
        )
        transcript_exon_rows_all = load_consensus_transcript_exon_rows(
            transcript_gtf_path,
            multiexon_only=False,
            assignment_tsv_path=assignment_tsv_path,
        )
        te_event_anno_map, _te_overlap_keys = _collect_te_overlap_event_annotations(
            transcript_exon_rows,
            te_bed_path=te_bed_path,
            min_overlap_bp=int(args.te_overlap_min_bp),
            min_overlap_frac=float(args.te_overlap_min_frac),
            boundary_flank_bp=int(args.splice_site_flank_bp),
        )
        annotation_df, annotation_path, teexon_path = write_transcript_exon_te_outputs(
            annotation_dir=annotation_dir,
            classify_dir=classify_dir,
            project_prefix=args.project,
            transcript_exon_rows=transcript_exon_rows,
            te_event_anno_map=te_event_anno_map,
            include_junction_columns=False,
            output_detail=bool(args.detail or args.debug),
            output_metaexon=not args.skip_hitindex,
        )
        log_message(
            "[INFO]",
            (
                "Transcript-exon annotation completed: "
                f"exons={len(annotation_df)}, TE-overlap candidates={_count_te_overlap_rows(annotation_df)}"
            ),
            color="info",
        )

        classify_tmp_dir = os.path.join(classify_dir, "tmp")
        input_buffer_bed = None
        mapping_path = None
        if not args.skip_hitindex:
            validate_gene_bed_matches_consensus(exon_bed, transcript_exon_rows_all)
            log_message("[INFO]", "Step 2/3: HITindex classification", bold=True, color="step")
            log_message("[INFO]", "Preparing HITindex exon model inputs.", color="info")
            os.makedirs(classify_tmp_dir, exist_ok=True)
            input_buffer_bed = metaexon_bed(exon_bed, annotation_dir, run_args, tmp_dir=classify_tmp_dir)
            if args.debug:
                log_message(
                    "[INFO]",
                    f"Generated metaexon BED is stored at {final_annotation_dir}.",
                    color="info",
                )
            else:
                log_message("[INFO]", "Generated metaexon BED.", color="info", console=False)

            mapping_path = os.path.join(annotation_dir, f"{args.project}.metaexon_exon_transcript.tsv")
            mapping_path, meta_n, mapped_n = export_metaexon_exon_transcript_map(
                metaexon_bed_path=os.path.join(annotation_dir, "metaexon.bed"),
                transcript_gtf_path=transcript_gtf_path,
                out_tsv_path=mapping_path,
            )
            annotation_df, annotation_path, teexon_path = write_transcript_exon_te_outputs(
                annotation_dir=annotation_dir,
                classify_dir=classify_dir,
                project_prefix=args.project,
                transcript_exon_rows=transcript_exon_rows,
                te_event_anno_map=te_event_anno_map,
                include_junction_columns=False,
                mapping_tsv_path=mapping_path,
                output_detail=bool(args.detail or args.debug),
                output_metaexon=True,
            )
            if args.debug:
                log_message(
                    "[INFO]",
                    (
                        "Metaexon-exon-transcript map exported: "
                        f"{os.path.join(final_annotation_dir, os.path.basename(mapping_path))} "
                        f"(metaexons={meta_n}, mapped_rows={mapped_n})"
                    ),
                    color="info",
                )
            else:
                log_message(
                    "[INFO]",
                    f"Metaexon-exon-transcript mapping completed: metaexons={meta_n}, mapped_rows={mapped_n}",
                    color="info",
                )
        replicate_tasks = []
        for sample in sample_list:
            for bam, replicate in zip(bamfiles_dict.get(sample, []), replicates_dict.get(sample, [])):
                replicate_tasks.append((sample, bam, replicate))
        _adapt_hitindex_reuse_options(args, replicate_tasks)
        _log_qual_parameters(args, sample_n=len(sample_list), replicate_n=len(replicate_tasks))
        if args.skip_hitindex and args.te_overlap_junction_evidence:
            log_message(
                "[INFO]",
                "Junction evidence disabled because --skip-hitindex was used.",
                color="info",
            )
        if args.hitindex_dir and not args.skip_hitindex:
            copied_n = seed_hitindex_outputs(
                source_hitindex_dir=args.hitindex_dir,
                staged_classify_dir=classify_dir,
                project=args.project,
                replicates=[replicate for _sample, _bam, replicate in replicate_tasks],
                calculate_afe_ale=args.calculate_afe_ale,
                te_overlap_junction_evidence=args.te_overlap_junction_evidence,
            )
            log_message(
                "[INFO]",
                f"HITindex reuse: source={args.hitindex_dir}, reusable_replicates={copied_n}",
                color="info",
            )

        max_jobs = int(args.njobs or args.threads)
        parallel_task_n = len(replicate_tasks) if not args.skip_hitindex else len(sample_list)
        worker_jobs = min(max_jobs, max(1, parallel_task_n))
        per_worker_threads = max(1, args.threads // worker_jobs)
        worker_log_dir = os.path.join(work_out_dir, "logs", "hitindex_workers")
        theano_cache_root = tempfile.gettempdir()

        def _process_replicate(sample, bam, replicate):
            import atexit
            import sys

            from util.common import write_logs as worker_write_logs

            local_args = copy.copy(run_args)
            local_args.threads = per_worker_threads
            worker_log_path = os.path.join(worker_log_dir, f"{_safe_log_name(replicate)}.log")
            replicate_started = time.monotonic()

            old_theano_flags = os.environ.get("THEANO_FLAGS")
            old_worker_log_file = getattr(worker_write_logs, "_LOG_FILE", None)
            unique_theano_dir = None
            if "theano" not in sys.modules:
                unique_theano_dir = tempfile.mkdtemp(prefix=f"TExTra_theano_{replicate}_", dir=theano_cache_root)
                old_flag_parts = [
                    part
                    for part in str(old_theano_flags or "").split(",")
                    if part.strip() and not part.strip().startswith("base_compiledir=")
                ]
                worker_theano_flags = ",".join(old_flag_parts + [f"base_compiledir={unique_theano_dir}"])
                os.environ["THEANO_FLAGS"] = worker_theano_flags
                # Register before Theano is imported; atexit runs Theano cleanup first, then removes this cache.
                atexit.register(lambda path=unique_theano_dir: shutil.rmtree(path, ignore_errors=True))
            worker_write_logs._LOG_FILE = None

            try:
                hitindex_dir = os.path.join(classify_dir, "HITindex")
                status = "computed"
                with _redirect_worker_output(worker_log_path, None):
                    if args.hitindex_dir and is_complete_hitindex_replicate(
                        hitindex_dir,
                        replicate,
                        calculate_afe_ale=args.calculate_afe_ale,
                    ):
                        log_message(
                            "[INFO]",
                            f"Skipping completed HITindex replicate from --hitindex-dir: {replicate}",
                            color="info",
                        )
                        status = "reused"
                    else:
                        classify_func(input_buffer_bed, [bam], [replicate], local_args)
                elapsed = time.monotonic() - replicate_started
                loss_text = (
                    "Finished [100%]: Average Loss = reused"
                    if status == "reused"
                    else _parse_average_loss(worker_log_path)
                )
                return replicate, status, loss_text, elapsed, worker_log_path
            finally:
                if old_theano_flags is None:
                    os.environ.pop("THEANO_FLAGS", None)
                else:
                    os.environ["THEANO_FLAGS"] = old_theano_flags
                worker_write_logs._LOG_FILE = old_worker_log_file

        if args.skip_hitindex:
            log_message("[INFO]", "Step 2/3: HITindex classification", bold=True, color="step")
            log_message("[INFO]", "Skipping HITindex as requested by --skip-hitindex.", color="info")
        else:
            log_message(
                "[INFO]",
                (
                    "Running HITindex classification: "
                    f"replicates={len(replicate_tasks)}, jobs={worker_jobs}, threads/job={per_worker_threads}"
                ),
                color="info",
            )
            log_message(
                "[INFO]",
                "HITindex model fitting may take several minutes; progress is reported when each replicate finishes.",
                color="info",
            )
            os.makedirs(worker_log_dir, exist_ok=True)
            hitindex_started = time.monotonic()
            hitindex_results = []
            hitindex_stream = Parallel(n_jobs=worker_jobs, verbose=0, return_as="generator_unordered")(
                delayed(_process_replicate)(sample, bam, replicate)
                for sample, bam, replicate in replicate_tasks
            )
            total_replicates = len(replicate_tasks)
            for done_n, result in enumerate(hitindex_stream, start=1):
                replicate, status, loss_text, elapsed, worker_log_path = result
                hitindex_results.append(result)
                if extend_log_path:
                    _append_extend_log_section(extend_log_path, f"HITINDEX-{replicate}", worker_log_path)
                if status == "reused":
                    message = (
                        f"HITindex replicate reused [{done_n}/{total_replicates}]: "
                        f"{replicate}, elapsed: {_format_elapsed(elapsed)} finished"
                    )
                else:
                    message = (
                        f"HITindex replicate completed [{done_n}/{total_replicates}]: "
                        f"{replicate}, {loss_text}, elapsed: {_format_elapsed(elapsed)} finished"
                    )
                log_message("[INFO]", message, color="info")
            reused_n = sum(1 for _replicate, status, *_rest in hitindex_results if status == "reused")
            computed_n = sum(1 for _replicate, status, *_rest in hitindex_results if status == "computed")
            log_message(
                "[INFO]",
                (
                    "HITindex completed: "
                    f"computed={computed_n}, reused={reused_n}, replicates={total_replicates}, "
                    f"elapsed: {_format_elapsed(time.monotonic() - hitindex_started)} finished"
                ),
                color="info",
            )
            if args.hitindex_dir and args.calculate_afe_ale:
                written_afe_ale = rebuild_combined_afe_ale_outputs_from_hitindex(
                    classify_dir,
                    args.project,
                    [replicate for _sample, _bam, replicate in replicate_tasks],
                    args,
                )
                log_message(
                    "[INFO]",
                    (
                        "Rebuilt combined AFE/ALE outputs from reusable HITindex PSI files: "
                        + ", ".join(f"{key}={path}" for key, path in sorted(written_afe_ale.items()))
                    ),
                    color="info",
                )
            if args.te_overlap_junction_evidence:
                finalize_te_overlap_support_table(classify_dir, args.project)
                raw_te_candidate_n = int(
                    pd.to_numeric(
                        annotation_df.get("te_overlap_pass_raw", pd.Series(dtype=int)),
                        errors="coerce",
                    ).fillna(0).astype(int).sum()
                )
                support_records = load_te_overlap_support_records(
                    classify_dir,
                    args.project,
                    require_records=raw_te_candidate_n > 0,
                )
                junction_summary = _build_junction_summary(support_records, te_event_anno_map)
                annotation_df, annotation_path, teexon_path = write_transcript_exon_te_outputs(
                    annotation_dir=annotation_dir,
                    classify_dir=classify_dir,
                    project_prefix=args.project,
                    transcript_exon_rows=transcript_exon_rows,
                    te_event_anno_map=te_event_anno_map,
                    include_junction_columns=True,
                    junction_summary=junction_summary,
                    junction_min_side_reads=float(getattr(args, "te_boundary_min_side_reads", 2.0)),
                    mapping_tsv_path=mapping_path,
                    output_detail=bool(args.detail or args.debug),
                    output_metaexon=True,
                )
                log_message(
                    "[INFO]",
                    (
                        "Junction evidence applied: "
                        f"TE-overlap retained={_count_te_overlap_rows(annotation_df)}"
                    ),
                    color="info",
                )
            export_project_exon_summary(classify_dir, args.project, annotation_df)
            if args.calculate_afe_ale:
                enrich_afe_ale_outputs(classify_dir, args.project, annotation_df)

        exon_events, _event_to_txs, _event_to_gene = _load_te_exon_event_tx_map(run_args.out_dir, args.project)
        structure_summary_df = _load_transcript_structure_position_table(run_args.out_dir, args.project)
        position_event_ids = set(exon_events)
        if os.path.isfile(teexon_path):
            try:
                teexon_for_position = pd.read_csv(teexon_path, sep="\t", usecols=lambda col: col == "metaexon_id")
                if "metaexon_id" in teexon_for_position.columns:
                    position_event_ids.update(
                        value
                        for value in teexon_for_position["metaexon_id"].fillna("").astype(str).str.strip()
                        if value
                    )
            except (OSError, ValueError, pd.errors.EmptyDataError, pd.errors.ParserError):
                pass
        position_summary_df, _position_detail_df = _load_hitindex_position_tables(
            run_args.out_dir,
            args.project,
            sorted(position_event_ids),
            structure_summary_df=structure_summary_df,
        )
        if args.detail or args.debug:
            write_project_exon_detail(
                classify_dir,
                args.project,
                annotation_df,
                position_summary_df,
                include_junction_detail=bool(args.te_overlap_junction_evidence),
            )
        add_position_summary_to_teexon(
            teexon_path,
            position_summary_df,
            keep_detail=bool(args.detail or args.debug),
            include_junction_detail=bool(args.te_overlap_junction_evidence),
        )
        if args.calculate_afe_ale:
            publish_afe_ale_outputs(
                classify_dir,
                args.project,
                keep_detail=bool(args.detail or args.debug),
                include_junction_detail=bool(args.te_overlap_junction_evidence),
            )
            for lock_name in [f"{args.project}.AFEPSI.lock", f"{args.project}.ALEPSI.lock"]:
                lock_path = os.path.join(classify_dir, lock_name)
                if os.path.exists(lock_path):
                    os.remove(lock_path)

        _log_detail_qual_summary(args, annotation_df, teexon_path, classify_dir)

        log_message("[INFO]", "Step 3/3: Publish qual outputs", bold=True, color="step")
        if args.debug:
            prune_debug_hitindex_outputs(
                classify_dir,
                args.project,
                calculate_afe_ale=bool(args.calculate_afe_ale),
                te_overlap_junction_evidence=bool(args.te_overlap_junction_evidence),
            )
            if os.path.isdir(annotation_dir):
                shutil.rmtree(annotation_dir, ignore_errors=True)
            if os.path.isdir(classify_tmp_dir):
                shutil.rmtree(classify_tmp_dir, ignore_errors=True)
        else:
            hitindex_dir = os.path.join(classify_dir, "HITindex")
            if os.path.isdir(hitindex_dir):
                shutil.rmtree(hitindex_dir, ignore_errors=True)
            if not args.detail:
                exon_detail_path = os.path.join(classify_dir, f"{args.project}.exon.detail.tsv")
                if os.path.isfile(exon_detail_path):
                    os.remove(exon_detail_path)
            if os.path.isdir(annotation_dir):
                shutil.rmtree(annotation_dir, ignore_errors=True)
            if os.path.isdir(classify_tmp_dir):
                shutil.rmtree(classify_tmp_dir, ignore_errors=True)

        result_files = {
            "transcript_exon_te_annotation": os.path.join(CLASSIFY_DIR, os.path.basename(annotation_path)),
            "te_overlap_exon": os.path.join(CLASSIFY_DIR, os.path.basename(teexon_path)),
        }
        if args.calculate_afe_ale:
            result_files["te_overlap_afepsi"] = os.path.join(CLASSIFY_DIR, f"{args.project}.TE_overlap.AFEPSI.tsv")
            result_files["te_overlap_alepsi"] = os.path.join(CLASSIFY_DIR, f"{args.project}.TE_overlap.ALEPSI.tsv")
            if args.detail or args.debug:
                result_files["te_overlap_afepsi_detail"] = os.path.join(
                    CLASSIFY_DIR,
                    f"{args.project}.TE_overlap.AFEPSI.detail.tsv",
                )
                result_files["te_overlap_alepsi_detail"] = os.path.join(
                    CLASSIFY_DIR,
                    f"{args.project}.TE_overlap.ALEPSI.detail.tsv",
                )
                result_files["afepsi_detail"] = os.path.join(CLASSIFY_DIR, f"{args.project}.AFEPSI.detail.tsv")
                result_files["alepsi_detail"] = os.path.join(CLASSIFY_DIR, f"{args.project}.ALEPSI.detail.tsv")
        if args.detail or args.debug:
            result_files["te_overlap_exon_detail"] = os.path.join(
                CLASSIFY_DIR,
                f"{args.project}.TE_overlap.exon.detail.tsv",
            )
            result_files["exon_detail"] = os.path.join(CLASSIFY_DIR, f"{args.project}.exon.detail.tsv")
        _write_qual_run_config(args, sample_list, replicate_tasks, result_files)

        publish_classify_dir(classify_dir, final_classify_dir)
        completed = True
        log_message("[SUCCESS]", f"Qual output directory: {final_classify_dir}", color="success")
        result_labels = [
            ("transcript_exon_te_annotation", "Transcript-exon annotation"),
            ("te_overlap_exon", "TE-overlap exon table"),
            ("te_overlap_afepsi", "TE-overlap AFEPSI"),
            ("te_overlap_alepsi", "TE-overlap ALEPSI"),
        ]
        if args.detail or args.debug:
            result_labels.extend(
                [
                    ("exon_detail", "Detail exon table"),
                    ("te_overlap_exon_detail", "Detail TE-overlap exon table"),
                    ("afepsi_detail", "Detail AFEPSI"),
                    ("alepsi_detail", "Detail ALEPSI"),
                    ("te_overlap_afepsi_detail", "Detail TE-overlap AFEPSI"),
                    ("te_overlap_alepsi_detail", "Detail TE-overlap ALEPSI"),
                ]
            )
        for key, label in result_labels:
            rel_path = result_files.get(key)
            if rel_path:
                log_message(
                    "[SUCCESS]",
                    f"{label}: {os.path.abspath(os.path.join(args.out_dir, rel_path))}",
                    color="success",
                )
        if args.debug and extend_log_path:
            log_message("[INFO]", f"Debug command/output log: {extend_log_path}", color="info")
    finally:
        if os.path.isdir(work_out_dir):
            if completed or not args.debug:
                shutil.rmtree(work_out_dir, ignore_errors=True)
            else:
                log_message("[WARNING]", f"Keeping failed qual staging directory: {work_out_dir}", color="warning")

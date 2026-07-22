"""Run HITindex metrics and combined PSI exports for qual classification."""

import copy
import os

import pandas as pd
from filelock import FileLock
from pandas.errors import EmptyDataError, ParserError

from util.qual.define_hitindex_constants import (
    HIT_IDENTITY_PARAMETERS,
    HITINDEX_MIN_JUNCTION_READS,
    HITINDEX_SPLIT_OVERLAP_BP,
    TE_JUNCTION_INTERNAL_GENE_PREFIX,
)
from util.qual.build_metaexons import (
    _build_te_overlap_unique_exon_rows,
    _build_union_bed,
    _build_unique_exon_bed,
    _copy_bed_with_gene_prefix,
    _flatten_unique_exon_support,
)
from util.qual.annotate_te_events import (
    _build_exon_event_key,
    _build_metaexon_te_annotation_map,
    _collect_exon_te_annotations,
    _collect_te_overlap_event_annotations,
)
from util.qual.load_transcript_exons import (
    _load_mapping_transcript_exon_rows,
)
from util.qual.calculate_hitindex import (
    calculateMetric,
    calculatePSI,
    call_terminal,
    edge_flagging,
    getExonReads,
    getExons,
    getJuncBed,
    readMetrics,
    running_genmodel,
    significance,
    writeMetrics,
)
from util.qual.extract_junctions import extract_junction
from util.common.define_layout import CONVERT_DIR
from util.common.write_logs import log_message


def HITindex_pipeline(input_buffer_bed, bamfiles, samples, output_hitindex_dir, args, combined_out_dir=None):
    """Run per-sample HITindex metrics and export combined PSI tables."""
    tmp_dir = os.path.join(output_hitindex_dir, 'tmp')
    os.makedirs(tmp_dir, exist_ok=True)
    if combined_out_dir is None:
        combined_out_dir = output_hitindex_dir
    os.makedirs(combined_out_dir, exist_ok=True)

    exon_outnames = []
    combined_afe = []
    combined_ale = []
    project_prefix = str(getattr(args, "project", "project"))
    te_bed_path = os.path.join(args.prep, CONVERT_DIR, 'TE_anno.bed')
    min_overlap_bp = int(getattr(args, "te_overlap_min_bp", 10))
    min_overlap_frac = float(getattr(args, "te_overlap_min_frac", 0.1))
    site_flank_bp = int(getattr(args, "splice_site_flank_bp", 10))
    beddict = getExons(input_buffer_bed)
    paramdict = dict(HIT_IDENTITY_PARAMETERS)
    per_exon_junction_records = []
    unique_exon_bed = None
    unique_exon_metric_bed = None
    unique_exon_beddict_template = None
    intersect_bed = input_buffer_bed
    unique_event_meta = {}
    te_event_anno_map = {}
    annotation_dir = None
    transcript_exon_rows = []
    te_overlap_junction_enabled = bool(getattr(args, "te_overlap_junction_evidence", False))
    if te_overlap_junction_enabled:
        annotation_dir = os.path.join(os.path.dirname(output_hitindex_dir), "annotation")
        mapping_tsv = os.path.join(annotation_dir, f"{project_prefix}.metaexon_exon_transcript.tsv")
        transcript_exon_rows = _load_mapping_transcript_exon_rows(mapping_tsv)
        te_event_anno_map, te_overlap_keys = _collect_te_overlap_event_annotations(
            transcript_exon_rows,
            te_bed_path=te_bed_path,
            min_overlap_bp=min_overlap_bp,
            min_overlap_frac=min_overlap_frac,
            boundary_flank_bp=site_flank_bp,
        )
        unique_te_overlap_rows, unique_event_meta = _build_te_overlap_unique_exon_rows(
            transcript_exon_rows,
            te_overlap_keys,
        )
        unique_exon_bed = _build_unique_exon_bed(
            unique_te_overlap_rows,
            os.path.join(tmp_dir, f"{project_prefix}.te_overlap_unique_exon_junction_check.bed"),
        )
        if unique_exon_bed is None:
            log_message(
                "[WARNING]",
                "TE-overlap junction evidence is enabled, but no TE-overlap unique-exon rows were found. "
                "Skip TE-overlap exon junction support export.",
                color="warning",
            )
        else:
            unique_exon_metric_bed = _copy_bed_with_gene_prefix(
                unique_exon_bed,
                os.path.join(tmp_dir, f"{project_prefix}.te_overlap_unique_exon_junction_check.prefixed.bed"),
                TE_JUNCTION_INTERNAL_GENE_PREFIX,
            )
            intersect_bed = _build_union_bed(
                input_buffer_bed,
                unique_exon_metric_bed,
                os.path.join(tmp_dir, f"{project_prefix}.hitindex_and_te_junction_union.bed"),
            )
            unique_exon_beddict_template = getExons(unique_exon_metric_bed)

    # Loop through BAM files and corresponding sample names.
    for bam, sample in zip(bamfiles, samples):
        # Calculating HITindex metrics.
        juncbam = extract_junction(bam, sample, tmp_dir, args)
        log_message("[INFO]", f"Processing HITindex metrics for sample: {sample}", color="info")
        getJuncBed(intersect_bed, juncbam, args.readtype, args.strand)
        startdict, enddict = getExonReads(juncbam, HITINDEX_SPLIT_OVERLAP_BP)
        HITdict = calculateMetric(beddict, startdict, enddict, HITINDEX_MIN_JUNCTION_READS)
        HITdict, param = edge_flagging(HITdict)
        exon_outname = f'{output_hitindex_dir}/{sample}.exon'
        writeMetrics(HITdict, param, exon_outname)
        seed = getattr(args, 'seed', None)
        running_genmodel(exon_outname, getattr(args, 'genmodel_iters', 100000), seed=seed)

        # Exon classification.
        HITcombo = readMetrics(exon_outname)
        HITcombo = significance(
            HITcombo,
            getattr(args, 'bootstrap_n', 1000),
            float(paramdict['HIThybrid']),
            exon_outname,
            seed=seed,
        )
        call_terminal(HITcombo, paramdict, exon_outname)
        if getattr(args, 'calculate_afe_ale', False):
            psi_outname = f'{output_hitindex_dir}/{sample}'
            calculatePSI(HITcombo, psi_outname)
            afe_path = psi_outname + '.AFEPSI'
            ale_path = psi_outname + '.ALEPSI'
            if os.path.isfile(afe_path):
                afe_df = pd.read_csv(afe_path, sep='\t')
                if not afe_df.empty:
                    afe_df['Sample'] = sample
                    combined_afe.append(afe_df)
            if os.path.isfile(ale_path):
                ale_df = pd.read_csv(ale_path, sep='\t')
                if not ale_df.empty:
                    ale_df['Sample'] = sample
                    combined_ale.append(ale_df)

        if unique_exon_metric_bed:
            tx_beddict = copy.deepcopy(unique_exon_beddict_template)
            tx_hitdict = calculateMetric(tx_beddict, startdict, enddict, HITINDEX_MIN_JUNCTION_READS)
            per_exon_junction_records.extend(
                _flatten_unique_exon_support(
                    tx_hitdict,
                    sample,
                    overlap_bp=HITINDEX_SPLIT_OVERLAP_BP,
                    read_threshold=HITINDEX_MIN_JUNCTION_READS,
                    te_event_anno_map=te_event_anno_map,
                    unique_event_meta=unique_event_meta,
                )
            )

        exon_outnames.append(exon_outname)

    metaexon_te_anno_map = {}
    if te_overlap_junction_enabled and te_event_anno_map:
        metaexon_te_anno_map = _build_metaexon_te_annotation_map(
            transcript_exon_rows,
            te_event_anno_map,
        )

    def _write_wide(df_list, value_col, fname):
        if not df_list:
            return
        all_df = pd.concat(df_list, ignore_index=True).copy()
        if all_df.empty:
            return
        metric_defaults = {
            "label": "no_overlap",
            "te_overlap_n": 0,
            "te_overlap_bp_max": 0,
            "te_overlap_frac_max": 0.0,
            "te_boundary_hit_any": 0,
            "te_overlap_pass_any": 0,
            "te_splice_site_repeat_TE": "",
            "te_other_overlap_TE": "",
        }
        all_df["__event_key"] = all_df.apply(lambda r: _build_exon_event_key(r["exon"], r["strand"]), axis=1)
        anno_map = _collect_exon_te_annotations(
            all_df,
            te_bed_path,
            min_overlap_bp=min_overlap_bp,
            min_overlap_frac=min_overlap_frac,
            boundary_flank_bp=site_flank_bp,
        )

        for col, default_v in metric_defaults.items():
            all_df[col] = all_df["__event_key"].map(
                lambda k: metaexon_te_anno_map.get(k, anno_map.get(k, {})).get(col, default_v)
            )

        index_cols = [
            'gene', 'exon', 'strand', 'label',
            'te_overlap_n', 'te_overlap_bp_max', 'te_overlap_frac_max',
            'te_boundary_hit_any', 'te_overlap_pass_any',
            'te_splice_site_repeat_TE', 'te_other_overlap_TE',
        ]
        wide = all_df.pivot_table(index=index_cols, columns='Sample', values=value_col, aggfunc='first')
        # Ensure all declared samples appear as columns, filling missing values with NA.
        wide = wide.reindex(columns=samples)
        out_path = os.path.join(combined_out_dir, fname)
        lock = FileLock(out_path + ".lock")
        with lock:
            if os.path.exists(out_path):  # merge with existing to avoid overwrite between groups
                prev = pd.read_csv(out_path, sep='\t')
                prev_samples = [c for c in prev.columns if c not in index_cols]
                merged = prev.set_index(index_cols).reindex(columns=prev_samples)
                # union sample columns
                sample_union = list(dict.fromkeys(prev_samples + list(wide.columns)))  # preserve order
                merged = merged.reindex(columns=sample_union)
                wide = wide.reindex(columns=sample_union)
                combined = pd.concat([merged, wide]).groupby(level=list(range(len(index_cols)))).first()
                combined_out = combined.reset_index()
                combined_out.to_csv(out_path, sep='\t', index=False)
            else:
                wide_out = wide.reset_index()
                wide_out.to_csv(out_path, sep='\t', index=False)

    _write_wide(combined_afe, 'AFEPSI', f"{project_prefix}.AFEPSI")
    _write_wide(combined_ale, 'ALEPSI', f"{project_prefix}.ALEPSI")

    if unique_exon_bed:
        support_cols = [
            "sample",
            "gene_id",
            "transcript_id",
            "exon",
            "chrom",
            "start",
            "end",
            "strand",
            "nleft",
            "nright",
            "total_junction_reads",
            "junction_supported",
            "junction_support_rule",
            "te_overlap_label",
            "te_overlap_n",
            "te_overlap_bp_max",
            "te_overlap_frac_max",
            "te_boundary_hit_any",
            "te_overlap_pass_any",
            "te_splice_site_repeat_TE",
            "te_other_overlap_TE",
        ]
        support_path = os.path.join(output_hitindex_dir, f"{project_prefix}.te_overlap_transcript_exon_junction_support.tsv")
        if per_exon_junction_records:
            support_df = pd.DataFrame(per_exon_junction_records, columns=support_cols)
            support_df = support_df.sort_values(
                by=["sample", "gene_id", "transcript_id", "chrom", "start", "end"],
                ascending=[True, True, True, True, True, True],
            )
        else:
            support_df = pd.DataFrame(columns=support_cols)
        lock = FileLock(support_path + ".lock")
        with lock:
            if os.path.exists(support_path):
                try:
                    prev_df = pd.read_csv(support_path, sep="\t")
                except (OSError, EmptyDataError, ParserError, UnicodeDecodeError):
                    prev_df = pd.DataFrame(columns=support_cols)
            else:
                prev_df = pd.DataFrame(columns=support_cols)

            for col in support_cols:
                if col not in prev_df.columns:
                    prev_df[col] = ""
                if col not in support_df.columns:
                    support_df[col] = ""
            prev_df = prev_df[support_cols]
            support_df = support_df[support_cols]

            merged_df = pd.concat([prev_df, support_df], ignore_index=True)
            if not merged_df.empty:
                dedup_cols = ["sample", "exon", "strand", "chrom", "start", "end"]
                merged_df = merged_df.drop_duplicates(subset=dedup_cols, keep="last")
                merged_df = merged_df.sort_values(
                    by=["sample", "gene_id", "transcript_id", "chrom", "start", "end"],
                    ascending=[True, True, True, True, True, True],
                )
            merged_df.to_csv(support_path, sep="\t", index=False)

    return exon_outnames

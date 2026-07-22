"""Write qual annotation outputs and summary enrichment tables."""

import collections
import glob
import os

import pandas as pd
from pandas.errors import EmptyDataError, ParserError

from util.qual.define_hitindex_constants import TE_BOUNDARY_MIN_SIDE_READS
from util.qual.annotate_te_events import _build_exon_event_key, _normalize_int
from util.qual.load_transcript_exons import _make_transcript_exon_id


def _infer_sample_id(replicate_label):
    value = str(replicate_label or "").strip()
    if "_rep" in value:
        prefix, suffix = value.rsplit("_rep", 1)
        if prefix and suffix.isdigit():
            return prefix
    return value


def _split_clean_tokens(value):
    tokens = []
    for token in str(value or "").split(","):
        token = token.strip()
        if token and token.lower() not in {"nan", "none", "na"}:
            tokens.append(token)
    return tokens


def _join_clean_tokens(values):
    return ",".join(sorted({str(v).strip() for v in values if str(v).strip()}))


def _join_split_tokens(values):
    tokens = []
    for value in values:
        tokens.extend(_split_clean_tokens(value))
    return _join_clean_tokens(tokens)


def _first_clean_value(values, default=""):
    for value in values:
        text = str(value).strip()
        if text and text.lower() not in {"nan", "none"}:
            return text
    return default


def _max_numeric_value(values, default=0):
    numeric = pd.to_numeric(pd.Series(list(values)), errors="coerce").dropna()
    if numeric.empty:
        return default
    value = numeric.max()
    return int(value) if float(value).is_integer() else float(value)


def _event_id_from_row(row):
    try:
        return (
            f"{row['exon_chrom']}:{int(row['exon_start'])}-"
            f"{int(row['exon_end'])}:{row['exon_strand']}"
        )
    except (KeyError, TypeError, ValueError):
        return ""


def _te_boundary_side_label(left_hit, right_hit, strand):
    if int(left_hit) and int(right_hit):
        return "both"
    if not int(left_hit) and not int(right_hit):
        return "none"
    if strand == "-":
        return "3'" if int(left_hit) else "5'"
    if strand == "+":
        return "5'" if int(left_hit) else "3'"
    return "unknown"


def _collapse_boundary_side(values):
    tokens = {token for value in values for token in _split_clean_tokens(value)}
    if not tokens:
        return "none"
    if "both" in tokens or {"5'", "3'"}.issubset(tokens):
        return "both"
    if "5'" in tokens:
        return "5'"
    if "3'" in tokens:
        return "3'"
    if "unknown" in tokens:
        return "unknown"
    return "none"


def _collapse_te_overlap_label(values):
    tokens = {str(value).strip() for value in values if str(value).strip()}
    return "TE_overlap" if "TE_overlap" in tokens else _first_clean_value(tokens, default="no_overlap")


def _build_transcript_exon_te_annotation_df(
    transcript_exon_rows,
    te_event_anno_map,
    junction_summary=None,
    include_junction_columns=True,
    junction_min_side_reads=TE_BOUNDARY_MIN_SIDE_READS,
):
    cols = [
        "metaexon_id",
        "transcript_exon_id",
        "gene_id",
        "transcript_id",
        "exon_id",
        "exon_number",
        "transcript_exon_count",
        "is_multiexon_transcript",
        "exon_chrom",
        "exon_start",
        "exon_end",
        "exon_strand",
        "gene_strand",
        "reference_gene_strand",
        "transcript_strand",
        "gene_exon_strand_consistent",
        "gene_assignment_strand_consistent",
        "gene_assignment_status",
        "te_gene_strand_relation",
        "te_overlap_label_raw",
        "te_overlap_pass_raw",
        "te_overlap_n",
        "te_overlap_bp_max",
        "te_overlap_frac_max",
        "te_boundary_hit_any",
        "te_left_boundary_hit_any",
        "te_right_boundary_hit_any",
        "te_boundary_side",
        "te_splice_site_repeat_TE",
        "te_other_overlap_TE",
        "te_overlap_label",
        "te_overlap_pass_any",
    ]
    junction_cols = [
        "junction_evidence_enabled",
        "junction_support_sample_n",
        "junction_support_sample_ids",
        "junction_te_side_reads_max",
        "junction_te_side_reads_mean_supported_samples",
        "junction_supported",
        "junction_pass",
        "te_overlap_degraded",
        "te_overlap_degrade_reason",
        "te_overlap_label_final",
        "te_overlap_pass_final",
    ]
    if include_junction_columns:
        cols.extend(junction_cols)

    rows = []
    for row in transcript_exon_rows:
        exon = f'{row["exon_chrom"]}:{int(row["exon_start"])}-{int(row["exon_end"])}'
        strand = str(row["exon_strand"])
        event_key = _build_exon_event_key(exon, strand)
        ann = te_event_anno_map.get(event_key, {})
        raw_pass = int(ann.get("te_overlap_pass_any", 0))
        raw_label = "TE_overlap" if raw_pass else "no_overlap"
        exon_id = f'{row["exon_chrom"]}:{int(row["exon_start"])}-{int(row["exon_end"])}:{strand}'
        left_hit = int(ann.get("te_left_boundary_hit_any", 0))
        right_hit = int(ann.get("te_right_boundary_hit_any", 0))
        boundary_side = _te_boundary_side_label(left_hit, right_hit, strand)

        rec = {
            "metaexon_id": str(row.get("metaexon_id", "")),
            "transcript_exon_id": str(row.get("transcript_exon_id", _make_transcript_exon_id(row))),
            "gene_id": str(row.get("gene_id", "")),
            "transcript_id": str(row.get("transcript_id", "")),
            "exon_id": exon_id,
            "exon_number": int(row.get("exon_number", 0)),
            "transcript_exon_count": int(row.get("transcript_exon_count", 0)),
            "is_multiexon_transcript": int(row.get("is_multiexon_transcript", 0)),
            "exon_chrom": str(row.get("exon_chrom", "")),
            "exon_start": int(row.get("exon_start", -1)),
            "exon_end": int(row.get("exon_end", -1)),
            "exon_strand": strand,
            "gene_strand": str(row.get("gene_strand", strand)),
            "reference_gene_strand": str(row.get("reference_gene_strand", "")),
            "transcript_strand": str(row.get("transcript_strand", strand)),
            "gene_exon_strand_consistent": int(row.get("gene_exon_strand_consistent", 1)),
            "gene_assignment_strand_consistent": str(row.get("gene_assignment_strand_consistent", "")),
            "gene_assignment_status": str(row.get("gene_assignment_status", "")),
            "te_gene_strand_relation": str(row.get("te_gene_strand_relation", "same")),
            "te_overlap_label_raw": raw_label,
            "te_overlap_pass_raw": raw_pass,
            "te_overlap_n": int(ann.get("te_overlap_n", 0)),
            "te_overlap_bp_max": int(ann.get("te_overlap_bp_max", 0)),
            "te_overlap_frac_max": float(ann.get("te_overlap_frac_max", 0.0)),
            "te_boundary_hit_any": int(ann.get("te_boundary_hit_any", 0)),
            "te_left_boundary_hit_any": left_hit,
            "te_right_boundary_hit_any": right_hit,
            "te_boundary_side": boundary_side,
            "te_splice_site_repeat_TE": str(ann.get("te_splice_site_repeat_TE", "")),
            "te_other_overlap_TE": str(ann.get("te_other_overlap_TE", "")),
            "te_overlap_label": raw_label,
            "te_overlap_pass_any": raw_pass,
        }

        if include_junction_columns:
            summary = (junction_summary or {}).get(event_key, {})
            junction_max = float(summary.get("junction_te_side_reads_max", 0.0))
            junction_pass = int(junction_max >= float(junction_min_side_reads))
            degraded = int(raw_pass == 1 and junction_pass == 0)
            final_pass = int(raw_pass == 1 and degraded == 0)
            rec.update(
                {
                    "junction_evidence_enabled": 1,
                    "te_boundary_side": boundary_side,
                    "junction_support_sample_n": int(summary.get("junction_support_sample_n", 0)),
                    "junction_support_sample_ids": str(summary.get("junction_support_sample_ids", "")),
                    "junction_te_side_reads_max": junction_max,
                    "junction_te_side_reads_mean_supported_samples": float(
                        summary.get("junction_te_side_reads_mean_supported_samples", 0.0)
                    ),
                    "junction_supported": int(summary.get("junction_supported", 0)),
                    "junction_pass": junction_pass,
                    "te_overlap_degraded": degraded,
                    "te_overlap_degrade_reason": (
                        f"te_side_junction_reads_max<{float(junction_min_side_reads):g}" if degraded else ""
                    ),
                    "te_overlap_label_final": "TE_overlap" if final_pass else "no_overlap",
                    "te_overlap_pass_final": final_pass,
                }
            )
            rec["te_overlap_label"] = rec["te_overlap_label_final"]
            rec["te_overlap_pass_any"] = rec["te_overlap_pass_final"]
        rows.append(rec)

    if rows:
        out_df = pd.DataFrame(rows, columns=cols).drop_duplicates()
        out_df = out_df.sort_values(
            by=["gene_id", "transcript_id", "exon_number", "exon_chrom", "exon_start", "exon_end"],
            ascending=[True, True, True, True, True, True],
        )
    else:
        out_df = pd.DataFrame(columns=cols)
    return out_df


def _select_transcript_exon_te_annotation_columns(
    annotation_df,
    *,
    include_detail=False,
    include_metaexon=True,
    include_junction_evidence=False,
):
    base_cols = [
        "exon_id",
        "gene_id",
        "transcript_id",
        "exon_number",
        "transcript_exon_count",
        "transcript_strand",
        "exon_chrom",
        "exon_start",
        "exon_end",
        "exon_strand",
        "te_overlap_label",
        "te_splice_site_repeat_TE",
    ]
    if include_metaexon:
        base_cols.insert(1, "metaexon_id")
    detail_cols = [
        "te_overlap_bp_max",
        "te_overlap_frac_max",
        "te_boundary_hit_any",
        "te_boundary_side",
    ]
    if include_junction_evidence:
        detail_cols.extend(
            [
                "junction_support_sample_n",
                "junction_support_sample_ids",
                "junction_te_side_reads_max",
                "junction_te_side_reads_mean_supported_samples",
            ]
        )
    selected = base_cols + (detail_cols if include_detail else [])
    return annotation_df[[col for col in selected if col in annotation_df.columns]].copy()


def _apply_metaexon_ids_to_annotation(annotation_df, mapping_tsv_path):
    if annotation_df.empty:
        return annotation_df
    if not os.path.isfile(mapping_tsv_path):
        raise RuntimeError(f"Metaexon-transcript mapping table is missing: {mapping_tsv_path}")
    try:
        mapping_df = pd.read_csv(mapping_tsv_path, sep="\t")
    except (OSError, EmptyDataError, ParserError, UnicodeDecodeError) as exc:
        raise RuntimeError(f"Failed to read metaexon-transcript mapping table: {mapping_tsv_path}") from exc
    required = {"metaexon_id", "transcript_id", "exon_chrom", "exon_start", "exon_end", "exon_strand"}
    if not required.issubset(mapping_df.columns):
        missing = sorted(required - set(mapping_df.columns))
        raise RuntimeError(
            f"Invalid metaexon-transcript mapping table: {mapping_tsv_path} missing column(s): {', '.join(missing)}"
        )
    work = mapping_df[list(required)].copy()
    work["exon_start"] = pd.to_numeric(work["exon_start"], errors="coerce")
    work["exon_end"] = pd.to_numeric(work["exon_end"], errors="coerce")
    work = work.dropna(subset=["exon_start", "exon_end"]).copy()
    work["exon_start"] = work["exon_start"].astype(int)
    work["exon_end"] = work["exon_end"].astype(int)
    work = work.drop_duplicates(
        subset=["transcript_id", "exon_chrom", "exon_start", "exon_end", "exon_strand"]
    )
    merged = annotation_df.drop(columns=["metaexon_id"], errors="ignore").merge(
        work,
        on=["transcript_id", "exon_chrom", "exon_start", "exon_end", "exon_strand"],
        how="left",
    )
    merged["metaexon_id"] = merged["metaexon_id"].fillna("").astype(str)
    ordered_cols = list(annotation_df.columns)
    if "metaexon_id" in ordered_cols:
        merged = merged[ordered_cols]
    return merged


def _build_junction_summary(per_exon_junction_records, te_event_anno_map):
    support_by_event = collections.defaultdict(lambda: collections.defaultdict(list))
    for row in per_exon_junction_records or []:
        exon = str(row.get("exon", "")).strip()
        strand = str(row.get("strand", "")).strip()
        if not exon or not strand:
            continue
        event_key = _build_exon_event_key(exon, strand)
        ann = te_event_anno_map.get(event_key, {})
        left_hit = int(ann.get("te_left_boundary_hit_any", 0))
        right_hit = int(ann.get("te_right_boundary_hit_any", 0))
        nleft = _normalize_int(row.get("nleft", 0), default=0)
        nright = _normalize_int(row.get("nright", 0), default=0)
        if left_hit and right_hit:
            te_side_reads = max(nleft, nright)
        elif left_hit:
            te_side_reads = nleft
        elif right_hit:
            te_side_reads = nright
        else:
            te_side_reads = max(nleft, nright)
        sample_id = _infer_sample_id(row.get("sample", ""))
        support_by_event[event_key][sample_id].append(int(te_side_reads))

    summary = {}
    for event_key, by_sample in support_by_event.items():
        sample_values = {}
        for sample_id, values in by_sample.items():
            supported_values = [int(v) for v in values if int(v) >= 1]
            sample_values[sample_id] = (sum(supported_values) / float(len(supported_values))) if supported_values else 0.0
        supported_samples = sorted([sample for sample, value in sample_values.items() if value >= 1])
        supported_values = [sample_values[sample] for sample in supported_samples]
        max_reads = max(supported_values) if supported_values else 0.0
        mean_supported = (sum(supported_values) / float(len(supported_values))) if supported_values else 0.0
        summary[event_key] = {
            "junction_support_sample_n": int(len(supported_samples)),
            "junction_support_sample_ids": ",".join(supported_samples),
            "junction_te_side_reads_max": float(max_reads),
            "junction_te_side_reads_mean_supported_samples": float(mean_supported),
            "junction_supported": int(len(supported_samples) >= 1),
        }
    return summary


def write_transcript_exon_te_outputs(
    annotation_dir,
    classify_dir,
    project_prefix,
    transcript_exon_rows,
    te_event_anno_map,
    include_junction_columns=True,
    junction_summary=None,
    junction_min_side_reads=TE_BOUNDARY_MIN_SIDE_READS,
    mapping_tsv_path=None,
    output_detail=False,
    output_metaexon=True,
):
    """Write transcript-exon TE annotation and TE-exon summary tables."""
    annotation_df = _build_transcript_exon_te_annotation_df(
        transcript_exon_rows,
        te_event_anno_map,
        junction_summary=junction_summary,
        include_junction_columns=include_junction_columns,
        junction_min_side_reads=junction_min_side_reads,
    )
    if mapping_tsv_path:
        annotation_df = _apply_metaexon_ids_to_annotation(annotation_df, mapping_tsv_path)
    os.makedirs(annotation_dir, exist_ok=True)
    os.makedirs(classify_dir, exist_ok=True)
    annotation_path = os.path.join(classify_dir, "transcript_exon_te_annotation.tsv")
    output_annotation_df = _select_transcript_exon_te_annotation_columns(
        annotation_df,
        include_detail=bool(output_detail),
        include_metaexon=bool(output_metaexon),
        include_junction_evidence=bool(include_junction_columns),
    )
    output_annotation_df.to_csv(annotation_path, sep="\t", index=False)

    pass_col = "te_overlap_pass_final" if "te_overlap_pass_final" in annotation_df.columns else "te_overlap_pass_any"
    teexon_df = annotation_df.loc[
        pd.to_numeric(annotation_df.get(pass_col, 0), errors="coerce").fillna(0).astype(int) == 1
    ].copy()
    teexon_path = os.path.join(classify_dir, f"{project_prefix}.TE_overlap.exon.tsv")
    teexon_df.to_csv(teexon_path, sep="\t", index=False)
    return annotation_df, annotation_path, teexon_path


def _aggregate_te_overlap_exon_table(te_df):
    if te_df.empty:
        return te_df
    work = te_df.copy()
    if "exon_id" not in work.columns:
        work["exon_id"] = work.apply(_event_id_from_row, axis=1)
    work["exon_id"] = work["exon_id"].fillna("").astype(str).str.strip()
    work = work.loc[work["exon_id"] != ""].copy()
    if work.empty:
        return work

    rows = []
    for exon_id, group in work.groupby("exon_id", sort=True):
        row = {
            "exon_id": exon_id,
            "gene_id": _join_split_tokens(group.get("gene_id", [])),
            "transcript_id": _join_split_tokens(group.get("transcript_id", [])),
            "transcript_exon_id": _join_split_tokens(group.get("transcript_exon_id", [])),
            "transcript_strand": _join_split_tokens(group.get("transcript_strand", [])),
            "exon_chrom": _first_clean_value(group.get("exon_chrom", [])),
            "exon_start": _max_numeric_value(group.get("exon_start", []), default=0),
            "exon_end": _max_numeric_value(group.get("exon_end", []), default=0),
            "exon_strand": _first_clean_value(group.get("exon_strand", [])),
            "te_overlap_label": _collapse_te_overlap_label(group.get("te_overlap_label", [])),
            "te_splice_site_repeat_TE": _join_split_tokens(group.get("te_splice_site_repeat_TE", [])),
            "te_overlap_bp_max": _max_numeric_value(group.get("te_overlap_bp_max", []), default=0),
            "te_overlap_frac_max": _max_numeric_value(group.get("te_overlap_frac_max", []), default=0.0),
            "te_boundary_side": _collapse_boundary_side(group.get("te_boundary_side", [])),
            "junction_te_side_reads_max": _max_numeric_value(
                group.get("junction_te_side_reads_max", []),
                default=0.0,
            ),
            "ID_position_summary": _first_clean_value(group.get("ID_position_summary", []), default="NA"),
            "ID_position_source": _first_clean_value(group.get("ID_position_source", []), default="NA"),
            "ID_position_roles": _join_split_tokens(group.get("ID_position_roles", [])) or "NA",
            "ID_position_confidence": _first_clean_value(group.get("ID_position_confidence", []), default="NA"),
            "transcript_structure_roles": _join_split_tokens(group.get("transcript_structure_roles", [])) or "NA",
            "position_evidence_relation": _first_clean_value(group.get("position_evidence_relation", []), default="NA"),
            "candidate_TE_event": _first_clean_value(group.get("candidate_TE_event", []), default="NA"),
            "candidate_TE_confidence": _first_clean_value(group.get("candidate_TE_confidence", []), default="NA"),
            "ID_position_evaluable_sample_n": _max_numeric_value(
                group.get("ID_position_evaluable_sample_n", []),
                default=0,
            ),
            "ID_position_evaluable_replicate_n": _max_numeric_value(
                group.get("ID_position_evaluable_replicate_n", []),
                default=0,
            ),
        }
        if "metaexon_id" in group.columns:
            row["metaexon_id"] = _join_split_tokens(group["metaexon_id"])
        rows.append(row)

    out = pd.DataFrame(rows)
    preferred = [
        "exon_id",
        "metaexon_id",
        "gene_id",
        "transcript_id",
        "transcript_exon_id",
        "transcript_strand",
        "exon_chrom",
        "exon_start",
        "exon_end",
        "exon_strand",
        "te_overlap_label",
        "te_splice_site_repeat_TE",
        "te_overlap_bp_max",
        "te_overlap_frac_max",
        "te_boundary_side",
        "ID_position_summary",
        "ID_position_source",
        "ID_position_roles",
        "ID_position_confidence",
        "transcript_structure_roles",
        "position_evidence_relation",
        "candidate_TE_event",
        "candidate_TE_confidence",
        "ID_position_evaluable_sample_n",
        "ID_position_evaluable_replicate_n",
    ]
    ordered = [col for col in preferred if col in out.columns]
    ordered.extend([col for col in out.columns if col not in ordered])
    return out[ordered]


def _select_te_overlap_exon_output_columns(df, *, include_detail=False, include_junction_detail=False):
    if df.empty:
        return df
    out = df.copy()
    rename_map = {}
    if include_detail:
        rename_map = {
            "ID_position_roles": "HITindex_structure_roles",
            "ID_position_evaluable_sample_n": "HITindex_evaluable_sample_n",
            "ID_position_evaluable_replicate_n": "HITindex_evaluable_replicate_n",
        }
        out = out.rename(columns=rename_map)
    default_cols = [
        "exon_id",
        "metaexon_id",
        "gene_id",
        "transcript_id",
        "transcript_exon_id",
        "te_overlap_label",
        "te_splice_site_repeat_TE",
        "te_boundary_side",
        "ID_position_summary",
        "ID_position_confidence",
        "candidate_TE_event",
        "candidate_TE_confidence",
    ]
    detail_cols = [
        "te_overlap_bp_max",
        "te_overlap_frac_max",
        "ID_position_source",
        "transcript_structure_roles",
        "HITindex_structure_roles",
        "HITindex_evaluable_sample_n",
        "HITindex_evaluable_replicate_n",
    ]
    if include_junction_detail:
        detail_cols.insert(2, "junction_te_side_reads_max")
    selected = default_cols + (detail_cols if include_detail else [])
    return out[[col for col in selected if col in out.columns]].copy()


def _merge_position_summary_columns(out, position_summary_df, summary_cols):
    if position_summary_df is not None and not position_summary_df.empty:
        summary = position_summary_df.copy()
        summary["event_id"] = summary["event_id"].astype(str)
        summary = summary.drop_duplicates(subset=["event_id"], keep="first")
        summary["__position_match_id"] = summary["event_id"]
        keep_cols = ["__position_match_id"] + [col for col in summary_cols if col in summary.columns]
        out["__position_match_id"] = out["event_id"]
        if "metaexon_id" in out.columns:
            meta_ids = out["metaexon_id"].fillna("").astype(str).str.strip()
            out.loc[meta_ids != "", "__position_match_id"] = meta_ids.loc[meta_ids != ""]
        out = out.drop(columns=summary_cols, errors="ignore").merge(
            summary[keep_cols],
            on="__position_match_id",
            how="left",
        )
        out = out.drop(columns=["__position_match_id"], errors="ignore")
    for col in summary_cols:
        if col not in out.columns:
            out[col] = 0 if col.endswith("_n") else "NA"
        elif col.endswith("_n"):
            out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0).astype(int)
        else:
            out[col] = out[col].fillna("NA").astype(str)
    return out


def write_project_exon_detail(
    classify_dir,
    project_prefix,
    annotation_df,
    position_summary_df,
    include_junction_detail=False,
):
    """Write all exon-level detail rows using the same schema as TE-overlap exon detail."""
    if annotation_df.empty:
        out_path = os.path.join(classify_dir, f"{project_prefix}.exon.detail.tsv")
        pd.DataFrame().to_csv(out_path, sep="\t", index=False)
        return out_path
    out = annotation_df.copy()
    out["event_id"] = out.apply(_event_id_from_row, axis=1)
    summary_cols = [
        "ID_position",
        "ID_position_summary",
        "ID_position_detail",
        "ID_position_source",
        "ID_position_roles",
        "ID_position_confidence",
        "transcript_structure_roles",
        "position_evidence_relation",
        "candidate_TE_event",
        "candidate_TE_confidence",
        "ID_position_evaluable_sample_n",
        "ID_position_evaluable_replicate_n",
    ]
    out = _merge_position_summary_columns(out, position_summary_df, summary_cols)
    out["exon_id"] = out["event_id"]
    out = out.drop(columns=["event_id", "te_exon_id"], errors="ignore")
    aggregated = _aggregate_te_overlap_exon_table(out)
    out_path = os.path.join(classify_dir, f"{project_prefix}.exon.detail.tsv")
    _select_te_overlap_exon_output_columns(
        aggregated,
        include_detail=True,
        include_junction_detail=bool(include_junction_detail),
    ).to_csv(out_path, sep="\t", index=False)
    return out_path


def add_position_summary_to_teexon(
    teexon_path,
    position_summary_df,
    keep_detail=False,
    include_junction_detail=False,
):
    """Add quant-required positional evidence columns to the default TE-exon table."""
    if not os.path.isfile(teexon_path):
        return teexon_path
    try:
        te_df = pd.read_csv(teexon_path, sep="\t")
    except (OSError, EmptyDataError, ParserError, UnicodeDecodeError) as exc:
        raise RuntimeError(f"Failed to read TE-exon table: {teexon_path}") from exc
    if te_df.empty:
        te_df.to_csv(teexon_path, sep="\t", index=False)
        return teexon_path
    required = {"exon_chrom", "exon_start", "exon_end", "exon_strand"}
    missing = sorted(required - set(te_df.columns))
    if missing:
        raise RuntimeError(f"TE-exon table missing column(s): {', '.join(missing)}")
    out = te_df.copy()
    out["event_id"] = (
        out["exon_chrom"].astype(str)
        + ":"
        + pd.to_numeric(out["exon_start"], errors="coerce").fillna(-1).astype(int).astype(str)
        + "-"
        + pd.to_numeric(out["exon_end"], errors="coerce").fillna(-1).astype(int).astype(str)
        + ":"
        + out["exon_strand"].astype(str)
    )
    summary_cols = [
        "ID_position",
        "ID_position_summary",
        "ID_position_detail",
        "ID_position_source",
        "ID_position_roles",
        "ID_position_confidence",
        "transcript_structure_roles",
        "position_evidence_relation",
        "candidate_TE_event",
        "candidate_TE_confidence",
        "ID_position_evaluable_sample_n",
        "ID_position_evaluable_replicate_n",
    ]
    out = _merge_position_summary_columns(out, position_summary_df, summary_cols)
    out["exon_id"] = out["event_id"]
    out = out.drop(columns=["event_id", "te_exon_id"], errors="ignore")
    lead_cols = [col for col in ["exon_id", "metaexon_id", "transcript_exon_id"] if col in out.columns]
    out = out[lead_cols + [col for col in out.columns if col not in lead_cols]]
    aggregated = _aggregate_te_overlap_exon_table(out)
    if keep_detail:
        root, ext = os.path.splitext(teexon_path)
        detail_path = f"{root}.detail{ext or '.tsv'}"
        _select_te_overlap_exon_output_columns(
            aggregated,
            include_detail=True,
            include_junction_detail=bool(include_junction_detail),
        ).to_csv(detail_path, sep="\t", index=False)
    _select_te_overlap_exon_output_columns(aggregated).to_csv(teexon_path, sep="\t", index=False)
    return teexon_path


def _build_metaexon_annotation_summary(annotation_df):
    if annotation_df.empty or "metaexon_id" not in annotation_df.columns:
        return {}
    work = annotation_df.copy()
    work["metaexon_id"] = work["metaexon_id"].fillna("").astype(str).str.strip()
    work = work.loc[work["metaexon_id"] != ""].copy()
    if work.empty:
        return {}
    out = {}
    for metaexon_id, group in work.groupby("metaexon_id", sort=True):
        raw_pass = pd.to_numeric(group.get("te_overlap_pass_raw", 0), errors="coerce").fillna(0).astype(int)
        final_col = "te_overlap_pass_final" if "te_overlap_pass_final" in group.columns else "te_overlap_pass_raw"
        final_pass = pd.to_numeric(group.get(final_col, 0), errors="coerce").fillna(0).astype(int)
        degraded = pd.to_numeric(group.get("te_overlap_degraded", 0), errors="coerce").fillna(0).astype(int) if "te_overlap_degraded" in group.columns else pd.Series([0] * len(group))
        te_rows = group.loc[raw_pass == 1].copy()
        out[str(metaexon_id)] = {
            "matched_transcript_exon_ids": _join_clean_tokens(group.get("transcript_exon_id", [])),
            "matched_transcript_ids": _join_clean_tokens(group.get("transcript_id", [])),
            "matched_te_exon_ids": _join_clean_tokens(te_rows.apply(_event_id_from_row, axis=1)) if not te_rows.empty else "",
            "matched_te_overlap_raw_n": int(raw_pass.sum()),
            "matched_te_overlap_final_n": int(final_pass.sum()),
            "matched_te_overlap_degraded_n": int(degraded.sum()),
            "te_overlap_label": "TE_overlap" if int(final_pass.max()) == 1 else "no_overlap",
            "te_overlap_n": int(pd.to_numeric(group.get("te_overlap_n", 0), errors="coerce").fillna(0).max()),
            "te_overlap_bp_max": int(pd.to_numeric(group.get("te_overlap_bp_max", 0), errors="coerce").fillna(0).max()),
            "te_overlap_frac_max": float(pd.to_numeric(group.get("te_overlap_frac_max", 0.0), errors="coerce").fillna(0.0).max()),
            "te_boundary_hit_any": int(pd.to_numeric(group.get("te_boundary_hit_any", 0), errors="coerce").fillna(0).max()),
            "te_boundary_side": _collapse_boundary_side(group.get("te_boundary_side", [])),
            "junction_te_side_reads_max": float(
                pd.to_numeric(group.get("junction_te_side_reads_max", 0.0), errors="coerce").fillna(0.0).max()
            ),
            "te_overlap_pass_any": int(final_pass.max()),
            "te_splice_site_repeat_TE": _join_clean_tokens(
                token
                for value in group.get("te_splice_site_repeat_TE", [])
                for token in _split_clean_tokens(value)
            ),
            "te_other_overlap_TE": _join_clean_tokens(
                token
                for value in group.get("te_other_overlap_TE", [])
                for token in _split_clean_tokens(value)
            ),
        }
    return out


def _add_metaexon_summary_columns(df, annotation_df, exon_col="exon", strand_col="strand"):
    if df.empty:
        return df
    summary_map = _build_metaexon_annotation_summary(annotation_df)
    if not summary_map:
        return df
    out = df.copy()
    out["metaexon_id"] = out[exon_col].astype(str) + ":" + out[strand_col].astype(str)
    summary_cols = [
        "matched_transcript_exon_ids",
        "matched_transcript_ids",
        "matched_te_exon_ids",
        "matched_te_overlap_raw_n",
        "matched_te_overlap_final_n",
        "matched_te_overlap_degraded_n",
        "te_overlap_label",
        "te_overlap_n",
        "te_overlap_bp_max",
        "te_overlap_frac_max",
        "te_boundary_hit_any",
        "te_boundary_side",
        "junction_te_side_reads_max",
        "te_overlap_pass_any",
        "te_splice_site_repeat_TE",
        "te_other_overlap_TE",
    ]
    for col in summary_cols:
        default = 0 if col.endswith("_n") or col in {"te_overlap_n", "te_overlap_bp_max", "te_boundary_hit_any", "te_overlap_pass_any"} else ""
        if col == "te_overlap_frac_max":
            default = 0.0
        if col == "junction_te_side_reads_max":
            default = 0.0
        out[col] = out["metaexon_id"].map(lambda k: summary_map.get(str(k), {}).get(col, default))
    return out


def export_project_exon_summary(classify_dir, project_prefix, annotation_df):
    """Combine per-replicate HITindex exon outputs into the project exon table."""
    hitindex_dir = os.path.join(classify_dir, "HITindex")
    out_path = os.path.join(classify_dir, f"{project_prefix}.exon.detail.tsv")
    rows = []
    if os.path.isdir(hitindex_dir):
        for exon_path in sorted(glob.glob(os.path.join(hitindex_dir, "*.exon"))):
            replicate = os.path.basename(exon_path).rsplit(".", 1)[0]
            try:
                df = pd.read_csv(exon_path, sep="\t")
            except (OSError, EmptyDataError, ParserError, UnicodeDecodeError) as exc:
                raise RuntimeError(f"Failed to read HITindex exon output: {exon_path}") from exc
            if df.empty:
                raise RuntimeError(f"HITindex exon output is empty: {exon_path}")
            df = df.copy()
            df.insert(0, "replicate", replicate)
            df.insert(0, "sample", _infer_sample_id(replicate))
            rows.append(df)
    if rows:
        out_df = pd.concat(rows, ignore_index=True)
        if {"exon", "strand"}.issubset(out_df.columns):
            out_df = _add_metaexon_summary_columns(out_df, annotation_df)
        else:
            raise RuntimeError(f"HITindex exon outputs in {hitindex_dir} must include exon and strand columns.")
    else:
        raise RuntimeError(f"No HITindex exon outputs were found in {hitindex_dir}.")
    out_df.to_csv(out_path, sep="\t", index=False)
    return out_path


def enrich_afe_ale_outputs(classify_dir, project_prefix, annotation_df):
    """Add metaexon TE annotation columns to combined AFE/ALE PSI outputs."""
    for suffix in ["AFEPSI", "ALEPSI"]:
        path = os.path.join(classify_dir, f"{project_prefix}.{suffix}")
        if not os.path.isfile(path):
            raise RuntimeError(f"Expected combined {suffix} output is missing: {path}")
        try:
            df = pd.read_csv(path, sep="\t")
        except (OSError, EmptyDataError, ParserError, UnicodeDecodeError) as exc:
            raise RuntimeError(f"Failed to read combined {suffix} output: {path}") from exc
        if not {"exon", "strand"}.issubset(df.columns):
            raise RuntimeError(f"Combined {suffix} output must include exon and strand columns: {path}")
        df = _add_metaexon_summary_columns(df, annotation_df)
        df.to_csv(path, sep="\t", index=False)


def publish_afe_ale_outputs(classify_dir, project_prefix, keep_detail=False, include_junction_detail=False):
    """Publish AFE/ALE outputs with distinct default and detail filenames."""
    for suffix in ["AFEPSI", "ALEPSI"]:
        path = os.path.join(classify_dir, f"{project_prefix}.{suffix}")
        if not os.path.isfile(path):
            continue
        try:
            df = pd.read_csv(path, sep="\t")
        except (OSError, EmptyDataError, ParserError, UnicodeDecodeError) as exc:
            raise RuntimeError(f"Failed to read combined {suffix} output: {path}") from exc
        pass_col = None
        for col in ["te_overlap_pass_any", "matched_te_overlap_final_n", "matched_te_overlap_raw_n"]:
            if col in df.columns:
                pass_col = col
                break
        if pass_col is None:
            raise RuntimeError(f"Combined {suffix} output lacks TE-overlap columns: {path}")
        published = _clarify_afe_ale_columns(df)
        if keep_detail:
            detail_path = os.path.join(classify_dir, f"{project_prefix}.{suffix}.detail.tsv")
            full_published = _aggregate_afe_ale_te_overlap_table(published)
            _select_afe_ale_te_overlap_output_columns(
                full_published,
                include_detail=True,
                include_junction_detail=bool(include_junction_detail),
            ).to_csv(detail_path, sep="\t", index=False)
        keep = pd.to_numeric(df[pass_col], errors="coerce").fillna(0).astype(float) > 0
        te_path = os.path.join(classify_dir, f"{project_prefix}.TE_overlap.{suffix}.tsv")
        te_overlap_published = _aggregate_afe_ale_te_overlap_table(published.loc[keep])
        if keep_detail:
            te_detail_path = os.path.join(classify_dir, f"{project_prefix}.TE_overlap.{suffix}.detail.tsv")
            _select_afe_ale_te_overlap_output_columns(
                te_overlap_published,
                include_detail=True,
                include_junction_detail=bool(include_junction_detail),
            ).to_csv(te_detail_path, sep="\t", index=False)
        _select_afe_ale_te_overlap_output_columns(te_overlap_published).to_csv(te_path, sep="\t", index=False)
        try:
            os.remove(path)
        except FileNotFoundError:
            pass
        try:
            os.remove(os.path.join(classify_dir, f"{project_prefix}.{suffix}.lock"))
        except FileNotFoundError:
            pass


def _clarify_afe_ale_columns(df):
    """Use explicit event-level column names in published AFE/ALE tables."""
    out = df.copy()
    rename_map = {}
    if "gene" in out.columns:
        rename_map["gene"] = "gene_id"
    if "exon" in out.columns:
        rename_map["exon"] = "exon_coord"
    if "strand" in out.columns:
        rename_map["strand"] = "exon_strand"
    if "label" in out.columns and "te_overlap_label" in out.columns:
        out = out.drop(columns=["label"])
    elif "label" in out.columns:
        rename_map["label"] = "te_overlap_label"
    out = out.rename(columns=rename_map)
    preferred = [
        "gene_id",
        "metaexon_id",
        "exon_coord",
        "exon_strand",
        "te_overlap_label",
    ]
    ordered = [col for col in preferred if col in out.columns]
    ordered.extend([col for col in out.columns if col not in ordered])
    return out[ordered]


def _aggregate_afe_ale_te_overlap_table(df):
    if df.empty:
        return df
    rows = []
    for _, row in df.iterrows():
        exon_ids = _split_clean_tokens(row.get("matched_te_exon_ids", ""))
        if not exon_ids:
            metaexon_id = str(row.get("metaexon_id", "")).strip()
            if metaexon_id:
                exon_ids = [metaexon_id]
        for exon_id in exon_ids:
            rec = row.to_dict()
            rec["exon_id"] = exon_id
            rows.append(rec)
    if not rows:
        return pd.DataFrame(columns=["exon_id"])

    expanded = pd.DataFrame(rows)
    metadata_cols = {
        "exon_id",
        "gene_id",
        "metaexon_id",
        "exon_coord",
        "exon_strand",
        "te_overlap_label",
        "te_overlap_n",
        "te_overlap_bp_max",
        "te_overlap_frac_max",
        "te_boundary_hit_any",
        "te_boundary_side",
        "junction_te_side_reads_max",
        "te_overlap_pass_any",
        "te_splice_site_repeat_TE",
        "te_other_overlap_TE",
        "matched_transcript_exon_ids",
        "matched_transcript_ids",
        "matched_te_exon_ids",
        "matched_te_overlap_raw_n",
        "matched_te_overlap_final_n",
        "matched_te_overlap_degraded_n",
    }
    measurement_cols = [col for col in expanded.columns if col not in metadata_cols]

    out_rows = []
    for exon_id, group in expanded.groupby("exon_id", sort=True):
        out = {
            "exon_id": exon_id,
            "gene_id": _join_split_tokens(group.get("gene_id", [])),
            "metaexon_id": _join_split_tokens(group.get("metaexon_id", [])),
            "transcript_id": _join_split_tokens(group.get("matched_transcript_ids", [])),
            "transcript_exon_id": _join_split_tokens(group.get("matched_transcript_exon_ids", [])),
            "te_overlap_label": _collapse_te_overlap_label(group.get("te_overlap_label", [])),
            "te_overlap_n": _max_numeric_value(group.get("te_overlap_n", []), default=0),
            "te_overlap_bp_max": _max_numeric_value(group.get("te_overlap_bp_max", []), default=0),
            "te_overlap_frac_max": _max_numeric_value(group.get("te_overlap_frac_max", []), default=0.0),
            "te_boundary_hit_any": _max_numeric_value(group.get("te_boundary_hit_any", []), default=0),
            "te_boundary_side": _collapse_boundary_side(group.get("te_boundary_side", [])),
            "junction_te_side_reads_max": _max_numeric_value(
                group.get("junction_te_side_reads_max", []),
                default=0.0,
            ),
            "te_splice_site_repeat_TE": _join_split_tokens(group.get("te_splice_site_repeat_TE", [])),
            "te_other_overlap_TE": _join_split_tokens(group.get("te_other_overlap_TE", [])),
        }
        for col in measurement_cols:
            out[col] = _join_split_tokens(group[col])
        out_rows.append(out)

    out_df = pd.DataFrame(out_rows)
    preferred = [
        "exon_id",
        "gene_id",
        "metaexon_id",
        "transcript_id",
        "transcript_exon_id",
        "te_overlap_label",
        "te_overlap_n",
        "te_overlap_bp_max",
        "te_overlap_frac_max",
        "te_boundary_hit_any",
        "te_boundary_side",
        "junction_te_side_reads_max",
        "te_splice_site_repeat_TE",
        "te_other_overlap_TE",
    ]
    ordered = [col for col in preferred if col in out_df.columns]
    ordered.extend([col for col in out_df.columns if col not in ordered])
    return out_df[ordered]


def _select_afe_ale_te_overlap_output_columns(df, *, include_detail=False, include_junction_detail=False):
    if df.empty:
        return df
    metadata_cols = {
        "exon_id",
        "metaexon_id",
        "gene_id",
        "transcript_id",
        "transcript_exon_id",
        "te_overlap_label",
        "te_splice_site_repeat_TE",
        "te_overlap_bp_max",
        "te_overlap_frac_max",
        "te_boundary_side",
        "junction_te_side_reads_max",
        "te_overlap_n",
        "te_boundary_hit_any",
        "te_overlap_pass_any",
        "te_other_overlap_TE",
    }
    measurement_cols = [col for col in df.columns if col not in metadata_cols]
    default_cols = [
        "exon_id",
        "metaexon_id",
        "gene_id",
        "transcript_id",
        "transcript_exon_id",
        "te_overlap_label",
        "te_splice_site_repeat_TE",
    ]
    detail_cols = [
        "te_overlap_bp_max",
        "te_overlap_frac_max",
        "te_boundary_side",
    ]
    if include_junction_detail:
        detail_cols.append("junction_te_side_reads_max")
    selected = default_cols + (detail_cols if include_detail else []) + measurement_cols
    return df[[col for col in selected if col in df.columns]].copy()

"""Annotate TE-overlap evidence for qual-mode exon events."""

import collections
import os

import pandas as pd
import pybedtools


def _normalize_int(value, default=0):
    try:
        if pd.isna(value):
            return int(default)
        return int(float(value))
    except (TypeError, ValueError):
        return int(default)


def _build_exon_event_key(exon, strand):
    return f"{str(exon)}|{str(strand)}"


def _format_te_set(values):
    if not values:
        return ""
    return ",".join(sorted(values))


def _collect_te_overlap_event_annotations(
    transcript_exon_rows,
    te_bed_path,
    min_overlap_bp=10,
    min_overlap_frac=0.1,
    boundary_flank_bp=10,
):
    if not transcript_exon_rows:
        return {}, set()
    event_rows = []
    for row in transcript_exon_rows:
        chrom = str(row.get("exon_chrom", ""))
        strand = str(row.get("exon_strand", ""))
        try:
            start = int(row.get("exon_start", -1))
            end = int(row.get("exon_end", -1))
        except (TypeError, ValueError):
            continue
        if not chrom or end <= start:
            continue
        event_rows.append(
            {
                "exon": f"{chrom}:{start}-{end}",
                "strand": strand,
            }
        )
    if not event_rows:
        return {}, set()

    event_df = pd.DataFrame(event_rows).drop_duplicates()
    anno_map = _collect_exon_te_annotations(
        event_df,
        te_bed_path=te_bed_path,
        min_overlap_bp=min_overlap_bp,
        min_overlap_frac=min_overlap_frac,
        boundary_flank_bp=boundary_flank_bp,
    )
    te_overlap_keys = {
        k for k, v in anno_map.items()
        if int(v.get("te_overlap_pass_any", 0)) == 1 or str(v.get("label", "")) == "TE_overlap"
    }
    return anno_map, te_overlap_keys


def _split_metaexon_id(metaexon_id):
    value = str(metaexon_id or "").strip()
    if not value or ":" not in value:
        return "", ""
    exon, strand = value.rsplit(":", 1)
    return exon, strand


def _update_te_string_set(out_set, value):
    for token in str(value or "").split(","):
        token = token.strip()
        if token:
            out_set.add(token)


def _build_metaexon_te_annotation_map(transcript_exon_rows, te_event_anno_map):
    """Aggregate updated transcript-exon TE annotations back to metaexon-level event keys."""
    if not transcript_exon_rows or not te_event_anno_map:
        return {}

    agg = {}
    for row in transcript_exon_rows:
        exon = f'{row["exon_chrom"]}:{int(row["exon_start"])}-{int(row["exon_end"])}'
        exon_strand = str(row["exon_strand"])
        transcript_key = _build_exon_event_key(exon, exon_strand)
        ann = te_event_anno_map.get(transcript_key)
        if ann is None:
            continue

        meta_exon, meta_strand = _split_metaexon_id(row.get("metaexon_id", ""))
        if not meta_exon or not meta_strand:
            continue
        meta_key = _build_exon_event_key(meta_exon, meta_strand)
        rec = agg.setdefault(
            meta_key,
            {
                "label": "no_overlap",
                "te_overlap_n": 0,
                "te_overlap_ids": set(),
                "te_overlap_bp_max": 0,
                "te_overlap_frac_max": 0.0,
                "te_boundary_hit_any": 0,
                "te_overlap_pass_any": 0,
                "te_splice_site_repeat_TE": set(),
                "te_other_overlap_TE": set(),
            },
        )
        pass_any = int(ann.get("te_overlap_pass_any", 0))
        rec["te_overlap_pass_any"] = max(int(rec["te_overlap_pass_any"]), pass_any)
        if pass_any:
            rec["label"] = "TE_overlap"
        _update_te_string_set(rec["te_overlap_ids"], ann.get("te_splice_site_repeat_TE", ""))
        _update_te_string_set(rec["te_overlap_ids"], ann.get("te_other_overlap_TE", ""))
        rec["te_overlap_bp_max"] = max(int(rec["te_overlap_bp_max"]), int(ann.get("te_overlap_bp_max", 0)))
        rec["te_overlap_frac_max"] = max(float(rec["te_overlap_frac_max"]), float(ann.get("te_overlap_frac_max", 0.0)))
        rec["te_boundary_hit_any"] = max(int(rec["te_boundary_hit_any"]), int(ann.get("te_boundary_hit_any", 0)))
        _update_te_string_set(rec["te_splice_site_repeat_TE"], ann.get("te_splice_site_repeat_TE", ""))
        _update_te_string_set(rec["te_other_overlap_TE"], ann.get("te_other_overlap_TE", ""))

    out = {}
    for meta_key, rec in agg.items():
        next_rec = dict(rec)
        next_rec["te_overlap_n"] = len(rec["te_overlap_ids"])
        next_rec.pop("te_overlap_ids", None)
        next_rec["te_splice_site_repeat_TE"] = _format_te_set(rec["te_splice_site_repeat_TE"])
        next_rec["te_other_overlap_TE"] = _format_te_set(rec["te_other_overlap_TE"])
        out[meta_key] = next_rec
    return out


def _parse_exon_span(exon_value):
    chrom, span = str(exon_value).split(":", 1)
    start_s, end_s = span.split("-", 1)
    return chrom, _normalize_int(start_s), _normalize_int(end_s)


def _collect_exon_te_annotations(df_all, te_bed_path, min_overlap_bp=10, min_overlap_frac=0.1, boundary_flank_bp=10):
    """
    Collect TE-overlap evidence for AFE/ALE exon events.

    Label rule is strict and unchanged in spirit:
      overlap_bp >= min_overlap_bp AND overlap_frac >= min_overlap_frac AND boundary_hit.
    """
    if df_all.empty or not os.path.isfile(te_bed_path):
        return {}

    rows = []
    for exon, strand in df_all[["exon", "strand"]].drop_duplicates().itertuples(index=False):
        if ":" not in str(exon):
            continue
        chrom, span = str(exon).split(":", 1)
        start_s, end_s = span.split("-", 1)
        start = _normalize_int(start_s)
        end = _normalize_int(end_s)
        if end <= start:
            continue
        key = _build_exon_event_key(exon, strand)
        rows.append((chrom, start, end, key, 0, str(strand)))
    if not rows:
        return {}

    event_bt = pybedtools.BedTool(rows)
    te_bt = pybedtools.BedTool(te_bed_path)
    overlaps = event_bt.intersect(te_bt, wa=True, wb=True, s=True)

    agg = collections.defaultdict(
        lambda: {
            "all_te": set(),
            "splice_te": set(),
            "other_te": set(),
            "pass_te": set(),
            "max_overlap_bp": 0,
            "max_overlap_frac": 0.0,
            "boundary_hit_any": 0,
            "overlap_pass_any": 0,
            "pass_left_boundary_hit_any": 0,
            "pass_right_boundary_hit_any": 0,
        }
    )

    for o in overlaps:
        fields = o.fields
        key = fields[3]
        event_start = _normalize_int(fields[1])
        event_end = _normalize_int(fields[2])
        te_start = _normalize_int(fields[7])
        te_end = _normalize_int(fields[8])
        te_id = str(fields[9])

        overlap_bp = min(event_end, te_end) - max(event_start, te_start)
        if overlap_bp <= 0:
            continue
        event_len = event_end - event_start
        overlap_frac = float(overlap_bp) / float(event_len) if event_len > 0 else 0.0
        left_boundary_hit = te_start <= event_start + int(boundary_flank_bp)
        right_boundary_hit = te_end >= event_end - int(boundary_flank_bp)
        boundary_hit = left_boundary_hit or right_boundary_hit
        overlap_core = (overlap_bp >= int(min_overlap_bp)) and (overlap_frac >= float(min_overlap_frac))
        overlap_pass = overlap_core and boundary_hit

        rec = agg[key]
        if overlap_core:
            rec["all_te"].add(te_id)
            rec["max_overlap_bp"] = max(rec["max_overlap_bp"], overlap_bp)
            rec["max_overlap_frac"] = max(rec["max_overlap_frac"], overlap_frac)
            if boundary_hit:
                rec["splice_te"].add(te_id)
                rec["boundary_hit_any"] = 1
            else:
                rec["other_te"].add(te_id)
            if overlap_pass:
                rec["pass_te"].add(te_id)
                rec["overlap_pass_any"] = 1
                if left_boundary_hit:
                    rec["pass_left_boundary_hit_any"] = 1
                if right_boundary_hit:
                    rec["pass_right_boundary_hit_any"] = 1

    out = {}
    for key, rec in agg.items():
        out[key] = {
            "label": "TE_overlap" if rec["overlap_pass_any"] else "no_overlap",
            "TE_info": _format_te_set(rec["pass_te"]),
            "te_overlap_n": len(rec["all_te"]),
            "te_overlap_bp_max": int(rec["max_overlap_bp"]),
            "te_overlap_frac_max": float(rec["max_overlap_frac"]),
            "te_boundary_hit_any": int(rec["boundary_hit_any"]),
            "te_overlap_pass_any": int(rec["overlap_pass_any"]),
            "te_splice_site_repeat_TE": _format_te_set(rec["splice_te"]),
            "te_other_overlap_TE": _format_te_set(rec["other_te"]),
            "te_left_boundary_hit_any": int(rec["pass_left_boundary_hit_any"]),
            "te_right_boundary_hit_any": int(rec["pass_right_boundary_hit_any"]),
        }
    return out

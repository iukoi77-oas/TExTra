"""Adapt qual-stage annotation and HITindex outputs for quant mode."""

import os
from collections import defaultdict

import pandas as pd

from util.common.define_layout import CLASSIFY_DIR
from util.common.write_logs import log_message

TE_CANDIDATE_BOUNDARY_TOLERANCE_BP = 10
HITINDEX_POSITION_ROLES = {
    "first": (("first",), "high"),
    "internal": (("internal",), "high"),
    "last": (("last",), "high"),
    "FirstInternal_high": (("first", "internal"), "high"),
    "FirstInternal_medium": (("first", "internal"), "medium"),
    "InternalLast_high": (("internal", "last"), "high"),
    "InternalLast_medium": (("internal", "last"), "medium"),
}
POSITION_ROLE_ORDER = {"first": 0, "internal": 1, "last": 2}
CONFIDENCE_ORDER = {"high": 0, "medium": 1, "low": 2, "structure_only": 3, "NA": 4}
CANDIDATE_EVENT_ORDER = {
    "candidate_TSS": 0,
    "candidate_TES": 1,
    "candidate_exonized": 2,
    "NA": 3,
}
_LOGGED_AMBIGUOUS_MULTI_GENE_REPORTS = set()


def _build_exon_event_id(chrom, start, end, strand):
    return f"{chrom}:{int(start)}-{int(end)}:{strand}"


def _parse_exon_coord(exon, strand=None):
    text = str(exon).strip()
    if not text:
        raise ValueError("empty exon coordinate")
    if text.count(":") >= 2:
        left, parsed_strand = text.rsplit(":", 1)
        chrom, span = left.split(":", 1)
        use_strand = parsed_strand
    else:
        chrom, span = text.split(":", 1)
        use_strand = strand
    start_s, end_s = span.split("-", 1)
    if not use_strand:
        raise ValueError(f"missing strand for exon coordinate: {exon}")
    return chrom, int(start_s), int(end_s), str(use_strand)


def _join_unique_values(values):
    vals = sorted({str(v).strip() for v in values if str(v).strip() and str(v).strip().lower() != "nan"})
    return ",".join(vals)


def _ambiguous_multi_gene_report_path(qual_dir, project):
    return os.path.join(
        qual_dir,
        CLASSIFY_DIR,
        f"{project}.ambiguous_multi_gene_exon_events.tsv",
    )


def _event_id_series(df):
    return df.apply(
        lambda row: _build_exon_event_id(
            row["exon_chrom"],
            row["exon_start"],
            row["exon_end"],
            row["exon_strand"],
        ),
        axis=1,
    )


def _find_ambiguous_multi_gene_events(df):
    if df.empty:
        return {}
    work = df.copy()
    if "event_id" not in work.columns:
        work["event_id"] = _event_id_series(work)
    ambiguous = {}
    for event_id, group in work.groupby("event_id", sort=True):
        genes = sorted({str(v).strip() for v in group["gene_id"] if str(v).strip()})
        if len(genes) <= 1:
            continue
        txs = sorted({str(v).strip() for v in group["transcript_id"] if str(v).strip()})
        ambiguous[str(event_id)] = {
            "exon_event": str(event_id),
            "gene_ids": ",".join(genes),
            "n_genes": len(genes),
            "transcript_ids": ",".join(txs),
            "n_transcripts": len(txs),
            "action": "excluded_from_quant",
            "reason": "exon_event_maps_to_multiple_genes",
        }
    return ambiguous


def _write_ambiguous_multi_gene_report(qual_dir, project, ambiguous_events):
    if not ambiguous_events:
        return None
    report_path = _ambiguous_multi_gene_report_path(qual_dir, project)
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    cols = [
        "exon_event",
        "gene_ids",
        "n_genes",
        "transcript_ids",
        "n_transcripts",
        "action",
        "reason",
    ]
    pd.DataFrame(
        [ambiguous_events[key] for key in sorted(ambiguous_events)],
        columns=cols,
    ).to_csv(report_path, sep="\t", index=False)
    return report_path


def _log_ambiguous_multi_gene_exclusion(ambiguous_events, report_path):
    if not ambiguous_events:
        return
    if report_path in _LOGGED_AMBIGUOUS_MULTI_GENE_REPORTS:
        return
    _LOGGED_AMBIGUOUS_MULTI_GENE_REPORTS.add(report_path)
    log_message(
        "[WARNING]",
        (
            f"Excluded {len(ambiguous_events)} ambiguous multi-gene TE-overlap exon event(s) "
            f"from quantification; report: {report_path}"
        ),
        color="warning",
    )


def _qual_annotation_table_path(qual_dir, project):
    new_path = os.path.join(qual_dir, CLASSIFY_DIR, "transcript_exon_te_annotation.tsv")
    if os.path.isfile(new_path):
        return new_path
    project_path = os.path.join(qual_dir, CLASSIFY_DIR, f"{project}.transcript_exon_te_annotation.tsv")
    if os.path.isfile(project_path):
        return project_path
    return os.path.join(
        qual_dir,
        CLASSIFY_DIR,
        "annotation",
        f"{project}.transcript_exon_te_annotation.tsv",
    )


def _qual_teexon_table_path(qual_dir, project):
    new_path = os.path.join(qual_dir, CLASSIFY_DIR, f"{project}.TE_overlap.exon.tsv")
    if os.path.isfile(new_path):
        return new_path
    return os.path.join(qual_dir, CLASSIFY_DIR, f"{project}.TEexon.tsv")


def _qual_exon_detail_table_path(qual_dir, project):
    new_path = os.path.join(qual_dir, CLASSIFY_DIR, f"{project}.exon.detail.tsv")
    if os.path.isfile(new_path):
        return new_path
    return os.path.join(qual_dir, CLASSIFY_DIR, f"{project}.exon")


def _load_te_exon_annotation_df(qual_dir, project):
    annotation_path = _qual_annotation_table_path(qual_dir, project)
    if not os.path.isfile(annotation_path):
        raise FileNotFoundError(f"qual TE-exon annotation table not found: {annotation_path}")
    try:
        df = pd.read_csv(annotation_path, sep="\t")
    except (OSError, ValueError, pd.errors.EmptyDataError, pd.errors.ParserError, UnicodeDecodeError) as exc:
        raise RuntimeError(f"Failed to read qual TE-exon annotation table: {annotation_path}") from exc
    required = {"transcript_id", "gene_id", "exon_chrom", "exon_start", "exon_end", "exon_strand"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise RuntimeError(
            f"Invalid qual TE-exon annotation table: {annotation_path} missing column(s): {', '.join(missing)}"
        )
    if df.empty:
        return df
    work = df.copy()
    if "te_overlap_pass_final" in work.columns:
        pass_col = "te_overlap_pass_final"
    elif "te_overlap_pass_any" in work.columns:
        pass_col = "te_overlap_pass_any"
    elif "te_overlap_pass_raw" in work.columns:
        pass_col = "te_overlap_pass_raw"
    elif "te_overlap_label" in work.columns:
        pass_col = "te_overlap_pass_any"
        work[pass_col] = (
            work["te_overlap_label"].fillna("").astype(str).str.strip().str.lower() == "te_overlap"
        ).astype(int)
    else:
        raise RuntimeError(
            f"Invalid qual TE-exon annotation table: {annotation_path} missing TE overlap pass/label column."
        )
    work[pass_col] = pd.to_numeric(work[pass_col], errors="coerce").fillna(0).astype(int)
    work = work.loc[work[pass_col] == 1].copy()
    if work.empty:
        return work
    work["transcript_id"] = work["transcript_id"].fillna("").astype(str).str.strip()
    if (work["transcript_id"] == "").any():
        raise RuntimeError(
            f"Invalid qual TE-exon annotation table: {annotation_path} has TE-overlap row(s) with empty transcript_id."
        )
    work["exon_start"] = pd.to_numeric(work["exon_start"], errors="coerce")
    work["exon_end"] = pd.to_numeric(work["exon_end"], errors="coerce")
    bad_coord = work["exon_start"].isna() | work["exon_end"].isna() | (work["exon_end"] <= work["exon_start"])
    if bad_coord.any():
        raise RuntimeError(
            f"Invalid qual TE-exon annotation table: {annotation_path} has malformed exon_start/exon_end."
        )
    work["exon_start"] = work["exon_start"].astype(int)
    work["exon_end"] = work["exon_end"].astype(int)
    if "te_exon_id" not in work.columns:
        if "exon_id" in work.columns:
            work["te_exon_id"] = work["exon_id"].astype(str)
        else:
            work["te_exon_id"] = (
                work["exon_chrom"].astype(str)
                + ":"
                + work["exon_start"].astype(str)
                + "-"
                + work["exon_end"].astype(str)
                + ":"
                + work["exon_strand"].astype(str)
            )
    work["te_exon_id"] = work["te_exon_id"].fillna("").astype(str).str.strip()
    work = work.loc[work["te_exon_id"] != ""].copy()
    return work


def _load_te_exon_event_tx_map(qual_dir, project):
    df = _load_te_exon_annotation_df(qual_dir, project)
    event_to_txs = defaultdict(set)
    event_to_gene = {}
    event_ids = set()
    ambiguous_events = _find_ambiguous_multi_gene_events(df)
    ambiguous_event_ids = set(ambiguous_events)
    report_path = _write_ambiguous_multi_gene_report(qual_dir, project, ambiguous_events)
    _log_ambiguous_multi_gene_exclusion(ambiguous_events, report_path)
    if df.empty:
        return [], event_to_txs, event_to_gene
    for _, row in df.iterrows():
        event_id = _build_exon_event_id(
            row.get("exon_chrom", ""),
            row.get("exon_start", -1),
            row.get("exon_end", -1),
            row.get("exon_strand", ""),
        )
        tx_id = str(row.get("transcript_id", "")).strip()
        gene_id = str(row.get("gene_id", "")).strip()
        if not event_id:
            continue
        if event_id in ambiguous_event_ids:
            continue
        event_ids.add(event_id)
        if tx_id:
            event_to_txs[event_id].add(tx_id)
        if gene_id:
            prev_gene = str(event_to_gene.get(event_id, "")).strip()
            if not prev_gene:
                event_to_gene[event_id] = gene_id
            elif prev_gene != gene_id:
                continue
    return sorted(event_ids), event_to_txs, event_to_gene


def _parse_te_interval_tokens(*values):
    intervals = []
    for value in values:
        for token in str(value).split(","):
            token = token.strip()
            if not token or token.lower() == "nan":
                continue
            parts = token.split("|")
            if len(parts) < 3:
                continue
            try:
                intervals.append((parts[0], int(float(parts[1])), int(float(parts[2]))))
            except ValueError:
                continue
    return intervals


def _te_overlaps_boundary(te_intervals, chrom, boundary, tolerance_bp=TE_CANDIDATE_BOUNDARY_TOLERANCE_BP):
    for te_chrom, te_start, te_end in te_intervals:
        if te_chrom == chrom and te_start <= boundary + tolerance_bp and te_end >= boundary - tolerance_bp:
            return True
    return False


def _te_covers_exon(te_intervals, chrom, start, end, tolerance_bp=TE_CANDIDATE_BOUNDARY_TOLERANCE_BP):
    for te_chrom, te_start, te_end in te_intervals:
        if te_chrom != chrom:
            continue
        if te_start <= start + tolerance_bp and te_end >= end - tolerance_bp:
            return True
    return False


def _roles_to_string(roles):
    clean = sorted({str(r).strip() for r in roles if str(r).strip() and str(r).strip() != "NA"},
                   key=lambda x: POSITION_ROLE_ORDER.get(x, 99))
    return ",".join(clean) if clean else "NA"


def _parse_roles_text(value):
    return {
        token.strip()
        for token in str(value).split(",")
        if token.strip() and token.strip() != "NA" and token.strip().lower() != "nan"
    }


def _parse_hitindex_position(position):
    raw = str(position).strip()
    if not raw or raw.lower() == "nan":
        return set(), "NA"
    if raw not in HITINDEX_POSITION_ROLES:
        return set(), "NA"
    roles, strength = HITINDEX_POSITION_ROLES[raw]
    return set(roles), strength


def _choose_candidate_summary(candidate_values):
    clean = [str(v).strip() for v in candidate_values if str(v).strip() and str(v).strip() != "NA"]
    if not clean:
        return "NA"
    return ",".join(sorted(set(clean), key=lambda x: CANDIDATE_EVENT_ORDER.get(x, 99)))


def _best_confidence(values):
    clean = [str(v).strip() for v in values if str(v).strip() and str(v).strip() != "NA"]
    if not clean:
        return "NA"
    return sorted(set(clean), key=lambda x: CONFIDENCE_ORDER.get(x, 99))[0]


def _format_role_confidence(role_to_conf):
    clean = {
        role: conf
        for role, conf in role_to_conf.items()
        if str(conf).strip() and str(conf).strip() != "NA"
    }
    if not clean:
        return "NA"
    confidences = set(clean.values())
    if len(confidences) == 1:
        return next(iter(confidences))
    return ";".join(f"{role}:{clean[role]}" for role in sorted(clean, key=lambda x: POSITION_ROLE_ORDER.get(x, 99)))


def _format_candidate_confidence(candidate_to_conf):
    clean = {
        candidate: conf
        for candidate, conf in candidate_to_conf.items()
        if str(conf).strip() and str(conf).strip() != "NA"
    }
    if not clean:
        return "NA"
    return ";".join(
        f"{candidate}:{clean[candidate]}"
        for candidate in sorted(clean, key=lambda x: CANDIDATE_EVENT_ORDER.get(x, 99))
    )


def _evidence_relation(hit_roles, structure_roles):
    hit_roles = set(hit_roles)
    structure_roles = set(structure_roles)
    if hit_roles and structure_roles:
        if hit_roles == structure_roles:
            return "consistent"
        if hit_roles & structure_roles:
            return "partially_consistent"
        return "conflict"
    if hit_roles:
        return "hitindex_only"
    if structure_roles:
        return "transcript_structure_only"
    return "NA"


def _classify_candidate_te_events(roles, chrom, start, end, strand, te_intervals):
    if not te_intervals or not roles:
        return "NA"
    candidates = []
    if "first" in roles:
        five_prime_boundary = start if strand == "+" else end
        if _te_overlaps_boundary(te_intervals, chrom, five_prime_boundary):
            candidates.append("candidate_TSS")
    if "last" in roles:
        three_prime_boundary = end if strand == "+" else start
        if _te_overlaps_boundary(te_intervals, chrom, three_prime_boundary):
            candidates.append("candidate_TES")
    if "internal" in roles and _te_covers_exon(te_intervals, chrom, start, end):
        candidates.append("candidate_exonized")
    return _choose_candidate_summary(candidates)


def _candidate_role_map(candidate_text):
    mapping = {}
    candidates = {
        token.strip()
        for token in str(candidate_text).split(",")
        if token.strip() and token.strip() != "NA" and token.strip().lower() != "nan"
    }
    if "candidate_TSS" in candidates:
        mapping["candidate_TSS"] = "first"
    if "candidate_TES" in candidates:
        mapping["candidate_TES"] = "last"
    if "candidate_exonized" in candidates:
        mapping["candidate_exonized"] = "internal"
    return mapping


def _transcript_structure_roles_from_row(row):
    try:
        exon_number = int(float(row.get("exon_number", "")))
        exon_count = int(float(row.get("transcript_exon_count", "")))
    except (TypeError, ValueError):
        return set()
    if exon_number <= 0 or exon_count <= 0:
        return set()
    if exon_count == 1:
        return {"first", "last"}
    strand = str(row.get("transcript_strand", "") or row.get("exon_strand", "")).strip()
    if strand == "-":
        if exon_number == 1:
            return {"last"}
        if exon_number == exon_count:
            return {"first"}
    else:
        if exon_number == 1:
            return {"first"}
        if exon_number == exon_count:
            return {"last"}
    return {"internal"}


def _load_transcript_structure_position_table(qual_dir, project):
    df = _load_te_exon_annotation_df(qual_dir, project)
    cols = ["event_id", "transcript_structure_roles", "candidate_TE_event_structure"]
    if df.empty:
        return pd.DataFrame(columns=cols)

    rows = []
    for event_id, group in df.groupby(
        df.apply(
            lambda row: _build_exon_event_id(
                row["exon_chrom"],
                row["exon_start"],
                row["exon_end"],
                row["exon_strand"],
            ),
            axis=1,
        ),
        sort=True,
    ):
        roles = set()
        for _, row in group.iterrows():
            roles.update(_transcript_structure_roles_from_row(row))
        chrom = str(group["exon_chrom"].iloc[0])
        start = int(pd.to_numeric(group["exon_start"], errors="coerce").iloc[0])
        end = int(pd.to_numeric(group["exon_end"], errors="coerce").iloc[0])
        strand = str(group["exon_strand"].iloc[0])
        te_intervals = _parse_te_interval_tokens(
            *(list(group.get("te_splice_site_repeat_TE", pd.Series([], dtype=str)))
              + list(group.get("te_other_overlap_TE", pd.Series([], dtype=str))))
        )
        rows.append(
            {
                "event_id": str(event_id),
                "transcript_structure_roles": _roles_to_string(roles),
                "candidate_TE_event_structure": _classify_candidate_te_events(
                    roles,
                    chrom,
                    start,
                    end,
                    strand,
                    te_intervals,
                ),
            }
        )
    return pd.DataFrame(rows, columns=cols)


def _load_hitindex_position_tables(qual_dir, project, event_ids, structure_summary_df=None):
    cols_summary = [
        "event_id",
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
    cols_detail = [
        "exon",
        "sample",
        "replicate",
        "gene_id",
        "transcript_id",
        "ID_position_raw",
        "ID_position_roles_hitindex",
        "ID_position_strength",
        "ID_position_roles_transcript_structure",
        "ID_position_source",
        "ID_position_confidence",
        "position_evidence_relation",
        "candidate_TE_event",
        "candidate_TE_confidence",
        "n_evaluable_replicates",
        "role_support_fraction",
        "HITindex",
        "nFE",
        "nINTERNAL",
        "nLE",
        "nSINGLE",
    ]
    structure_by_event = {}
    if structure_summary_df is not None and not structure_summary_df.empty:
        for _, row in structure_summary_df.iterrows():
            event_id = str(row.get("event_id", ""))
            structure_by_event[event_id] = {
                "roles": _parse_roles_text(row.get("transcript_structure_roles", "NA")),
                "candidate": str(row.get("candidate_TE_event_structure", "NA") or "NA"),
            }

    teexon_path = _qual_teexon_table_path(qual_dir, project)
    summary_required = {
        "ID_position_summary",
        "ID_position_confidence",
        "candidate_TE_event",
        "candidate_TE_confidence",
    }
    if os.path.isfile(teexon_path):
        try:
            teexon_df = pd.read_csv(teexon_path, sep="\t")
        except (OSError, ValueError, pd.errors.EmptyDataError, pd.errors.ParserError, UnicodeDecodeError) as exc:
            raise RuntimeError(f"Failed to read qual TE-exon table: {teexon_path}") from exc
        if summary_required.issubset(teexon_df.columns):
            work = teexon_df.copy()
            if "event_id" not in work.columns:
                if "exon_id" in work.columns:
                    work["event_id"] = work["exon_id"].fillna("").astype(str).str.strip()
                else:
                    required_event_cols = {"exon_chrom", "exon_start", "exon_end", "exon_strand"}
                    missing = sorted(required_event_cols - set(work.columns))
                    if missing:
                        raise RuntimeError(
                            f"Invalid qual TE-exon table: {teexon_path} missing column(s): {', '.join(missing)}"
                        )
                    work["event_id"] = work.apply(
                        lambda row: _build_exon_event_id(
                            row["exon_chrom"],
                            row["exon_start"],
                            row["exon_end"],
                            row["exon_strand"],
                        ),
                        axis=1,
                    )
            rows = []
            by_event = {str(row.get("event_id", "")): row.to_dict() for _, row in work.iterrows()}
            for event_id in sorted(set(str(e) for e in event_ids)):
                row = by_event.get(str(event_id), {})
                position_roles = str(
                    row.get(
                        "HITindex_structure_roles",
                        row.get("ID_position_roles", row.get("ID_position_summary", "NA")),
                    )
                    or "NA"
                )
                rows.append(
                    {
                        "event_id": event_id,
                        "ID_position": str(row.get("ID_position", row.get("ID_position_summary", "NA")) or "NA"),
                        "ID_position_summary": str(row.get("ID_position_summary", "NA") or "NA"),
                        "ID_position_detail": str(row.get("ID_position_detail", "NA") or "NA"),
                        "ID_position_source": str(row.get("ID_position_source", "NA") or "NA"),
                        "ID_position_roles": position_roles,
                        "ID_position_confidence": str(row.get("ID_position_confidence", "NA") or "NA"),
                        "transcript_structure_roles": str(row.get("transcript_structure_roles", "NA") or "NA"),
                        "position_evidence_relation": str(row.get("position_evidence_relation", "NA") or "NA"),
                        "candidate_TE_event": str(row.get("candidate_TE_event", "NA") or "NA"),
                        "candidate_TE_confidence": str(row.get("candidate_TE_confidence", "NA") or "NA"),
                        "ID_position_evaluable_sample_n": int(
                            pd.to_numeric(
                                pd.Series(
                                    [
                                        row.get(
                                            "HITindex_evaluable_sample_n",
                                            row.get("ID_position_evaluable_sample_n", 0),
                                        )
                                    ]
                                ),
                                errors="coerce",
                            ).fillna(0).iloc[0]
                        ),
                        "ID_position_evaluable_replicate_n": int(
                            pd.to_numeric(
                                pd.Series(
                                    [
                                        row.get(
                                            "HITindex_evaluable_replicate_n",
                                            row.get("ID_position_evaluable_replicate_n", 0),
                                        )
                                    ]
                                ),
                                errors="coerce",
                            ).fillna(0).iloc[0]
                        ),
                    }
                )
            return pd.DataFrame(rows, columns=cols_summary), pd.DataFrame(columns=cols_detail)

    exon_path = _qual_exon_detail_table_path(qual_dir, project)
    if not os.path.isfile(exon_path):
        summary_rows = []
        for event_id in sorted(set(str(e) for e in event_ids)):
            structure = structure_by_event.get(event_id, {"roles": set(), "candidate": "NA"})
            roles = structure["roles"]
            candidate = structure["candidate"]
            summary_rows.append(
                {
                    "event_id": event_id,
                    "ID_position": _roles_to_string(roles),
                    "ID_position_summary": _roles_to_string(roles),
                    "ID_position_detail": "NA",
                    "ID_position_source": "transcript_structure" if roles else "NA",
                    "ID_position_roles": _roles_to_string(roles),
                    "ID_position_confidence": "structure_only" if roles else "NA",
                    "transcript_structure_roles": _roles_to_string(roles),
                    "position_evidence_relation": "transcript_structure_only" if roles else "NA",
                    "candidate_TE_event": candidate,
                    "candidate_TE_confidence": (
                        _format_candidate_confidence({c: "structure_only" for c in _candidate_role_map(candidate)})
                        if candidate != "NA"
                        else "NA"
                    ),
                    "ID_position_evaluable_sample_n": 0,
                    "ID_position_evaluable_replicate_n": 0,
                }
            )
        return pd.DataFrame(summary_rows, columns=cols_summary), pd.DataFrame(columns=cols_detail)

    df = pd.read_csv(exon_path, sep="\t")
    required = {"sample", "replicate", "exon", "strand", "ID_position"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise RuntimeError(f"HITindex exon table missing required column(s): {', '.join(missing)}")

    event_set = set(str(e) for e in event_ids)
    records = []
    warned_unknown = set()
    malformed_coord_n = 0
    for _, row in df.iterrows():
        try:
            chrom, start, end, strand = _parse_exon_coord(row.get("exon", ""), row.get("strand", ""))
        except ValueError:
            malformed_coord_n += 1
            continue
        event_id = _build_exon_event_id(chrom, start, end, strand)
        if event_id not in event_set:
            continue

        raw_position = str(row.get("ID_position", "")).strip()
        if raw_position and raw_position not in HITINDEX_POSITION_ROLES and raw_position not in warned_unknown:
            warned_unknown.add(raw_position)
            log_message(
                "[WARNING]",
                f"Unknown HITindex ID_position '{raw_position}' encountered in {exon_path}.",
                color="warning",
            )
        hit_roles, hit_strength = _parse_hitindex_position(raw_position)
        structure_roles = structure_by_event.get(event_id, {"roles": set()})["roles"]
        te_intervals = _parse_te_interval_tokens(
            row.get("te_splice_site_repeat_TE", ""),
            row.get("te_other_overlap_TE", ""),
        )
        candidate = _classify_candidate_te_events(
            hit_roles,
            chrom,
            start,
            end,
            strand,
            te_intervals,
        )
        records.append(
            {
                "exon": event_id,
                "sample": str(row.get("sample", "")),
                "replicate": str(row.get("replicate", "")),
                "gene_id": str(row.get("gene", "")),
                "transcript_id": _join_unique_values(str(row.get("matched_transcript_ids", "")).split(",")),
                "ID_position_raw": raw_position or "NA",
                "ID_position_roles_hitindex": _roles_to_string(hit_roles),
                "ID_position_strength": hit_strength,
                "ID_position_roles_transcript_structure": _roles_to_string(structure_roles),
                "ID_position_source": "hitindex" if hit_roles else "NA",
                "ID_position_confidence": "NA",
                "position_evidence_relation": _evidence_relation(hit_roles, structure_roles),
                "candidate_TE_event": candidate,
                "candidate_TE_confidence": "NA",
                "n_evaluable_replicates": pd.NA,
                "role_support_fraction": "NA",
                "HITindex": row.get("HITindex", pd.NA),
                "nFE": row.get("nFE", pd.NA),
                "nINTERNAL": row.get("nINTERNAL", pd.NA),
                "nLE": row.get("nLE", pd.NA),
                "nSINGLE": row.get("nSINGLE", pd.NA),
            }
        )
    if malformed_coord_n:
        log_message(
            "[WARNING]",
            f"Skipped {malformed_coord_n} malformed exon coordinate row(s) while reading HITindex exon table: {exon_path}",
            color="warning",
        )

    detail_df = pd.DataFrame(records, columns=cols_detail)
    if detail_df.empty:
        summary_rows = []
        for event_id in sorted(set(str(e) for e in event_ids)):
            structure = structure_by_event.get(event_id, {"roles": set(), "candidate": "NA"})
            roles = structure["roles"]
            candidate = structure["candidate"]
            summary_rows.append(
                {
                    "event_id": event_id,
                    "ID_position": _roles_to_string(roles),
                    "ID_position_summary": _roles_to_string(roles),
                    "ID_position_detail": "NA",
                    "ID_position_source": "transcript_structure" if roles else "NA",
                    "ID_position_roles": _roles_to_string(roles),
                    "ID_position_confidence": "structure_only" if roles else "NA",
                    "transcript_structure_roles": _roles_to_string(roles),
                    "position_evidence_relation": "transcript_structure_only" if roles else "NA",
                    "candidate_TE_event": candidate,
                    "candidate_TE_confidence": (
                        _format_candidate_confidence({c: "structure_only" for c in _candidate_role_map(candidate)})
                        if candidate != "NA"
                        else "NA"
                    ),
                    "ID_position_evaluable_sample_n": 0,
                    "ID_position_evaluable_replicate_n": 0,
                }
            )
        return pd.DataFrame(summary_rows, columns=cols_summary), detail_df

    sample_summary = []
    for (event_id, sample), group in detail_df.groupby(["exon", "sample"], sort=True):
        structure_roles = structure_by_event.get(str(event_id), {"roles": set()})["roles"]
        evaluable = group.loc[group["ID_position_roles_hitindex"].astype(str) != "NA"].copy()
        n_eval = int(evaluable["replicate"].astype(str).nunique())
        role_to_conf = {}
        role_support_fraction = {}
        for role in ["first", "internal", "last"]:
            role_rows = evaluable[
                evaluable["ID_position_roles_hitindex"].astype(str).apply(
                    lambda text, r=role: r in _parse_roles_text(text)
                )
            ]
            n_role = int(role_rows["replicate"].astype(str).nunique())
            if n_eval == 0 or n_role == 0:
                continue
            support_fraction = n_role / n_eval
            strengths = set(role_rows["ID_position_strength"].astype(str))
            has_high = "high" in strengths
            has_medium = "medium" in strengths
            structure_known = bool(structure_roles)
            consistent = (role in structure_roles) if structure_known else True
            if support_fraction >= 0.5 and has_high and consistent:
                conf = "high"
            elif (
                (support_fraction >= 0.5 and has_medium and consistent)
                or (support_fraction < 0.5 and has_high and consistent)
                or (support_fraction >= 0.5 and has_high and not consistent)
            ):
                conf = "medium"
            else:
                conf = "low"
            role_to_conf[role] = conf
            role_support_fraction[role] = round(support_fraction, 6)

        sample_roles = set(role_to_conf)
        sample_candidate = _choose_candidate_summary(group["candidate_TE_event"])
        candidate_to_conf = {}
        for candidate, role in _candidate_role_map(sample_candidate).items():
            if role in role_to_conf:
                candidate_to_conf[candidate] = role_to_conf[role]

        fraction_text = (
            ";".join(
                f"{role}:{role_support_fraction[role]:.6f}"
                for role in sorted(role_support_fraction, key=lambda x: POSITION_ROLE_ORDER.get(x, 99))
            )
            if role_support_fraction
            else "NA"
        )
        idx = detail_df.index[
            (detail_df["exon"].astype(str) == str(event_id))
            & (detail_df["sample"].astype(str) == str(sample))
        ]
        detail_df.loc[idx, "ID_position_confidence"] = _format_role_confidence(role_to_conf)
        detail_df.loc[idx, "candidate_TE_confidence"] = _format_candidate_confidence(candidate_to_conf)
        detail_df.loc[idx, "n_evaluable_replicates"] = n_eval
        detail_df.loc[idx, "role_support_fraction"] = fraction_text
        sample_summary.append(
            {
                "event_id": event_id,
                "sample": sample,
                "ID_position_roles": _roles_to_string(sample_roles),
                "ID_position_detail": _join_unique_values(group["ID_position_raw"]),
                "ID_position_confidence": _format_role_confidence(role_to_conf),
                "candidate_TE_event": sample_candidate,
                "candidate_TE_confidence": _format_candidate_confidence(candidate_to_conf),
                "n_evaluable_replicates": n_eval,
            }
        )
    sample_summary_df = pd.DataFrame(sample_summary)

    summary_rows = []
    all_event_ids = sorted(set(str(e) for e in event_ids))
    for event_id in all_event_ids:
        group = sample_summary_df.loc[sample_summary_df["event_id"].astype(str) == str(event_id)].copy()
        structure_roles = structure_by_event.get(str(event_id), {"roles": set(), "candidate": "NA"})["roles"]
        structure_candidate = structure_by_event.get(str(event_id), {"roles": set(), "candidate": "NA"})["candidate"]
        role_to_conf = {}
        for role in ["first", "internal", "last"]:
            role_conf = []
            for value, conf_text in zip(group.get("ID_position_roles", []), group.get("ID_position_confidence", [])):
                if role not in _parse_roles_text(value):
                    continue
                if ":" in str(conf_text):
                    for token in str(conf_text).split(";"):
                        key, _, val = token.partition(":")
                        if key == role:
                            role_conf.append(val)
                else:
                    role_conf.append(str(conf_text))
            best = _best_confidence(role_conf)
            if best != "NA":
                role_to_conf[role] = best

        hit_roles = set(role_to_conf)
        source = "hitindex" if hit_roles else ("transcript_structure" if structure_roles else "NA")
        final_roles = hit_roles if hit_roles else structure_roles
        final_confidence = _format_role_confidence(role_to_conf) if hit_roles else ("structure_only" if structure_roles else "NA")
        relation = _evidence_relation(hit_roles, structure_roles) if hit_roles else ("transcript_structure_only" if structure_roles else "NA")
        detail = _join_unique_values(detail_df.loc[detail_df["exon"] == event_id, "ID_position_raw"])
        candidate_to_conf = {}
        for cand_text, cand_conf_text in zip(group.get("candidate_TE_event", []), group.get("candidate_TE_confidence", [])):
            for candidate in str(cand_text).split(","):
                candidate = candidate.strip()
                if not candidate or candidate == "NA":
                    continue
                if ":" in str(cand_conf_text):
                    for token in str(cand_conf_text).split(";"):
                        key, _, val = token.partition(":")
                        if key == candidate:
                            candidate_to_conf.setdefault(candidate, []).append(val)
                else:
                    candidate_to_conf.setdefault(candidate, []).append(str(cand_conf_text))
        if not candidate_to_conf and not hit_roles and structure_candidate != "NA":
            for candidate in _candidate_role_map(structure_candidate):
                candidate_to_conf[candidate] = ["structure_only"]
        candidate_best = {
            candidate: _best_confidence(conf_values)
            for candidate, conf_values in candidate_to_conf.items()
        }
        candidate = _choose_candidate_summary(candidate_best.keys())
        candidate_confidence = _format_candidate_confidence(candidate_best)
        summary_rows.append(
            {
                "event_id": event_id,
                "ID_position": _roles_to_string(final_roles),
                "ID_position_summary": _roles_to_string(final_roles),
                "ID_position_detail": detail,
                "ID_position_source": source,
                "ID_position_roles": _roles_to_string(final_roles),
                "ID_position_confidence": final_confidence,
                "transcript_structure_roles": _roles_to_string(structure_roles),
                "position_evidence_relation": relation,
                "candidate_TE_event": candidate,
                "candidate_TE_confidence": candidate_confidence,
                "ID_position_evaluable_sample_n": int(
                    (pd.to_numeric(group.get("n_evaluable_replicates", pd.Series(dtype=float)), errors="coerce").fillna(0) > 0).sum()
                ) if not group.empty else 0,
                "ID_position_evaluable_replicate_n": int(
                    pd.to_numeric(group.get("n_evaluable_replicates", pd.Series(dtype=float)), errors="coerce").fillna(0).sum()
                ) if not group.empty else 0,
            }
        )

    return pd.DataFrame(summary_rows, columns=cols_summary), detail_df


def _load_exon_event_annotation_table(qual_dir, project, hitindex_summary_df=None):
    te_df = _load_te_exon_annotation_df(qual_dir, project)
    cols = [
        "metaexon",
        "event_id",
        "te_exon_id",
        "transcript_id",
        "gene",
        "gene_id",
        "chrom",
        "start",
        "end",
        "strand",
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
        "event_source",
        "te_overlap_label_raw",
        "te_overlap_label",
        "te_overlap_label_final",
        "te_overlap_n",
        "te_overlap_bp_max",
        "te_overlap_frac_max",
        "te_boundary_side",
        "te_boundary_hit_any",
        "junction_te_side_reads_max",
        "te_overlap_pass_raw",
        "te_overlap_pass_any",
        "te_overlap_pass_final",
        "te_splice_site_repeat_TE",
        "te_other_overlap_TE",
    ]
    if te_df.empty:
        return pd.DataFrame(columns=cols)

    work = te_df.copy()
    ambiguous_events = _find_ambiguous_multi_gene_events(work)
    ambiguous_event_ids = set(ambiguous_events)
    report_path = _write_ambiguous_multi_gene_report(qual_dir, project, ambiguous_events)
    _log_ambiguous_multi_gene_exclusion(ambiguous_events, report_path)

    def _col_or_default(col, default):
        if col in work.columns:
            return work[col]
        return pd.Series([default] * len(work), index=work.index)

    def _first_existing_col(candidates, default):
        for col in candidates:
            if col in work.columns:
                return work[col]
        return pd.Series([default] * len(work), index=work.index)

    work["event_id"] = _event_id_series(work)
    if ambiguous_event_ids:
        work = work.loc[~work["event_id"].astype(str).isin(ambiguous_event_ids)].copy()
        if work.empty:
            return pd.DataFrame(columns=cols)
    work["te_overlap_label_raw"] = _first_existing_col(
        ["te_overlap_label_raw", "te_overlap_label"],
        "TE_overlap",
    ).fillna("TE_overlap").astype(str)
    work["te_overlap_label_final"] = _first_existing_col(
        ["te_overlap_label_final", "te_overlap_label", "te_overlap_label_raw"],
        "TE_overlap",
    ).fillna("TE_overlap").astype(str)
    work["te_overlap_pass_raw"] = pd.to_numeric(
        _first_existing_col(["te_overlap_pass_raw", "te_overlap_pass_any"], 1),
        errors="coerce",
    ).fillna(1).astype(int)
    work["te_overlap_pass_final"] = pd.to_numeric(
        _first_existing_col(["te_overlap_pass_final", "te_overlap_pass_any", "te_overlap_pass_raw"], 1),
        errors="coerce",
    ).fillna(1).astype(int)
    work["te_overlap_n"] = pd.to_numeric(_col_or_default("te_overlap_n", 0), errors="coerce").fillna(0).astype(int)
    work["te_overlap_bp_max"] = pd.to_numeric(_col_or_default("te_overlap_bp_max", 0), errors="coerce").fillna(0).astype(int)
    work["te_overlap_frac_max"] = pd.to_numeric(_col_or_default("te_overlap_frac_max", 0.0), errors="coerce").fillna(0.0)
    work["te_boundary_side"] = _col_or_default("te_boundary_side", "none").fillna("none").astype(str)
    work["te_boundary_hit_any"] = pd.to_numeric(_col_or_default("te_boundary_hit_any", 0), errors="coerce").fillna(0).astype(int)
    work["junction_te_side_reads_max"] = pd.to_numeric(
        _col_or_default("junction_te_side_reads_max", 0.0),
        errors="coerce",
    ).fillna(0.0)
    work["te_splice_site_repeat_TE"] = _col_or_default("te_splice_site_repeat_TE", "").fillna("").astype(str)
    work["te_other_overlap_TE"] = _col_or_default("te_other_overlap_TE", "").fillna("").astype(str)

    hitindex_by_event = {}
    if hitindex_summary_df is not None and not hitindex_summary_df.empty:
        for _, row in hitindex_summary_df.iterrows():
            hitindex_by_event[str(row.get("event_id", ""))] = row.to_dict()

    grouped_rows = []
    for event_id, group in work.groupby("event_id", sort=True):
        genes = sorted({str(v).strip() for v in group["gene_id"] if str(v).strip()})
        hit = hitindex_by_event.get(str(event_id), {})
        position_summary = str(hit.get("ID_position_summary", "NA") or "NA")
        position_detail = str(hit.get("ID_position_detail", "NA") or "NA")
        position_source = str(hit.get("ID_position_source", "NA") or "NA")
        position_roles = str(hit.get("ID_position_roles", position_summary) or "NA")
        position_confidence = str(hit.get("ID_position_confidence", "NA") or "NA")
        structure_roles = str(hit.get("transcript_structure_roles", "NA") or "NA")
        evidence_relation = str(hit.get("position_evidence_relation", "NA") or "NA")
        candidate_te_event = str(hit.get("candidate_TE_event", "NA") or "NA")
        candidate_te_confidence = str(hit.get("candidate_TE_confidence", "NA") or "NA")
        grouped_rows.append(
            {
                "metaexon": str(event_id),
                "event_id": str(event_id),
                "te_exon_id": _join_unique_values(group["te_exon_id"]),
                "transcript_id": _join_unique_values(group["transcript_id"]),
                "gene": genes[0] if genes else "",
                "gene_id": genes[0] if genes else "",
                "chrom": str(group["exon_chrom"].iloc[0]),
                "start": int(pd.to_numeric(group["exon_start"], errors="coerce").iloc[0]),
                "end": int(pd.to_numeric(group["exon_end"], errors="coerce").iloc[0]),
                "strand": str(group["exon_strand"].iloc[0]),
                "ID_position": position_summary,
                "ID_position_summary": position_summary,
                "ID_position_detail": position_detail,
                "ID_position_source": position_source,
                "ID_position_roles": position_roles,
                "ID_position_confidence": position_confidence,
                "transcript_structure_roles": structure_roles,
                "position_evidence_relation": evidence_relation,
                "candidate_TE_event": candidate_te_event,
                "candidate_TE_confidence": candidate_te_confidence,
                "ID_position_evaluable_sample_n": int(hit.get("ID_position_evaluable_sample_n", 0) or 0),
                "ID_position_evaluable_replicate_n": int(hit.get("ID_position_evaluable_replicate_n", 0) or 0),
                "event_source": "te_exon",
                "te_overlap_label_raw": "TE_overlap" if int(group["te_overlap_pass_raw"].max()) == 1 else "no_overlap",
                "te_overlap_label": "TE_overlap" if int(group["te_overlap_pass_final"].max()) == 1 else "no_overlap",
                "te_overlap_label_final": "TE_overlap" if int(group["te_overlap_pass_final"].max()) == 1 else "no_overlap",
                "te_overlap_n": int(group["te_overlap_n"].max()),
                "te_overlap_bp_max": int(group["te_overlap_bp_max"].max()),
                "te_overlap_frac_max": float(group["te_overlap_frac_max"].max()),
                "te_boundary_side": _join_unique_values(
                    token
                    for value in group["te_boundary_side"]
                    for token in str(value).split(",")
                ) or "none",
                "te_boundary_hit_any": int(group["te_boundary_hit_any"].max()),
                "junction_te_side_reads_max": float(group["junction_te_side_reads_max"].max()),
                "te_overlap_pass_raw": int(group["te_overlap_pass_raw"].max()),
                "te_overlap_pass_any": int(group["te_overlap_pass_final"].max()),
                "te_overlap_pass_final": int(group["te_overlap_pass_final"].max()),
                "te_splice_site_repeat_TE": _join_unique_values(
                    token
                    for value in group["te_splice_site_repeat_TE"]
                    for token in str(value).split(",")
                ),
                "te_other_overlap_TE": _join_unique_values(
                    token
                    for value in group["te_other_overlap_TE"]
                    for token in str(value).split(",")
                ),
            }
        )
    return pd.DataFrame(grouped_rows, columns=cols).sort_values(by=["event_id"]).reset_index(drop=True)

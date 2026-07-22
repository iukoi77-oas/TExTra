"""Write TE-derived exon usage tables and clean quant outputs."""

import os
import shutil
from collections import defaultdict

import pandas as pd

from util.common.write_logs import log_message

EXPRESSED_ANY_MIN_USAGE = 0.1

DEFAULT_EXON_USAGE_COLUMNS = [
    "exon_id",
    "metaexon_id",
    "gene_id",
    "transcript_id",
    "te_overlap_label",
    "te_boundary_side",
    "te_splice_site_repeat_TE",
    "ID_position_summary",
    "ID_position_confidence",
    "candidate_TE_event",
    "candidate_TE_confidence",
]

DETAIL_EXON_USAGE_EXTRA_COLUMNS = [
    "te_overlap_bp_max",
    "te_overlap_frac_max",
    "junction_te_side_reads_max",
    "support_sample_n",
    "support_sample_ratio",
    "ID_position_source",
    "transcript_structure_roles",
    "HITindex_structure_roles",
    "HITindex_evaluable_sample_n",
    "HITindex_evaluable_replicate_n",
]


def _build_exon_event_usage_and_matrix(
    transcript_df,
    exon_events,
    event_to_txs,
    event_to_gene,
    tx_to_gene,
    event_anno_df,
):
    if transcript_df.empty:
        usage_cols = [
            "sample",
            "exon_id",
            "gene_id",
            "event_class",
            "te_class",
            "usage_tx",
            "supporting_transcripts",
            "competing_transcripts",
            "usage_definition",
        ]
        return pd.DataFrame(columns=usage_cols), pd.DataFrame({"exon": list(exon_events)})

    gene_to_txs = defaultdict(set)
    for tx_id, gene_id in tx_to_gene.items():
        gene_to_txs[str(gene_id)].add(str(tx_id))

    event_anno = {}
    if not event_anno_df.empty:
        for _, row in event_anno_df.iterrows():
            event_anno[str(row["metaexon"])] = row.to_dict()

    usage_rows = []
    usage_matrix_rows = []
    samples = sorted(transcript_df["sample"].astype(str).unique())
    for sample in samples:
        sub = transcript_df[transcript_df["sample"].astype(str) == sample].copy()
        sub["tpm"] = pd.to_numeric(sub["tpm"], errors="coerce").fillna(0.0)
        tx_val = sub.set_index("transcript_id")["tpm"].to_dict()
        gene_tot = sub.groupby("gene_id", as_index=False)["tpm"].sum().set_index("gene_id")["tpm"].to_dict()

        for event_id in exon_events:
            support = sorted(event_to_txs.get(event_id, set()))
            gene_id = str(event_to_gene.get(event_id, "") or "")
            if not gene_id and support:
                gset = sorted({str(tx_to_gene.get(tx, "")) for tx in support if str(tx_to_gene.get(tx, ""))})
                if len(gset) == 1:
                    gene_id = gset[0]

            num = float(sum(float(tx_val.get(tx, 0.0)) for tx in support))
            den = float(gene_tot.get(gene_id, 0.0)) if gene_id else 0.0
            usage = (num / den) if den > 0 else 0.0

            te_class = ""
            ann = event_anno.get(event_id, {})
            if ann:
                te_class = str(ann.get("te_splice_site_repeat_TE", "") or ann.get("te_other_overlap_TE", "") or "")

            competing = sorted(set(gene_to_txs.get(gene_id, set())) - set(support)) if gene_id else []
            event_class = str(ann.get("event_source", "exon") or "exon") if ann else "exon"
            usage_rows.append(
                {
                    "sample": sample,
                    "exon_id": event_id,
                    "gene_id": gene_id,
                    "event_class": event_class,
                    "te_class": te_class,
                    "usage_tx": usage,
                    "supporting_transcripts": ",".join(support),
                    "competing_transcripts": ",".join(competing),
                    "usage_definition": "sum_tpm_supporting_transcripts / sum_tpm_gene_transcripts",
                }
            )
            usage_matrix_rows.append({"exon": event_id, "sample": sample, "usage": usage})

    usage_df = pd.DataFrame(usage_rows)
    usage_long = pd.DataFrame(usage_matrix_rows)
    if usage_long.empty:
        usage_wide = pd.DataFrame({"exon": list(exon_events)})
    else:
        usage_wide = (
            usage_long.pivot_table(index="exon", columns="sample", values="usage", aggfunc="mean", fill_value=0.0)
            .reset_index()
        )
    return usage_df, usage_wide


def _log_missing_supporting_transcripts(transcript_df, event_to_txs):
    if transcript_df.empty or not event_to_txs:
        return
    expected_txs = set()
    for txs in event_to_txs.values():
        expected_txs.update(str(tx) for tx in txs if str(tx).strip())
    if not expected_txs:
        return

    for sample, sub in transcript_df.groupby("sample", sort=True):
        present_txs = set(sub["transcript_id"].astype(str))
        missing_txs = sorted(expected_txs - present_txs)
        if missing_txs:
            log_message(
                "[WARNING]",
                (
                    f"Sample '{sample}' is missing {len(missing_txs)} TE-supporting transcript(s) "
                    "in RSEM isoforms output; their usage contribution is treated as 0."
                ),
                color="warning",
            )


def _write_project_exon_usage_tables(
    usage_wide,
    sample_list,
    out_dir,
    project,
    event_anno_df,
    write_detail=False,
):
    # enrich with annotation columns for mode3 compatibility
    if event_anno_df.empty:
        anno = pd.DataFrame({"metaexon": usage_wide["exon"].astype(str)})
    else:
        anno = event_anno_df.copy()
    if "metaexon" not in anno.columns:
        anno["metaexon"] = usage_wide["exon"].astype(str)

    merged = usage_wide.merge(anno, left_on="exon", right_on="metaexon", how="left")
    merged["exon"] = merged["exon"].astype(str)
    if "gene_id" not in merged.columns:
        if "gene" in merged.columns:
            merged["gene_id"] = merged["gene"].fillna("").astype(str)
        else:
            merged["gene_id"] = ""

    defaults = {
        "ID_position": "NA",
        "ID_position_summary": "NA",
        "ID_position_detail": "NA",
        "ID_position_source": "NA",
        "ID_position_roles": "NA",
        "ID_position_confidence": "NA",
        "transcript_structure_roles": "NA",
        "position_evidence_relation": "NA",
        "candidate_TE_event": "NA",
        "candidate_TE_confidence": "NA",
        "ID_position_evaluable_sample_n": 0,
        "ID_position_evaluable_replicate_n": 0,
        "te_overlap_label": "no_overlap",
        "te_overlap_bp_max": 0,
        "te_overlap_frac_max": 0.0,
        "te_boundary_side": "none",
        "te_boundary_hit_any": 0,
        "junction_te_side_reads_max": 0.0,
        "te_overlap_pass_any": 0,
        "te_splice_site_repeat_TE": "",
        "te_other_overlap_TE": "",
    }
    for col, val in defaults.items():
        if col not in merged.columns:
            merged[col] = val
        merged[col] = merged[col].fillna(val)
    merged["ID_position"] = merged["ID_position"].astype(str)
    merged["ID_position_summary"] = merged["ID_position_summary"].astype(str)
    merged["ID_position_detail"] = merged["ID_position_detail"].astype(str)
    merged["ID_position_source"] = merged["ID_position_source"].astype(str)
    merged["ID_position_roles"] = merged["ID_position_roles"].astype(str)
    merged["ID_position_confidence"] = merged["ID_position_confidence"].astype(str)
    merged["transcript_structure_roles"] = merged["transcript_structure_roles"].astype(str)
    merged["position_evidence_relation"] = merged["position_evidence_relation"].astype(str)
    merged["candidate_TE_event"] = merged["candidate_TE_event"].astype(str)
    merged["candidate_TE_confidence"] = merged["candidate_TE_confidence"].astype(str)

    # ensure sample columns always exist
    for sample in sample_list:
        if sample not in merged.columns:
            merged[sample] = 0.0
        merged[sample] = pd.to_numeric(merged[sample], errors="coerce").fillna(0.0)

    if sample_list:
        merged["expressed_any"] = (
            merged[sample_list].max(axis=1) >= float(EXPRESSED_ANY_MIN_USAGE)
        ).astype(int)
        support_matrix = merged[sample_list].ge(float(EXPRESSED_ANY_MIN_USAGE))
        merged["support_sample_n"] = support_matrix.sum(axis=1).astype(int)
        merged["support_sample_ratio"] = (
            merged["support_sample_n"].astype(float) / float(len(sample_list))
        )
    else:
        merged["expressed_any"] = 0
        merged["support_sample_n"] = 0
        merged["support_sample_ratio"] = 0.0

    merged["exon_id"] = merged["exon"].astype(str)
    merged["metaexon_id"] = merged["metaexon"].fillna(merged["exon"]).astype(str)
    merged["HITindex_structure_roles"] = merged["ID_position_roles"].astype(str)
    merged["HITindex_evaluable_sample_n"] = pd.to_numeric(
        merged["ID_position_evaluable_sample_n"],
        errors="coerce",
    ).fillna(0).astype(int)
    merged["HITindex_evaluable_replicate_n"] = pd.to_numeric(
        merged["ID_position_evaluable_replicate_n"],
        errors="coerce",
    ).fillna(0).astype(int)

    sample_cols = [sample for sample in sample_list if sample in merged.columns]
    default_cols = ["exon_id"] + sample_cols + [c for c in DEFAULT_EXON_USAGE_COLUMNS if c != "exon_id"]
    for col in default_cols:
        if col not in merged.columns:
            merged[col] = ""
    default_df = merged[default_cols].copy()

    project_out = os.path.join(out_dir, f"{project}.TE_overlap.exon_usage.tsv")
    default_df.to_csv(project_out, sep="\t", index=False)

    detail_out = None
    if write_detail:
        detail_cols = default_cols + [c for c in DETAIL_EXON_USAGE_EXTRA_COLUMNS if c not in default_cols]
        for col in detail_cols:
            if col not in merged.columns:
                merged[col] = ""
        detail_out = os.path.join(out_dir, f"{project}.TE_overlap.exon_usage.detail.tsv")
        merged[detail_cols].to_csv(detail_out, sep="\t", index=False)

    return project_out, detail_out


def _cleanup_quant_outputs(count_dir, project, keep_gene_abundance=False, keep_detail=False, keep_dirs=None):
    keep_dirs = set(keep_dirs or [])
    keep_exact = {
        f"{project}.TE_overlap.exon_usage.tsv",
        "transcript_abundance.tsv",
    }
    if keep_detail:
        keep_exact.add(f"{project}.TE_overlap.exon_usage.detail.tsv")
    if keep_gene_abundance:
        keep_exact.add("gene_abundance.tsv")

    for name in os.listdir(count_dir):
        path = os.path.join(count_dir, name)
        if os.path.isdir(path):
            if name in keep_dirs:
                continue
            shutil.rmtree(path, ignore_errors=True)
            continue
        if name in keep_exact:
            continue
        try:
            os.remove(path)
        except FileNotFoundError:
            pass

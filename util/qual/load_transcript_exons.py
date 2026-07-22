"""Load consensus transcript-exon structures for qual mode."""

import collections
import os

import pandas as pd
from pandas.errors import EmptyDataError, ParserError


def _make_transcript_exon_id(row):
    exon_id = str(row.get("exon_id", "")).strip()
    transcript_id = str(row.get("transcript_id", "")).strip()
    if exon_id:
        return f"{transcript_id}:{exon_id}" if transcript_id else exon_id
    chrom = str(row.get("exon_chrom", row.get("chrom", "")))
    start = str(row.get("exon_start", row.get("start0", "")))
    end = str(row.get("exon_end", row.get("end0", "")))
    strand = str(row.get("exon_strand", row.get("strand", "")))
    return f"{transcript_id}:{chrom}:{start}-{end}:{strand}" if transcript_id else f"{chrom}:{start}-{end}:{strand}"


def _load_transcript_assignment_rows(assignment_tsv_path):
    if not assignment_tsv_path or not os.path.isfile(assignment_tsv_path):
        return {}
    try:
        df = pd.read_csv(assignment_tsv_path, sep="\t")
    except (OSError, EmptyDataError, ParserError, UnicodeDecodeError) as exc:
        raise RuntimeError(f"Failed to read transcript-gene assignment table: {assignment_tsv_path}") from exc
    if "transcript_id" not in df.columns:
        raise RuntimeError(f"Invalid transcript-gene assignment table: {assignment_tsv_path} missing transcript_id column.")
    df = df.copy()
    df["transcript_id"] = df["transcript_id"].fillna("").astype(str).str.strip()
    df = df.loc[df["transcript_id"] != ""].copy()
    return {
        str(row["transcript_id"]): {str(k): row[k] for k in df.columns}
        for _, row in df.iterrows()
    }


def _load_mapping_transcript_exon_rows(mapping_tsv_path):
    if not mapping_tsv_path or not os.path.isfile(mapping_tsv_path):
        return []
    try:
        df = pd.read_csv(mapping_tsv_path, sep="\t")
    except (OSError, EmptyDataError, ParserError, UnicodeDecodeError):
        return []

    required_cols = ["transcript_id", "gene_id", "exon_chrom", "exon_start", "exon_end", "exon_strand"]
    if not set(required_cols).issubset(df.columns):
        return []

    select_cols = list(required_cols)
    if "metaexon_id" in df.columns:
        select_cols.append("metaexon_id")
    valid = df[select_cols].copy()
    valid["transcript_id"] = valid["transcript_id"].fillna("").astype(str).str.strip()
    valid = valid.loc[valid["transcript_id"] != ""].copy()
    if valid.empty:
        return []

    valid["gene_id"] = valid["gene_id"].fillna("").astype(str).str.strip()
    valid["exon_chrom"] = valid["exon_chrom"].fillna("").astype(str).str.strip()
    valid["exon_strand"] = valid["exon_strand"].fillna("").astype(str).str.strip()
    valid["exon_start"] = pd.to_numeric(valid["exon_start"], errors="coerce")
    valid["exon_end"] = pd.to_numeric(valid["exon_end"], errors="coerce")
    valid = valid.dropna(subset=["exon_start", "exon_end"]).copy()
    valid["exon_start"] = valid["exon_start"].astype(int)
    valid["exon_end"] = valid["exon_end"].astype(int)
    valid = valid.loc[(valid["exon_chrom"] != "") & (valid["exon_end"] > valid["exon_start"])].copy()
    if valid.empty:
        return []

    dedup_cols = ["transcript_id", "gene_id", "exon_chrom", "exon_start", "exon_end", "exon_strand"]
    if "metaexon_id" in valid.columns:
        dedup_cols = ["metaexon_id"] + dedup_cols
    valid = valid.drop_duplicates(subset=dedup_cols)
    return valid.to_dict("records")


def load_consensus_transcript_exon_rows(transcript_gtf_path, multiexon_only=True, assignment_tsv_path=None):
    """Load consensus transcript exons and annotate exon order/count per transcript."""
    assignment_by_tx = _load_transcript_assignment_rows(assignment_tsv_path)
    raw_rows = []
    with open(transcript_gtf_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            try:
                start0 = int(parts[3]) - 1
                end0 = int(parts[4])
            except ValueError:
                continue
            if end0 <= start0:
                continue
            attrs = {}
            for token in str(parts[8]).split(";"):
                token = token.strip()
                if not token or " " not in token:
                    continue
                key, value = token.split(" ", 1)
                attrs[key] = value.strip().strip('"')
            transcript_id = str(attrs.get("transcript_id", "")).strip()
            if not transcript_id:
                continue
            assignment = assignment_by_tx.get(transcript_id, {})
            gene_id = str(attrs.get("gene_id", "")).strip()
            exon_id = str(attrs.get("exon_id") or attrs.get("exon_number") or "").strip()
            exon_number_raw = str(attrs.get("exon_number", "")).strip()
            transcript_strand = str(parts[6])
            reference_gene_strand = str(
                assignment.get("reference_strand", attrs.get("reference_strand", "")) or ""
            ).strip()
            if reference_gene_strand not in {"+", "-"}:
                reference_gene_strand = ""
            gene_strand = reference_gene_strand or transcript_strand
            strand_consistent = str(
                assignment.get("strand_consistent", attrs.get("strand_consistent", "unknown")) or "unknown"
            ).strip()
            raw_rows.append(
                {
                    "gene_id": gene_id,
                    "transcript_id": transcript_id,
                    "exon_id": exon_id,
                    "exon_number_raw": exon_number_raw,
                    "exon_chrom": str(parts[0]),
                    "exon_start": int(start0),
                    "exon_end": int(end0),
                    "exon_strand": str(parts[6]),
                    "gene_strand": gene_strand,
                    "transcript_strand": transcript_strand,
                    "reference_gene_strand": reference_gene_strand,
                    "gene_assignment_strand_consistent": strand_consistent,
                    "gene_assignment_status": str(
                        assignment.get("assignment_status", attrs.get("gene_assignment_status", "")) or ""
                    ),
                }
            )

    if not raw_rows:
        return []

    by_tx = collections.defaultdict(list)
    for row in raw_rows:
        by_tx[row["transcript_id"]].append(row)

    out_rows = []
    for transcript_id, rows in by_tx.items():
        tx_count = len(rows)
        if multiexon_only and tx_count < 2:
            continue
        strand = str(rows[0].get("exon_strand", ""))
        reverse = strand == "-"
        sorted_rows = sorted(rows, key=lambda r: (int(r["exon_start"]), int(r["exon_end"])), reverse=reverse)
        for idx, row in enumerate(sorted_rows, start=1):
            next_row = dict(row)
            try:
                exon_number = int(str(row.get("exon_number_raw", "")).strip())
            except ValueError:
                exon_number = int(idx)
            if not next_row.get("exon_id"):
                next_row["exon_id"] = f"exon_{exon_number}"
            next_row["exon_number"] = int(exon_number)
            next_row["transcript_exon_count"] = int(tx_count)
            next_row["is_multiexon_transcript"] = int(tx_count >= 2)
            next_row["transcript_exon_id"] = _make_transcript_exon_id(next_row)
            reference_gene_strand = str(next_row.get("reference_gene_strand", ""))
            next_row["gene_exon_strand_consistent"] = int(
                reference_gene_strand in {"+", "-"}
                and reference_gene_strand == str(next_row.get("exon_strand", ""))
            )
            if reference_gene_strand in {"+", "-"}:
                next_row["te_gene_strand_relation"] = (
                    "same" if reference_gene_strand == str(next_row.get("exon_strand", "")) else "opposite"
                )
            else:
                next_row["te_gene_strand_relation"] = "unknown"
            out_rows.append(next_row)
    return out_rows

"""Manage qual staging, resume, support, and published outputs."""

import os
import shutil
import tempfile
from collections import Counter, defaultdict

import pandas as pd

from util.common.define_layout import CLASSIFY_DIR
from util.common.write_logs import log_message


def publish_classify_dir(staged_classify_dir, final_classify_dir):
    """Replace the final classify directory only after a staged run succeeds."""
    final_parent = os.path.dirname(final_classify_dir)
    os.makedirs(final_parent, exist_ok=True)
    backup_dir = None
    if os.path.exists(final_classify_dir):
        backup_dir = tempfile.mkdtemp(prefix=f".{CLASSIFY_DIR}.backup.", dir=final_parent)
        os.rmdir(backup_dir)
        os.rename(final_classify_dir, backup_dir)
    try:
        os.rename(staged_classify_dir, final_classify_dir)
    except OSError:
        if backup_dir is not None and os.path.exists(backup_dir) and not os.path.exists(final_classify_dir):
            os.rename(backup_dir, final_classify_dir)
        raise
    if backup_dir is not None:
        shutil.rmtree(backup_dir, ignore_errors=True)


def cleanup_hitindex_dir(classify_dir, calculate_afe_ale=False):
    """Keep only per-sample HITindex outputs requested by current qual switches."""
    hitindex_dir = os.path.join(classify_dir, "HITindex")
    if not os.path.isdir(hitindex_dir):
        return
    keep_suffixes = [".exon"]
    if calculate_afe_ale:
        keep_suffixes.extend([".AFEPSI", ".ALEPSI"])
    keep_suffixes = tuple(keep_suffixes)
    for name in os.listdir(hitindex_dir):
        path = os.path.join(hitindex_dir, name)
        if os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=True)
            continue
        if name.endswith(keep_suffixes):
            continue
        try:
            os.remove(path)
        except FileNotFoundError:
            pass


def prune_debug_hitindex_outputs(classify_dir, project, calculate_afe_ale=False, te_overlap_junction_evidence=False):
    """Keep only debug-approved reusable HITindex outputs."""
    hitindex_dir = os.path.join(classify_dir, "HITindex")
    if not os.path.isdir(hitindex_dir):
        return
    support_name = f"{project}.te_overlap_transcript_exon_junction_support.tsv"
    keep_suffixes = [".exon"]
    if calculate_afe_ale:
        keep_suffixes.extend([".AFEPSI", ".ALEPSI"])
    keep_suffixes = tuple(keep_suffixes)
    for name in os.listdir(hitindex_dir):
        path = os.path.join(hitindex_dir, name)
        if os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=True)
            continue
        keep = name.endswith(keep_suffixes)
        if te_overlap_junction_evidence and name == support_name:
            keep = True
        if keep:
            continue
        try:
            os.remove(path)
        except FileNotFoundError:
            pass


def is_complete_hitindex_replicate(hitindex_dir, replicate, calculate_afe_ale=False):
    """Return True when a replicate has all required reusable HITindex outputs."""
    exon_path = os.path.join(hitindex_dir, f"{replicate}.exon")
    if not os.path.isfile(exon_path):
        return False
    try:
        exon_df = pd.read_csv(exon_path, sep="\t", nrows=1)
    except (OSError, ValueError, pd.errors.EmptyDataError, pd.errors.ParserError):
        return False
    required_cols = {"exon", "gene", "strand", "HITindex", "ID", "ID_position"}
    if not required_cols.issubset(exon_df.columns):
        return False
    if calculate_afe_ale:
        for suffix in ["AFEPSI", "ALEPSI"]:
            if not os.path.isfile(os.path.join(hitindex_dir, f"{replicate}.{suffix}")):
                return False
    return True


def inspect_junction_support_reuse(source_hitindex_dir, project):
    """Return whether reusable project-level junction support is available."""
    support_path = os.path.join(
        source_hitindex_dir,
        f"{project}.te_overlap_transcript_exon_junction_support.tsv",
    )
    if not os.path.isfile(support_path):
        return False, f"missing {os.path.basename(support_path)}"
    try:
        support_df = pd.read_csv(support_path, sep="\t")
    except (OSError, ValueError, pd.errors.EmptyDataError, pd.errors.ParserError) as exc:
        return False, f"failed to read {os.path.basename(support_path)}: {exc}"
    required_cols = {"sample", "exon", "strand", "nleft", "nright"}
    missing = sorted(required_cols - set(support_df.columns))
    if missing:
        return False, f"{os.path.basename(support_path)} missing column(s): {', '.join(missing)}"
    if support_df.empty:
        return False, f"{os.path.basename(support_path)} is empty"
    return True, support_path


def inspect_afe_ale_reuse(source_hitindex_dir, replicates):
    """Return AFE/ALE reuse status across all expected replicates."""
    required_cols = {
        "AFEPSI": {"gene", "exon", "strand", "AFEPSI"},
        "ALEPSI": {"gene", "exon", "strand", "ALEPSI"},
    }
    files_seen = 0
    problems = []
    for replicate in replicates:
        for suffix, required in required_cols.items():
            path = os.path.join(source_hitindex_dir, f"{replicate}.{suffix}")
            if not os.path.isfile(path):
                problems.append(f"{replicate}.{suffix}: missing")
                continue
            files_seen += 1
            try:
                df = pd.read_csv(path, sep="\t", nrows=1)
            except (OSError, ValueError, pd.errors.EmptyDataError, pd.errors.ParserError) as exc:
                problems.append(f"{replicate}.{suffix}: failed to read ({exc})")
                continue
            missing = sorted(required - set(df.columns))
            if missing:
                problems.append(f"{replicate}.{suffix}: missing column(s) {', '.join(missing)}")
            elif df.empty:
                problems.append(f"{replicate}.{suffix}: empty")
    if files_seen == 0:
        return "absent", "no per-replicate AFE/ALE files were found"
    if problems:
        preview = "; ".join(problems[:4])
        if len(problems) > 4:
            preview += f"; ... ({len(problems)} issue(s) total)"
        return "invalid", preview
    return "complete", f"replicates={len(replicates)}"


def seed_hitindex_outputs(
    source_hitindex_dir,
    staged_classify_dir,
    project,
    replicates,
    calculate_afe_ale=False,
    te_overlap_junction_evidence=False,
):
    """Copy completed HITindex outputs from an explicit source directory."""
    if not os.path.isdir(source_hitindex_dir):
        return 0
    target_hitindex_dir = os.path.join(staged_classify_dir, "HITindex")
    os.makedirs(target_hitindex_dir, exist_ok=True)
    copied = 0
    for replicate in replicates:
        if not is_complete_hitindex_replicate(
            source_hitindex_dir,
            replicate,
            calculate_afe_ale=calculate_afe_ale,
        ):
            continue
        for suffix in [".exon", ".AFEPSI", ".ALEPSI"]:
            if suffix in {".AFEPSI", ".ALEPSI"} and not calculate_afe_ale:
                continue
            src = os.path.join(source_hitindex_dir, f"{replicate}{suffix}")
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(target_hitindex_dir, os.path.basename(src)))
        copied += 1
    if te_overlap_junction_evidence:
        support_name = f"{project}.te_overlap_transcript_exon_junction_support.tsv"
        src_support = os.path.join(source_hitindex_dir, support_name)
        if os.path.isfile(src_support):
            shutil.copy2(src_support, os.path.join(target_hitindex_dir, support_name))
    return copied


def finalize_te_overlap_support_table(classify_dir, project):
    """Sort and report the TE-overlap transcript-exon junction support table."""
    support_path = os.path.join(
        classify_dir,
        "HITindex",
        f"{project}.te_overlap_transcript_exon_junction_support.tsv",
    )
    if not os.path.isfile(support_path):
        return
    try:
        df = pd.read_csv(support_path, sep="\t")
    except (OSError, ValueError, pd.errors.EmptyDataError, pd.errors.ParserError):
        return
    required_cols = {"exon", "strand", "nleft", "nright"}
    if not required_cols.issubset(df.columns):
        return

    work = df.copy()
    work["exon"] = work["exon"].astype(str)
    work["strand"] = work["strand"].astype(str)
    work["nleft"] = pd.to_numeric(work["nleft"], errors="coerce").fillna(0).astype(int)
    work["nright"] = pd.to_numeric(work["nright"], errors="coerce").fillna(0).astype(int)
    work["__event_key"] = work["exon"] + "|" + work["strand"]

    event_total = int(work["__event_key"].nunique())
    nonzero_events = int(
        work.loc[(work["nleft"] > 0) | (work["nright"] > 0), "__event_key"].nunique()
    )
    zero_events = int(max(0, event_total - nonzero_events))
    if zero_events > 0:
        log_message(
            "[INFO]",
            (
                "TE-overlap exon junction support summary: "
                f"retained {event_total} exon(s), including {zero_events} exon(s) "
                "with nleft=0 and nright=0 in all samples."
            ),
            color="info",
        )

    full = work.drop(columns=["__event_key"], errors="ignore")
    if not full.empty:
        sort_cols = [c for c in ["sample", "gene_id", "transcript_id", "chrom", "start", "end"] if c in full.columns]
        if sort_cols:
            full = full.sort_values(by=sort_cols, ascending=[True] * len(sort_cols))
    keep_cols = [
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
        "te_overlap_label",
        "te_overlap_bp_max",
        "te_overlap_frac_max",
        "te_boundary_hit_any",
        "te_splice_site_repeat_TE",
    ]
    full = full[[col for col in keep_cols if col in full.columns]]
    full.to_csv(support_path, sep="\t", index=False)


def load_te_overlap_support_records(classify_dir, project, require_records=False):
    """Load TE-overlap support records produced by per-replicate HITindex runs."""
    support_path = os.path.join(
        classify_dir,
        "HITindex",
        f"{project}.te_overlap_transcript_exon_junction_support.tsv",
    )
    if not os.path.isfile(support_path):
        if require_records:
            raise RuntimeError(
                f"TE-overlap junction evidence is enabled, but support table is missing: {support_path}"
            )
        return []
    try:
        support_df = pd.read_csv(support_path, sep="\t")
    except (OSError, ValueError, pd.errors.EmptyDataError, pd.errors.ParserError) as exc:
        raise RuntimeError(f"Failed to read TE-overlap junction support table: {support_path}") from exc
    required_cols = {"sample", "exon", "strand", "nleft", "nright"}
    missing = sorted(required_cols - set(support_df.columns))
    if missing:
        raise RuntimeError(
            f"Invalid TE-overlap junction support table: {support_path} missing column(s): {', '.join(missing)}"
        )
    if support_df.empty and require_records:
        raise RuntimeError(
            f"TE-overlap junction evidence is enabled, but support table is empty: {support_path}"
        )
    return support_df.to_dict("records")


def load_metaexon_rows(metaexon_bed_path):
    """Load metaexon BED rows for transcript-exon overlap mapping."""
    rows = []
    with open(metaexon_bed_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom, start_s, end_s, name, _score, strand = parts[:6]
            try:
                start0 = int(start_s)
                end0 = int(end_s)
            except ValueError:
                continue
            if end0 <= start0:
                continue
            name_parts = str(name).split(";")
            gene_id = name_parts[1] if len(name_parts) > 1 else ""
            rows.append(
                {
                    "metaexon_id": f"{chrom}:{start0}-{end0}:{strand}",
                    "gene_id": str(gene_id),
                    "chrom": str(chrom),
                    "start0": int(start0),
                    "end0": int(end0),
                    "strand": str(strand),
                }
            )
    return rows


def build_exon_index(exon_rows, bin_size=50000):
    """Build a simple genomic bin index for transcript exon rows."""
    index = defaultdict(list)
    for row in exon_rows:
        b0 = int(row["start0"] // bin_size)
        b1 = int((row["end0"] - 1) // bin_size)
        for b in range(b0, b1 + 1):
            index[(row["chrom"], row["strand"], b)].append(row)
    return index


def export_metaexon_exon_transcript_map(metaexon_bed_path, transcript_gtf_path, out_tsv_path):
    """Export overlap relationships between metaexons and consensus transcript exons."""
    from util.qual.load_transcript_exons import load_consensus_transcript_exon_rows

    meta_rows = load_metaexon_rows(metaexon_bed_path)
    exon_rows = [
        {
            "transcript_id": str(row.get("transcript_id", "")),
            "gene_id": str(row.get("gene_id", "")),
            "exon_id": str(row.get("exon_id", "")),
            "chrom": str(row.get("exon_chrom", "")),
            "start0": int(row.get("exon_start", -1)),
            "end0": int(row.get("exon_end", -1)),
            "strand": str(row.get("exon_strand", "")),
        }
        for row in load_consensus_transcript_exon_rows(transcript_gtf_path, multiexon_only=False)
    ]
    exon_index = build_exon_index(exon_rows)

    out_rows = []
    for meta in meta_rows:
        b0 = int(meta["start0"] // 50000)
        b1 = int((meta["end0"] - 1) // 50000)
        matched = 0
        seen = set()
        for b in range(b0, b1 + 1):
            cands = exon_index.get((meta["chrom"], meta["strand"], b), [])
            for exon in cands:
                if meta["gene_id"] and exon["gene_id"] and (meta["gene_id"] != exon["gene_id"]):
                    continue
                ov_bp = min(meta["end0"], exon["end0"]) - max(meta["start0"], exon["start0"])
                if ov_bp <= 0:
                    continue
                key = (
                    meta["metaexon_id"],
                    exon["transcript_id"],
                    exon["exon_id"],
                    exon["chrom"],
                    exon["start0"],
                    exon["end0"],
                    exon["strand"],
                )
                if key in seen:
                    continue
                seen.add(key)
                matched += 1
                out_rows.append(
                    {
                        "metaexon_id": meta["metaexon_id"],
                        "gene_id": meta["gene_id"],
                        "metaexon_chrom": meta["chrom"],
                        "metaexon_start": meta["start0"],
                        "metaexon_end": meta["end0"],
                        "metaexon_strand": meta["strand"],
                        "transcript_id": exon["transcript_id"],
                        "exon_id": exon["exon_id"],
                        "exon_chrom": exon["chrom"],
                        "exon_start": exon["start0"],
                        "exon_end": exon["end0"],
                        "exon_strand": exon["strand"],
                        "overlap_bp": int(ov_bp),
                    }
                )
        if matched == 0:
            out_rows.append(
                {
                    "metaexon_id": meta["metaexon_id"],
                    "gene_id": meta["gene_id"],
                    "metaexon_chrom": meta["chrom"],
                    "metaexon_start": meta["start0"],
                    "metaexon_end": meta["end0"],
                    "metaexon_strand": meta["strand"],
                    "transcript_id": "",
                    "exon_id": "",
                    "exon_chrom": "",
                    "exon_start": "",
                    "exon_end": "",
                    "exon_strand": "",
                    "overlap_bp": 0,
                }
            )

    out_df = pd.DataFrame(
        out_rows,
        columns=[
            "metaexon_id",
            "gene_id",
            "metaexon_chrom",
            "metaexon_start",
            "metaexon_end",
            "metaexon_strand",
            "transcript_id",
            "exon_id",
            "exon_chrom",
            "exon_start",
            "exon_end",
            "exon_strand",
            "overlap_bp",
        ],
    )
    out_df = out_df.sort_values(
        by=["metaexon_id", "transcript_id", "exon_start", "exon_end"],
        ascending=[True, True, True, True],
    )
    out_df.to_csv(out_tsv_path, sep="\t", index=False)
    return out_tsv_path, int(len(meta_rows)), int((out_df["transcript_id"] != "").sum())

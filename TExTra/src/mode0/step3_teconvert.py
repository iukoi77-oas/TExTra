"""Annotation conversion helpers for TE and consensus transcript BED files."""

import os
import pandas as pd
import shutil

from util.common.write_logs import log_message
from util.prep.convert_annotations import gtf_to_bed, merge_and_sort_bed
from util.common.define_layout import CONVERT_DIR


def convert_func(gtf, args):
    """Convert TE annotation and assembled gene GTF into BED files."""
    log_message("[INFO]", "Step 3/3: Annotation conversion", bold=True, color="step")

    convert_dir = os.path.join(args.out_dir, CONVERT_DIR)
    if os.path.exists(convert_dir):
        shutil.rmtree(convert_dir)
    os.makedirs(convert_dir)

    log_message("[INFO]", "Converting TE annotation.", color="info")
    log_message("[INFO]", f"TE annotation input: {args.te}", color="info", console=getattr(args, "debug", False))
    file_extension = os.path.splitext(args.te)[1].lower()
    te_bed_path = os.path.join(convert_dir, "TE_anno_1.bed")
    if file_extension in {".gtf", ".gff"}:
        gtf_to_bed(args.te, te_bed_path, 'TE')
    elif file_extension in {".out", ".txt", ".bed"}:
        rmsk_to_bed(args.te, te_bed_path, file_extension)
    else:
        log_message("[ERROR]", f"Please check te annotation file format (gtf/gff/out/txt/bed).", color="error")
        raise RuntimeError("Unsupported TE annotation file format. Expected gtf/gff/out/txt/bed.")
    
    extend_bed_path = None
    if args.extend:
        log_message("[INFO]", "Converting extended TE annotation.", color="info")
        log_message("[INFO]", f"Extended TE annotation input: {args.extend}", color="info", console=getattr(args, "debug", False))
        extend_extension = os.path.splitext(args.extend)[1].lower()
        extend_bed_path = os.path.join(convert_dir, "TE_anno_2.bed")

        if extend_extension in {".gtf", ".gff"}:
            gtf_to_bed(args.extend, extend_bed_path, 'TE')
        elif extend_extension in {".out", ".txt", ".bed"}:
            rmsk_to_bed(args.extend, extend_bed_path, extend_extension)
        else:
            log_message("[ERROR]", f"Please check extended TE annotation file format (gtf/gff/out/txt/bed).", color="error")
            raise RuntimeError("Unsupported extended TE annotation file format. Expected gtf/gff/out/txt/bed.")

    final_bed_path = os.path.join(convert_dir, "TE_anno.bed")
    merge_and_sort_bed(te_bed_path, extend_bed_path, final_bed_path)
    os.remove(te_bed_path)
    if extend_bed_path and os.path.exists(extend_bed_path):
        os.remove(extend_bed_path)

    log_message("[SUCCESS]", f"TE BED: {os.path.abspath(final_bed_path)}", color="success")

    log_message("[INFO]", "Converting consensus transcript annotation.", color="info")
    log_message("[INFO]", f"Consensus transcript annotation input: {gtf}", color="info", console=getattr(args, "debug", False))
    file_extension = os.path.splitext(gtf)[1].lower()
    gene_bed_path = os.path.join(convert_dir, "gene_anno.bed")
    if file_extension not in {".gtf", ".gff"}:
        raise RuntimeError("Unsupported consensus transcript annotation format. Expected gtf/gff.")
    gtf_to_bed(gtf, gene_bed_path, 'gene')
    log_message("[SUCCESS]", f"Gene BED: {os.path.abspath(gene_bed_path)}", color="success")


def rmsk_to_bed(rmsk_path, bed_output_path, file_extension):
    """Convert RepeatMasker-like formats (.txt/.out/.bed) to normalized BED6."""
    skip_classes = {"simple_repeat", "low_complexity", "satellite"}
    dropped_missing_cf = 0
    dropped_simple = 0

    def ensure_nonempty(df):
        if df.empty:
            raise RuntimeError(
                f"Converted TE BED is empty after parsing/filtering {rmsk_path}. "
                "Please check TE annotation format and filtering of Simple_repeat/Low_complexity/Satellite."
            )

    def split_two_columns(series, sep):
        split_cols = series.astype(str).str.split(sep, n=1, expand=True)
        if split_cols.shape[1] < 2:
            split_cols[1] = None
        return split_cols.iloc[:, :2]

    if file_extension == ".txt":
        df = pd.read_csv(rmsk_path, sep="\t", header=None, comment="#")

        df_bed = df.iloc[:, [5, 6, 7, 10, 1, 9, 11, 12]].copy()
        df_bed.columns = ["chrom", "start", "end", "te_name", "score", "strand", "class", "family"]
        df_bed["te_name"] = df_bed["te_name"].fillna("Unknown")

        before = len(df_bed)
        df_bed = df_bed.dropna(subset=["class", "family"])
        dropped_missing_cf += before - len(df_bed)
        before = len(df_bed)
        df_bed = df_bed[~df_bed["class"].astype(str).str.lower().isin(skip_classes)]
        dropped_simple += before - len(df_bed)
        df_bed["formatted"] = df_bed["te_name"] + ":" + df_bed["family"] + ":" + df_bed["class"]

    elif file_extension == ".out":
        records = []
        with open(rmsk_path) as handle:
            for line in handle:
                fields = line.strip().split()
                if len(fields) < 11:
                    continue
                try:
                    float(fields[0])
                except ValueError:
                    continue
                strand = "-" if fields[8] == "C" else fields[8]
                records.append(
                    {
                        "chrom": fields[4],
                        "start": fields[5],
                        "end": fields[6],
                        "te_name": fields[9],
                        "score": fields[0],
                        "strand": strand,
                        "classification": fields[10],
                    }
                )
        df_bed = pd.DataFrame.from_records(
            records,
            columns=["chrom", "start", "end", "te_name", "score", "strand", "classification"],
        )
        df_bed["start"] = pd.to_numeric(df_bed["start"], errors="coerce") - 1
        df_bed["end"] = pd.to_numeric(df_bed["end"], errors="coerce")
        df_bed = df_bed.dropna(subset=["start", "end"])
        df_bed["start"] = df_bed["start"].astype(int)
        df_bed["end"] = df_bed["end"].astype(int)
        ensure_nonempty(df_bed)
        df_bed["te_name"] = df_bed["te_name"].fillna("Unknown")
        df_bed[["class", "family"]] = split_two_columns(df_bed["classification"], "/")
        before = len(df_bed)
        df_bed = df_bed.dropna(subset=["class", "family"])
        dropped_missing_cf += before - len(df_bed)
        before = len(df_bed)
        df_bed = df_bed[~df_bed["class"].astype(str).str.lower().isin(skip_classes)]
        dropped_simple += before - len(df_bed)
        df_bed["formatted"] = df_bed["te_name"] + ":" + df_bed["family"] + ":" + df_bed["class"]
    
    elif file_extension == ".bed":
        df_bed = pd.read_csv(rmsk_path, sep="\t", header=None)
        if df_bed.shape[1] < 6:
            raise RuntimeError(f"BED TE annotation requires at least 6 columns: {rmsk_path}")
        df_bed = df_bed.iloc[:, :6].copy()
        df_bed.columns = ["chrom", "start", "end", "raw_formatted", "score", "strand"]

        def normalize_te_format(formatted):
            """transform 'SINE/Alu|AluSx' to 'AluSx:Alu:SINE' """
            if ":" in formatted:
                parts = formatted.split(":")
                if len(parts) >= 3 and parts[1] and parts[2]:
                    return f"{parts[0]}:{parts[1]}:{parts[2]}"
            if "|" in formatted and "/" in formatted:
                te_class, rest = formatted.split("/", 1)
                te_family, te_name = rest.split("|", 1)
                return f"{te_name}:{te_family}:{te_class}"
            return None

        df_bed["formatted"] = df_bed["raw_formatted"].apply(normalize_te_format)
        before = len(df_bed)
        df_bed = df_bed.dropna(subset=["formatted"])
        dropped_missing_cf += before - len(df_bed)
        split_cols = df_bed["formatted"].astype(str).str.split(":", n=2, expand=True)
        while split_cols.shape[1] < 3:
            split_cols[split_cols.shape[1]] = None
        df_bed["te_name"] = split_cols[0]
        df_bed["te_name"] = df_bed["te_name"].fillna("Unknown")
        df_bed["family"] = split_cols[1]
        df_bed["class"] = split_cols[2]
        before = len(df_bed)
        df_bed = df_bed[~df_bed["class"].astype(str).str.lower().isin(skip_classes)]
        dropped_simple += before - len(df_bed)

    else:
        raise ValueError(f"Unsupported file extension: {file_extension}")

    ensure_nonempty(df_bed)

    df_bed["formatted_id"] = df_bed.apply(
        lambda feature: f"{feature.chrom}|{feature.start}|{feature.end}|{feature.formatted}|{feature.score}|{feature.strand}", axis=1
    )

    df_bed = df_bed[["chrom", "start", "end", "formatted_id", "score", "strand"]]

    df_bed.to_csv(bed_output_path, sep="\t", header=False, index=False)
    if dropped_missing_cf > 0:
        log_message("[WARNING]", f"Dropped {dropped_missing_cf} TE records with missing Class/Family.", color="warning")
    if dropped_simple > 0:
        log_message("[WARNING]", f"Filtered {dropped_simple} TE records from Simple_repeat/Low_complexity/Satellite.", color="warning")

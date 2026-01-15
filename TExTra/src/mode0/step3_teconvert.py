import subprocess
import os
import glob
import pandas as pd
import shutil

from util.logger import *
from util.toolkit import *

def convent_func(gtf, args):
    log_message("[INFO]", "Step 3: Convert TE and gene annotation to bed", bold=True, color="step")

    # Create the directory for storing assembly results
    convert_dir = os.path.join(args.out_dir, "convert")
    if os.path.exists(convert_dir):
        shutil.rmtree(convert_dir)
    os.makedirs(convert_dir)

    clock = get_current_time()
    print(f"{clock}\tConvert TE annotation({args.te})...")
    file_extension = os.path.splitext(args.te)[1].lower()
    te_bed_path = os.path.join(convert_dir, "TE_anno_1.bed")
    if file_extension in {".gtf", ".gff"}:
        gtf_to_bed(args.te, te_bed_path, 'TE')
    elif file_extension in {".out", ".txt", ".bed"}:
        rmsk_to_bed(args.te, te_bed_path, file_extension)
    else:
        log_message("[ERROR]", f"Please check te annotation file format (gtf/gff/out/txt/bed).", color="error")
        sys.exit(1)
    
    extend_bed_path = None
    if args.extend:
        print(f"{clock}\tConvert extended TE annotation ({args.extend})...")
        extend_extension = os.path.splitext(args.extend)[1].lower()
        extend_bed_path = os.path.join(convert_dir, "TE_anno_2.bed")

        if extend_extension in {".gtf", ".gff"}:
            gtf_to_bed(args.extend, extend_bed_path, 'TE')
        elif extend_extension in {".out", ".txt", ".bed"}:
            rmsk_to_bed(args.extend, extend_bed_path, extend_extension)
        else:
            log_message("[ERROR]", f"Please check extended TE annotation file format (gtf/gff/out/txt/bed).", color="error")
            sys.exit(1)

    # merge BED
    final_bed_path = os.path.join(convert_dir, "TE_anno.bed")
    merge_and_sort_bed(te_bed_path, extend_bed_path, final_bed_path)
    os.remove(te_bed_path)
    if extend_bed_path and os.path.exists(extend_bed_path):
        os.remove(extend_bed_path)

    log_message("[SUCCESS]", f"Converted TE BED is stored at {final_bed_path}", color="success")

    clock = get_current_time()
    print(f"{clock}\tConvert novel gene annotation({gtf})...")
    file_extension = os.path.splitext(gtf)[1].lower()
    gene_bed_path = os.path.join(convert_dir, "gene_anno.bed")
    if file_extension in {".gtf", ".gff"}:
        gtf_to_bed(gtf, gene_bed_path, 'gene')
    log_message("[SUCCESS]", f"Converted gene BED is stored at {gene_bed_path}", color="success")


def rmsk_to_bed(rmsk_path, bed_output_path, file_extension):
    if file_extension == ".txt":
        df = pd.read_csv(rmsk_path, sep="\t", header=None, comment="#")

        df_bed = df.iloc[:, [5, 6, 7, 10, 1, 9, 11, 12]].copy()
        df_bed.columns = ["chrom", "start", "end", "te_name", "score", "strand", "class", "family"]

        df_bed = df_bed.dropna(subset=["family"])
        df_bed["formatted"] = df_bed["te_name"] + ":" + df_bed["family"] + ":" + df_bed["class"]

    elif file_extension == ".out":
        df = pd.read_csv(rmsk_path, sep="\t", header=None)
        df_bed = df.iloc[:, [4, 5, 6, 9, 0, 8, 10]].copy()
        df_bed.columns = ["chrom", "start", "end", "te_name", "score", "strand", "classification"]
        df_bed[["class", "family"]] = df_bed["classification"].str.split("/", expand=True, n=1)
        df_bed = df_bed.dropna(subset=["family"])
        df_bed["formatted"] = df_bed["te_name"] + ":" + df_bed["family"] + ":" + df_bed["class"]
    
    elif file_extension == ".bed":
        df_bed = pd.read_csv(rmsk_path, sep="\t", header=None)
        df_bed.columns = ["chrom", "start", "end", "raw_formatted", "score", "strand"]

        def normalize_te_format(formatted):
            """transform 'SINE/Alu|AluSx' to 'AluSx:Alu:SINE' """
            if "|" in formatted and "/" in formatted:
                te_class, rest = formatted.split("/", 1)
                te_family, te_name = rest.split("|", 1)
                return f"{te_name}:{te_family}:{te_class}"
            return formatted 

        df_bed["formatted_id"] = df_bed["raw_formatted"].apply(normalize_te_format)

    else:
        raise ValueError(f"Unsupported file extension: {file_extension}")

    df_bed["formatted_id"] = df_bed.apply(
        lambda feature: f"{feature.chrom}|{feature.start}|{feature.end}|{feature.formatted}|{feature.score}|{feature.strand}", axis=1
    )

    df_bed = df_bed[["chrom", "start", "end", "formatted_id", "score", "strand"]]

    df_bed.to_csv(bed_output_path, sep="\t", header=False, index=False)
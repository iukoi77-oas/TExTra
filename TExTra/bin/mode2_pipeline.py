import os
import sys
import glob
import argparse
import pandas as pd
import tempfile
import subprocess as sp
from itertools import chain
from joblib import Parallel, delayed
from collections import Counter
from tempfile import NamedTemporaryFile

# Get the script directory and add the parent directory to the path
script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.insert(0, os.path.abspath(parent_dir))

from util.logger import log_message
from util.toolkit import generate_metaexon_gtf, generate_TEexon
from TExTra.src.mode2.step1_count import count_func
from TExTra.src.mode2.step2_draw import draw_func

def parse_arguments(args_list):
    parser = argparse.ArgumentParser(
        prog="TExTra quant",
        description="TExTra pipeline for quantitative analysis (quantifying TE exonization)",
        add_help=True, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input and output parameters
    parser.add_argument("-o", "--out_dir", type=str, required=True, help="Output directory for results (must match Mode1 output path).")
    
    # Computational parameters
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use (default: 4).")
    parser.add_argument("-s", "--samples", required=True, help="Samples' name, separated by commas")

    # Genome and annotation
    parser.add_argument("-r", "--te", type=str, required=True, help="Path to TE annotation GTF file.")
    parser.add_argument("--strand", type=str, choices=["none", "rf", "fr", "r", "f"], default="none", 
                        help="Strand-specific RNA-seq library type (options: none, rf, fr, r, f).")
    parser.add_argument("--readtype", choices=["paired", "single"], default="paired", help="Read type: paired or single")
    
    # Other parameters
    # parser.add_argument("--mode1", type=str, required=True, help="Path to Mode1 output (should match --out_dir).")
    parser.add_argument("--prep", required=True, help="Path to TExTra prep output")
    parser.add_argument("--qual", required=True, help="Path to TExTra qual output")
    parser.add_argument("-e","--EM", help = "Run estimation-maximization on TE counts given number of times (optional, specify 0 if no EM desired; default=auto)", type=str, default = "auto")
    parser.add_argument("--bw", action="store_true", help="Generate BigWig files (optional; default: False).")
    parser.add_argument("--normlib", choices=["RPM", "None"], default="None", 
                        help="Normalize bedgraphs by library size (optional; required if --bw is enabled).")
    parser.add_argument("--project", type=str, default="project",
                        help="Project name (default: project)")
    parser.add_argument("--keep-temp", action="store_true", help="Keep intermediate files (default: False).")
    parser.add_argument("--tempfolder", type=str, default=None, help="Path to temporary folder for intermediate files (optional).")
    
    return parser.parse_args(args_list)

def main(args_list=None):
    args_list = sys.argv[2:]
    args = parse_arguments(args_list)
    
    # Ensure output directory exists
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Parse sample list
    sample_list = args.samples.split(",")
    alignment_dir = os.path.join(args.prep, 'alignment')
    bamfiles_dict = {sample: [] for sample in sample_list}
    replicates_dict = {sample: [] for sample in sample_list}

    for folder in os.listdir(alignment_dir):
        folder_path = os.path.join(alignment_dir, folder)
        
        if os.path.isdir(folder_path):  # 只处理目录
            for sample in sample_list:
                # 生成通配符路径，例如 group1_rep1_*_accepted_hits.bam
                bam_pattern = os.path.join(folder_path, f"{sample}_*_accepted_hits.bam")
                matched_bams = glob.glob(bam_pattern)  # 获取所有匹配的 BAM 文件
                
                if matched_bams:  # 如果找到匹配的 BAM 文件
                    bamfiles_dict[sample].extend(matched_bams)
                    replicates_dict[sample].append(folder)
    
    # Step 1: Handle exon files with multiple samples
    log_message("[INFO]", "Step 1: Handling exon files with multiple samples", bold=True, color="step")
    candidate_dir = os.path.join(args.qual, "candidate")
    count_dir = os.path.join(args.out_dir, "quantification")
    os.makedirs(count_dir, exist_ok=True)
    gtf_out = generate_metaexon_gtf(sample_list, candidate_dir, count_dir, args.project)
    log_message("[SUCCESS]", f"GTF written to {gtf_out}.", color="success")
    
    # log_message("[INFO]", "Step 1: Cleaning TE annotation", bold=True, color="step")
    # annotation = TE_annotation(args.te, clean_out)
    # annotation.clean_func()
    
    # Step 2: Quantify TE and exon expression
    log_message("[INFO]", "Step 2: Quantifying the expression of exonization TE", bold=True, color="step")
    # Count the uniq reads for exon
    output_exon_file = os.path.join(count_dir, f"{args.project}_exon.txt")
    # output_exon_file = os.path.join(count_dir, f"{args.project}_te.txt")

    featurecounts_cmd = [
        "featureCounts",
        "-a", gtf_out,
        "-o", output_exon_file,
        "-t", "exon",        
        "-f", 
        "-T", str(args.threads),
        "-O", "-M", "--fraction"       
    ]

    if args.readtype == "paired":
        featurecounts_cmd += ["-p", "-B", "-C"]  
    else:
        pass
    
    if args.strand == "fr" or args.strand == "f":
        featurecounts_cmd += ["-s", str(1)]  
    elif args.strand == "rf" or args.strand == "r":
        featurecounts_cmd += ["-s", str(2)]
    else:
        featurecounts_cmd += ["-s", str(0)] 
    
    bamfiles = list(chain.from_iterable(bamfiles_dict.values()))
    featurecounts_cmd.extend(bamfiles) 
    try:
        sp.run(featurecounts_cmd, check=True)
        print("FeatureCounts ran successfully. Output saved to:", output_exon_file)
    except sp.CalledProcessError as e:
        print("Error running featureCounts:", e)

    # Count the uniq and multi reads for TE
    convert_out = os.path.join(args.prep, "convert")
    te_temp = tempfile.NamedTemporaryFile(delete=False, dir = convert_out, prefix="TE_anno" +  ".tmp")
    rmsk_bed = os.path.join(convert_out, "TE_anno.bed")
    awk_cmd = f'''awk '
        BEGIN {{ OFS = "\\t" }}
        /^$/ || /^#/ {{ next }}
        {{
            split($4, name_parts, "|");
            rgb = ($6 == "+") ? "120,91,12" : "94,137,255";
            print $1, $2, $3, 
                name_parts  [1] "|" name_parts  [2] "|" name_parts  [3] "|" name_parts  [4] "|.|" $6,
                ".", $6, $2, $3, rgb
        }}' {rmsk_bed} > {te_temp.name}'''
    sp.run(awk_cmd, shell=True, check=True)
    Parallel(n_jobs=args.threads)(delayed(count_func)(te_temp.name, args, sample, rep) for sample in sample_list for rep in replicates_dict[sample])
    os.unlink(te_temp.name)

    # Process TEexon
    quant_dir = os.path.join(args.out_dir, "quantification")
    te_files = glob.glob(os.path.join(quant_dir, "*_TEcounts.txt"))
    teexon_df = generate_TEexon(gtf_out, output_exon_file, te_files)

    output_teexon_file = os.path.join(quant_dir, f"{args.project}_TEexons.txt")
    teexon_df.to_csv(output_teexon_file, sep='\t')
    log_message("[SUCCESS]", f"The quantification of TEexon has been saved at {output_teexon_file}.", color="success")
    
    # # Step 3: Merge biological replicate count files
    # log_message("[INFO]", "Step 3: Merging the count files of different biological replicates", bold=True, color="step")
    # quant_dir = os.path.join(args.out_dir, "quantification")
    # output_file = os.path.join(quant_dir, "exonTE_total.txt")
    
    # exon_te_files = glob.glob(os.path.join(quant_dir, "*_exonTE.txt"))
    # if not exon_te_files:
    #     log_message("[ERROR]", "No _exonTE.txt files found.")
    #     sys.exit(1)
    
    # data_frames = [pd.read_csv(file, sep="\t")[["TE_ID", "gene", "Sample", "tot_count"]] for file in exon_te_files]
    # merged_df = pd.concat(data_frames)
    # pivot_df = merged_df.pivot_table(index=["TE_ID", "gene"], columns="Sample", values="tot_count", aggfunc="first").reset_index()
    # pivot_df.to_csv(output_file, sep="\t", index=False)
    
    # log_message("[SUCCESS]", f"Files merged into {output_file}", color="success")

    # Step 4: Generate BigWig files (if enabled)
    if args.bw:
        replicates = list(chain.from_iterable(replicates_dict.values()))
        log_message("[INFO]", "Step 4: Generating BigWig files", bold=True, color="step")
        Parallel(n_jobs=args.threads)(delayed(draw_func)(args, rep) for rep in replicates)
        log_message("[SUCCESS]", "The BigWig files have been completed.", color="success")

if __name__ == "__main__":
    main()
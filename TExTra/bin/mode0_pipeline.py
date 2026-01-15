import sys
import os
import argparse
import pandas as pd

# Get the script directory and add the parent directory to the path
script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.insert(0, os.path.abspath(parent_dir))

from src.mode0.step1_alignment import align_func
from src.mode0.step2_assembly import assemble_func
from src.mode0.step3_teconvert import convent_func

def parse_arguments(args_list):
    parser = argparse.ArgumentParser(
        prog="TExTra prep",
        description="TExTra pipeline for qualitative analysis (distinguishing TE-exon relationships)",
        add_help=True, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("-i", "--input", required=True, help="Path to input TSV file containing mate1, mate2, and group columns")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use")
    parser.add_argument("-g", "--genome", required=True, help="Path to genome FASTA file")
    parser.add_argument("-G", "--gene", required=True, help="Path to GTF file with gene annotations")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for results")
    parser.add_argument("-r", "--te", required=True, help="Path to TE annotation GTF file")
    parser.add_argument("--index", help="Path to genome index directory (optional)")
    parser.add_argument("--extend", help="Path to an optional extended TE annotation file (optional)")
    parser.add_argument("--is_short", action="store_true", default=True, help="Indicate if short-read sequencing is used (default: True)")
    parser.add_argument("--is_long", action="store_true", help="Indicate if long-read sequencing is used (currently not fully supported)")
    parser.add_argument("--strand", choices=["none", "rf", "fr", "r", "f"], default="none", help="Strand specificity: none, rf (reverse-forward), fr (forward-reverse), r (reverse, single-end), f (forward, single-end)")
    parser.add_argument("--readtype", choices=["paired", "single"], default="paired", help="Read type: paired or single")
    parser.add_argument("--keep-temp", action="store_true", help="Keep intermediate files (default: False).")
    parser.add_argument(
        "--taco-disable", 
        action="store_true", 
        help="Disable the use of Taco for GTF merging and use Stringtie instead (default: False)"
    )
    parser.add_argument(
        "--de-novo-disable", 
        action="store_true", 
        help="Disable the use of de novo for Assembly (default: False)"
    )
    parser.add_argument("--best", action="store_true", help="Use optimal parameters to merge assemblies.")

    return parser.parse_args(args_list)

def main(args_list=None):
    args_list = sys.argv[2:]
    args = parse_arguments(args_list)
    
    # Notify user if --is_long is used
    if args.is_long:
        print("Warning: --is_long option is not fully supported yet.")
    
    # Read input data
    input_df = pd.read_csv(args.input, sep="\t", header=None)

    # Assign meaningful column names
    input_df.columns = ["samples"] + [f"rep{i+1}" for i in range(input_df.shape[1] - 1)]
    input_df = input_df.replace(["NA", "NaN", "nan", ""], pd.NA)

    # Extract replicate files (all remaining columns) and convert to strings
    samples = input_df.iloc[:, 0]
    replicates = input_df.iloc[:, 1:].astype(str)

    # Lists to store processed mate1, mate2, and group values
    mate1_list, mate2_list, group_list = [], [], []

    # Iterate over each sample and its corresponding replicate files
    # for sample, files in zip(samples, replicates.values):
    #     for i, file in enumerate(files, 1):  # Enumerate to get replicate index (1-based)
    #         if pd.isna(file) or file == "nan":  # Skip empty values
    #             continue
    #         if "," in file:  # If paired-end sequencing, split mate1 and mate2
    #             mate1, mate2 = file.split(",")
    #         else:  # If single-end sequencing, set mate2 as "-"
    #             mate1, mate2 = file, "-"
            
    #         # Format group as sample_rep1, sample_rep2, etc.
    #         group = f"{sample}_rep{i}"
            
    #         # Store values in respective lists
    #         mate1_list.append(mate1)
    #         mate2_list.append(mate2)
    #         group_list.append(group)

    for _, row in input_df.iterrows():
        sample = row.iloc[0] 

        rep_idx = 0  # >>> NEW: real replicate index per sample

        for file in row.iloc[1:]:
            if pd.isna(file):
                continue

            rep_idx += 1  # >>> MODIFIED: increment only for real replicates

            if "," in str(file):
                mate1, mate2 = file.split(",", 1)
            else:
                mate1, mate2 = file, "-"

            # =========================
            # >>> MODIFIED: biologically correct replicate naming
            # =========================
            run_id = f"{sample}_rep{rep_idx}"

            mate1_list.append(mate1)
            mate2_list.append(mate2)
            group_list.append(run_id)

    # Create the final DataFrame
    input_info = pd.DataFrame({
        "mate1": mate1_list,
        "mate2": mate2_list,
        "sample": group_list
    })
    
    # Execute pipeline steps
    bamfiles, samples = align_func(input_info, args)
    novel_gtf = assemble_func(bamfiles, samples, args)
    convent_func(novel_gtf, args)


if __name__ == "__main__":
    main()
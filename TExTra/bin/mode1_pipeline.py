import sys
import os
import glob
import argparse
import shutil
import pandas as pd
from joblib import Parallel, delayed

# Get the script directory and add the parent directory to the path
script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.insert(0, os.path.abspath(parent_dir))

from util.toolkit import *
from util.toolkit import *
from util.HITindex import *

from TExTra.src.mode1.step1_classify import classify_func
from TExTra.src.mode1.step2_candidate import candidate_func

def parse_arguments(args_list):
    parser = argparse.ArgumentParser(
        prog="TExTra qual",
        description="TExTra pipeline for qualitative analysis (distinguishing TE-exon relationships)",
        add_help=True, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # parser.add_argument("-i", "--input", required=True, help="Path to input TSV file containing mate1, mate2, and group columns")
    parser.add_argument("-s", "--samples", required=True, help="Samples' name, separated by commas")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use")
    parser.add_argument("--prep", required=True, help="Path to TExTra prep output")
    # parser.add_argument("-g", "--genome", required=True, help="Path to genome FASTA file")
    # parser.add_argument("-G", "--gene", required=True, help="Path to GTF file with gene annotations")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for results")
    # parser.add_argument("-r", "--te", required=True, help="Path to TE annotation GTF file")
    # parser.add_argument("--index_dir", help="Path to genome index directory (optional)")
    parser.add_argument("--is_short", action="store_true", default=True, help="Indicate if short-read sequencing is used (default: True)")
    parser.add_argument("--is_long", action="store_true", help="Indicate if long-read sequencing is used (currently not fully supported)")
    parser.add_argument("--strand", choices=["none", "rf", "fr", "r", "f"], default="none", help="Strand specificity: none, rf (reverse-forward), fr (forward-reverse), r (reverse, single-end), f (forward, single-end)")
    parser.add_argument("--ss3buffer", type=int, default=20, help="3' splice site buffer size")
    parser.add_argument("--ss5buffer", type=int, default=50, help="5' splice site buffer size")
    parser.add_argument("--readtype", choices=["paired", "single"], default="paired", help="Read type: paired or single")
    parser.add_argument("--threshold", type=float, default=0.75, help="Threshold for candidate selection")
    parser.add_argument("--keep-temp", action="store_true", help="Keep intermediate files.")
    
    return parser.parse_args(args_list)

def main(args_list=None):
    args_list = sys.argv[2:]
    args = parse_arguments(args_list)
    
    # Notify user if --is_long is used
    if args.is_long:
        print("Warning: --is_long option is not fully supported yet.")
    
    # # Read input data
    # input_info = pd.read_csv(args.input, header=None, sep="\t", usecols=[0,1,2], names=['mate1', 'mate2', 'group'])
    
    # # Execute pipeline steps
    # bamfiles, samples = align_func(input_info, args)
    # novel_gtf, merge_novel_ref_gtf = assemble_func(bamfiles, samples, args)
    sample_list = [s.strip() for s in args.samples.split(",")]

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
    
    exon_bed = os.path.join(args.prep, 'convert', 'gene_anno.bed')
    novel_gtf = os.path.join(args.prep, 'assembly', 'novel_transcripts.gtf')

    log_message("[INFO]", "Step 3: Classify exons", bold=True, color="step")

    # Create the directory for storing classify results
    classify_dir = os.path.join(args.out_dir, "classify")
    # if os.path.exists(classify_dir):
    #     shutil.rmtree(classify_dir)
    # os.makedirs(classify_dir)
    clock = get_current_time()
    print(f"{clock}\tCreated 'classify' directory at {classify_dir}")
    annotation_dir = os.path.join(classify_dir, "annotation")
    os.makedirs(annotation_dir, exist_ok=True)

    print('Step3-1: Generate metaexon BED')
    clock = get_current_time()
    print(f"{clock}\tProcessing merged BED file...")
    input_buffer_bed = metaexon_bed(exon_bed, annotation_dir, args)
    log_message("[SUCCESS]", f"Generated metaexon BED is stored at {annotation_dir}.", color="success")

    hitindex_dir = os.path.join(classify_dir, "HITindex")

    log_message("[INFO]", "Step 4: Identify candidate TE exonization", bold=True, color="step")
    
    # Create the directory for storing candidate results
    candidate_dir = os.path.join(args.out_dir, "candidate")
    if os.path.exists(candidate_dir):
        shutil.rmtree(candidate_dir)
    os.makedirs(candidate_dir)
    clock = get_current_time()
    print(f"{clock}\tCreated 'candidate' directory at {candidate_dir}")

    for sample in sample_list:
        group_dir = os.path.join(candidate_dir, sample)
        tmp_dir = os.path.join(group_dir, "tmp")
        os.makedirs(tmp_dir, exist_ok=True)

    def process_sample(sample):
        bamfiles = bamfiles_dict[sample]
        replicates = replicates_dict[sample]

        # classified_exons = [
        #     f"{hitindex_dir}/{rep}.exon" for rep in replicates
        # ]

        classified_exons = classify_func(input_buffer_bed, bamfiles, replicates, args)
        candidate_func(classified_exons, bamfiles, replicates, sample, novel_gtf, args)

    Parallel(n_jobs=args.threads, verbose=5)(delayed(process_sample)(sample) for sample in sample_list)

if __name__ == "__main__":
    main()
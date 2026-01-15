import os
import pybedtools
import shutil
from util.toolkit import *
from util.logger import *

def candidate_func(classified_exons, bamfiles, samples, group, merge_novel_ref_gtf, args):
    log_message("[INFO]", "Step 4: Identify candidate TE exonization", bold=True, color="step")
    
    # Create the directory for storing candidate results
    candidate_dir = os.path.join(args.out_dir, "candidate")
    # if os.path.exists(candidate_dir):
    #     shutil.rmtree(candidate_dir)
    # os.makedirs(candidate_dir)
    # clock = get_current_time()
    # print(f"{clock}\tCreated 'candidate' directory at {candidate_dir}")
    tmp_dir = os.path.join(candidate_dir, f"{group}/tmp")
    # os.makedirs(tmp_dir, exist_ok=True)
    # print(tmp_dir)
    if not os.path.exists(tmp_dir):
        raise OSError(f"Directory {tmp_dir} was not created successfully.")

    # print('Step4-1: Convert TE GTF to BED')
    # clock = get_current_time()
    # print(f"{clock}\tProcessing TE annotation GTF file...")
    # te_output_bed = os.path.join(tmp_dir, "TE_converted.bed")
    # gtf_to_bed(args.te, te_output_bed, 'TE')
    # log_message("[SUCCESS]", f"Converted TE BED is stored at {te_output_bed}.", color="success")

    print('Step4-1: Filter exonerated TEs')
    clock = get_current_time()
    print(f"{clock}\tProcessing biological replicates of classified exons...")
    merged_classified_exons = os.path.join(tmp_dir, "merged_classified_exons.bed")
    merged_exon(classified_exons, samples, merged_classified_exons, args.threshold)
    log_message("[SUCCESS]", f"Merged HIT exons is stored at {merged_classified_exons}.", color="success")
    clock = get_current_time()
    print(f"{clock}\tScanning for TEs undergoing exonization...")
    te_output_bed = os.path.join(args.prep, 'convert', 'TE_anno.bed')
    scan_TE_exon(bamfiles, samples, te_output_bed, merge_novel_ref_gtf, merged_classified_exons, tmp_dir, args)
    log_message("[SUCCESS]", f"Exonization TEs are saved to {candidate_dir}.", color="success")
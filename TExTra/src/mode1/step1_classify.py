import os
import shutil
from util.logger import *
from util.HITindex import *
from util.toolkit import *
import warnings
warnings.simplefilter('ignore')

def classify_func(buffer_exon_bed, bamfiles, samples, args):
    # log_message("[INFO]", "Step 3: Classify exons", bold=True, color="step")

    # # Create the directory for storing classify results
    # classify_dir = os.path.join(args.out_dir, "classify")
    # if os.path.exists(classify_dir):
    #     shutil.rmtree(classify_dir)
    # os.makedirs(classify_dir)
    # clock = get_current_time()
    # print(f"{clock}\tCreated 'classify' directory at {classify_dir}")
    # annotation_dir = os.path.join(classify_dir, "annotation")
    # os.makedirs(annotation_dir, exist_ok=True)
    # tmp_dir = os.path.join(annotation_dir, "tmp")
    # os.makedirs(tmp_dir, exist_ok=True)

    # print('Step3-1: Convert GTF to BED')
    # clock = get_current_time()
    # print(f"{clock}\tProcessing novel GTF file...")
    # novel_exon_gtf = os.path.join(tmp_dir, "novel_exon.gtf")
    # filter_exon_from_gtf(merged_gtf, novel_exon_gtf)
    # novel_output_bed = os.path.join(tmp_dir, "novel_exon.bed")
    # gtf_to_bed(novel_exon_gtf, novel_output_bed, 'gene')
    # print(f"{clock}\tProcessing genome GTF file...")
    # genome_exon_gtf = os.path.join(tmp_dir, "genome_exon.gtf")
    # filter_exon_from_gtf(args.gene, genome_exon_gtf)
    # genome_output_bed = os.path.join(tmp_dir, "genome_exon.bed")
    # gtf_to_bed(genome_exon_gtf, genome_output_bed, 'gene')
    # print(f"{clock}\tMerging gene GTF files...")
    # output_gtf_bed = os.path.join(annotation_dir, "raw_exon.bed")
    # merge_and_sort_bed(novel_output_bed, genome_output_bed, output_gtf_bed)
    # log_message("[SUCCESS]", f"Step3-1 is done.", color="success")
    
    annotation_dir = os.path.join(args.out_dir, "classify", "annotation")

    # print('Step3-2: Generate metaexon BED')
    # clock = get_current_time()
    # print(f"{clock}\tProcessing merged BED file...")
    # input_buffer_bed = os.path.join(annotation_dir, f'buffer_ss3-{args.ss3buffer}_ss5-{args.ss5buffer}_metaexon.bed')
    # log_message("[SUCCESS]", f"Generated metaexon BED is stored at {annotation_dir}.", color="success")

    print('Step3-3: Compute HITindex of exons and classify')
    clock = get_current_time()
    print(f"{clock}\tRunning the HITindex pipeline...")
    classify_dir = os.path.join(args.out_dir, "classify")
    output_hitindex_dir = os.path.join(classify_dir, "HITindex")
    os.makedirs(output_hitindex_dir, exist_ok=True)
    exon_outnames = HITindex_pipeline(buffer_exon_bed, bamfiles, samples, output_hitindex_dir, args)
    log_message("[SUCCESS]", f"The results of the HITindex-based calculations and classifications are stored at {output_hitindex_dir}.", color="success")

    return exon_outnames
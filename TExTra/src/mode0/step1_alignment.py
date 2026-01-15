import subprocess
import os
import pandas as pd
import shutil

from joblib import Parallel, delayed
from util.logger import *

def align_func(input_info, args):
    """
    Perform STAR alignment with parallel processing.
    """
    log_message("[INFO]", "Step 1: Mapping input files to genomes", bold=True, color="step")

    alignment_dir = os.path.join(args.out_dir, "alignment")
    os.makedirs(alignment_dir, exist_ok=True)
    clock = get_current_time()
    print(f"{clock}\tCreated 'alignment' directory at {alignment_dir}")
    
    # if all(input_info.iloc[:, 0:2].apply(lambda x: all([file.endswith('.bam') or file == '-' for file in x]), axis=1)):
    #     log_message("[WARNING]", "BAM files already exist. Skip Step1-1 and Step1-2. Sorting BAM files...", color="warning")
    #     bam_files = [file for file in input_info[['mate1']].values.flatten() if isinstance(file, str) and file.endswith('.bam')]
    #     samples = input_info['sample']
    #     sorted_bams = Parallel(n_jobs=args.threads)(
    #         delayed(sort_bam)(bam, sample, alignment_dir, args) for bam, sample in zip(bam_files, samples) if bam.endswith('.bam')
    #     )
    #     return process_bam(sorted_bams, samples, args), samples

    input_info['mate1'] = input_info['mate1'].astype(str).str.strip()
    input_info['mate2'] = input_info['mate2'].astype(str).str.strip()

    def is_bam_file(file):
        if pd.isna(file):
            return True
        file_str = str(file).strip() 
        if file_str in ('', '-', 'nan'):
            return True
        result = file_str.endswith('.bam')
        return result

    if all(input_info[['mate1', 'mate2']].apply(lambda row: all(is_bam_file(f) for f in row), axis=1)):
        log_message("[WARNING]", "BAM files already exist. Skip STAR alignment...", color="warning")
        bam_files = []
        for _, row in input_info.iterrows():
            for f in [row['mate1'], row['mate2']]:
                if is_bam_file(f) and f != '-' and not pd.isna(f):
                    f_clean = str(f).strip()
                    if f_clean != '':
                        bam_files.append(f_clean)

        samples = input_info['sample'].tolist()
        sorted_bams = Parallel(n_jobs=args.threads)(
            delayed(sort_bam)(bam, sample, alignment_dir, args) for bam, sample in zip(bam_files, samples)
        )
        return process_bam(sorted_bams, samples, args), samples
    
    # Generate genome index if needed
    if not args.index:
        args.index = os.path.join(alignment_dir, 'index')
    if not os.path.exists(os.path.join(args.index, 'SAindex')):
        clock = get_current_time()
        print(f"{clock}\tCreating index using STAR...")
        subprocess.call(["STAR", "--runThreadN", str(args.threads), "--runMode", "genomeGenerate", 
                         "--genomeDir", args.index, "--genomeFastaFiles", args.genome, 
                         "--sjdbGTFfile", args.gene, "--sjdbOverhang", "100"], 
                         stdout=subprocess.DEVNULL)
    
    # Run STAR alignment in parallel
    print('Step1-2: Performing alignment')
    clock = get_current_time()
    print(f"{clock}\tRunning STAR alignment...")
    bam_files = Parallel(n_jobs=args.threads)(
        delayed(run_star)(row, alignment_dir, args) for _, row in input_info.iterrows()
    )
    log_message("[SUCCESS]", "Step1-2 is done.", color="success")

    sample_names = input_info['sample'].tolist()
    output_bams = process_bam(bam_files, sample_names, args)
    print(bam_files)
    print(sample_names)

    # If keep_temp is False, delete unnecessary files
    if not args.keep_temp:
        delete_temp_files(bam_files, alignment_dir, input_info['sample'])

    return output_bams, input_info['sample']


def delete_temp_files(bam_files, alignment_dir, samples):
    """
    Delete intermediate files after alignment if keep_temp is False.
    """
    for bam, sample in zip(bam_files, samples):
        sample_folder = os.path.dirname(bam)
        for file in os.listdir(sample_folder):
            file_path = os.path.join(sample_folder, file)
            if not (file.endswith("accepted_hits.bam") or file == "Aligned.sortedByCoord.out.bam"):
                try:
                    if os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                        log_message("[INFO]", f"Deleted temporary folder: {file_path}", color="info")
                    else:
                        os.remove(file_path)
                        log_message("[INFO]", f"Deleted temporary file: {file_path}", color="info")
                except Exception as e:
                    log_message("[ERROR]", f"Error deleting {file_path}: {str(e)}", color="error")


def sort_bam(bam_file, sample, alignment_dir, args):
    sample_folder = os.path.join(alignment_dir, sample)
    os.makedirs(sample_folder, exist_ok=True)
    sorted_bam = os.path.join(sample_folder, "Aligned.sortedByCoord.out.bam")
    subprocess.run(["samtools", "sort", "-@", str(args.threads), "-o", sorted_bam, bam_file], check=True)
    
    if not args.keep_temp:
        delete_temp_files([bam_file], alignment_dir, [sample])
    
    return sorted_bam

def run_star(row, alignment_dir, args):
    mate1 = row['mate1']
    mate2 = row['mate2']
    sample = row['sample']
    
    mate2 = '' if mate2 == '-' else mate2
    output_prefix = os.path.join(alignment_dir, f"{sample}/")
    os.makedirs(output_prefix, exist_ok=True)
    
    star_cmd = ["STAR", "--genomeDir", args.index, "--runThreadN", str(args.threads), 
                "--twopassMode", "Basic", "--readFilesCommand", "zcat", "--readFilesIn", mate1, mate2,
                "--outSAMtype", "BAM", "SortedByCoordinate", "--outFileNamePrefix", output_prefix,
                "--chimOutType", "Junctions", "--chimSegmentMin", "10", "--outBAMsortingBinsN", "200"]
    subprocess.call(star_cmd, stdout=subprocess.DEVNULL)
    
    if not args.keep_temp:
        delete_temp_files([os.path.join(output_prefix, f) for f in os.listdir(output_prefix) 
                           if f != "Aligned.sortedByCoord.out.bam"], alignment_dir, [sample])
    
    return os.path.join(output_prefix, 'Aligned.sortedByCoord.out.bam')


def process_bam(bam_files, samples, args):
    print('Step1-3: Filtering BAM files')
    clock = get_current_time()
    print(f"{clock}\tProcessing BAM files...")
    output_bams = Parallel(n_jobs=args.threads)(
        delayed(filter_bam)(bam, sample, args) for bam, sample in zip(bam_files, samples)
    )
    log_message("[SUCCESS]", "Step1-3 is done.", color="success")
    return output_bams

def filter_bam(bam_file, sample, args):
    bam_dir = os.path.dirname(bam_file)
    output_bam = os.path.join(bam_dir, f"{sample}_accepted_hits.bam")
    subprocess.call(["samtools", "view", "-@", str(args.threads), "-b", "-F", "4", bam_file], stdout=open(output_bam, 'w'))
    return output_bam
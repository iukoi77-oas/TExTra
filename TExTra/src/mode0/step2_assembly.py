import subprocess
import os
import glob
import pandas as pd
import shutil

from util.logger import *
from util.toolkit import *
from joblib import Parallel, delayed

def assemble_func(bamfiles, samples, args):
    log_message("[INFO]", "Step 2: Assemble transcripts", bold=True, color="step")

    # Create the directory for storing assembly results
    assembly_dir = os.path.join(args.out_dir, "assembly")
    if os.path.exists(assembly_dir):
        shutil.rmtree(assembly_dir)
    os.makedirs(assembly_dir)
    clock = get_current_time()
    print(f"{clock}\tCreated 'assembly' directory at {assembly_dir}")
    tmp_dir = os.path.join(assembly_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    print('Step2-1: Running Stringtie')
    clock = get_current_time()
    print(f"{clock}\tProcessing BAM files...")
    # Construct the Stringtie command based on different options
    if args.is_long:
        stringtie_cmd = ['stringtie', '-L', '-a', '3']  # For long gene models
    elif args.strand == 'rf':
        stringtie_cmd = ['stringtie', '--rf']  # For reverse-forward strand
    elif args.strand == 'fr':
        stringtie_cmd = ['stringtie', '--fr']  # For forward-reverse strand
    else:
        stringtie_cmd = ['stringtie']  # Default Stringtie command

    # Parallelizing Stringtie for each BAM file
    output_gtfs = Parallel(n_jobs=args.threads)(
        delayed(run_stringtie)(bam, sample, stringtie_cmd, tmp_dir, args) for bam, sample in zip(bamfiles, samples)
    )

    # Filter out None values in case of errors during the StringTie process
    output_gtfs = [gtf for gtf in output_gtfs if gtf is not None]

    # Write all the output_gtfs to a txt file
    output_gtfs_txt = os.path.join(assembly_dir, "output_gtfs.txt")
    with open(output_gtfs_txt, 'w') as f:
        for gtf in output_gtfs:
            f.write(f"{os.path.abspath(gtf)}\n")
    
    log_message("[SUCCESS]", f"Step2-1 is done. All output GTF files have been written to {output_gtfs_txt}", color="success")

    print('Step2-2-1: Running GTFs Merge')
    clock = get_current_time()
    # Check if --taco-disable flag is set; if not, use taco for merge
    if args.taco_disable:
        print(f"{clock}\tJoining GTF files by Stringtie merge...")
        merge_command = [
            'stringtie', '--merge', '-i',
            '-p', str(args.threads),  # Number of threads
            '-o', os.path.join(tmp_dir, "stringtie_merged.gtf")
        ] + output_gtfs

        # Use optimal parameters to merge assemblies
        if args.best:
            if args.de_novo_disable:
                merge_command.extend(['-F', '3.48', '-m', '1000', '-f', '0.311']) 
            else:
                merge_command.extend(['-F', '7.19', '-m', '1000', '-f', '0.345'])

        output_merge_gtf = os.path.join(tmp_dir, "stringtie_merged.gtf")
    else:
        # Use taco tool from the util directory (assuming taco is in parent util folder)
        print(f"{clock}\tJoining GTF files by TACO...")
        # Search for the taco folder dynamically (look for any folder containing "taco")
        taco_dirs = glob.glob(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), 'util', 'taco*'))
        taco_dirs = [d for d in taco_dirs if os.path.isdir(d)]
        
        if not taco_dirs:
            log_message("[ERROR]", f"Taco tool not found in the expected location({os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), 'util')}).", color="error")
            return None
        
        # Assuming we want the first taco directory found
        taco_dir = taco_dirs[0]
        taco_path = os.path.join(taco_dir, 'taco_run')

        if not os.path.exists(taco_path):
            log_message("[ERROR]", f"Taco script not found in {taco_dir}.", color="error")
            return None
        
        print(f"Taco tool found at {taco_path}")
        
        # Run taco with the combined input GTF file
        merge_command = [
            f'{taco_path}', '-p', str(args.threads)
        ]

        # Use optimal parameters to merge assemblies
        if args.best:
            if args.de_novo_disable:
                merge_command.extend(['--filter-min-expr', '11.7', "--filter-min-length", "500", "--isoform-frac", "0.380"])
            else:
                merge_command.extend(['--filter-min-expr', '24.2', "--filter-min-length", "500", "--isoform-frac", "0.345"])

        # Run taco
        output_gtfs_txt =os.path.abspath(output_gtfs_txt)
        merge_command.extend([output_gtfs_txt])
        output_merge_gtf = os.path.join(tmp_dir, 'output/assembly.gtf')

    try:
        subprocess.run(merge_command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, cwd=tmp_dir)
        log_message("[SUCCESS]", "Step2-1 is done.", color="success")
    except subprocess.CalledProcessError as e:
        log_message("[ERROR]", f"Error during GTF merge. Error: {e.stderr.decode()}", color="error")

    print('Step2-2: Running gffcompare')
    clock = get_current_time()
    print(f"{clock}\tComparing GTF files...")
    output_gffcompare_prefix = os.path.join(tmp_dir, "gffcompare")
    subprocess.run([
        'gffcompare', '-r', str(args.gene),
        '-o', str(output_gffcompare_prefix),
        str(output_merge_gtf)
    ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    log_message("[SUCCESS]", "Step2-2 is done.", color="success")

    print('Step2-3: Filtering transcripts')
    clock = get_current_time()
    print(f"{clock}\tProcessing gffcompared GTF file...")
    output_final_gtf = os.path.join(assembly_dir, "novel_transcripts.gtf")
    input_tmap = glob.glob(os.path.join(os.path.dirname(output_merge_gtf), 'gffcompare*tmap'))
    input_gtf = glob.glob(os.path.join(tmp_dir, 'gffcompare*gtf'))
    modify_gtf_with_mapping(input_tmap[0], input_gtf[0], output_final_gtf)
    log_message("[SUCCESS]", "Step2-3 is done.", color="success")
        
    return output_final_gtf


def run_stringtie(bam, sample, stringtie_cmd, tmp_dir, args):
    output_gtf = os.path.join(tmp_dir, f"{sample}.gtf")
    
    stringtie_command = [
        *stringtie_cmd, '-u',
        '-p', str(args.threads),
        '-M', str(1.0),
        '-o', output_gtf
    ]

    if args.de_novo_disable and hasattr(args, 'gene'):
        stringtie_command.extend(['-G', args.gene])
    
    stringtie_command.append(bam)

    try:
        subprocess.run(stringtie_command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        return output_gtf
    except subprocess.CalledProcessError as e:
        log_message("[ERROR]", f"Error during StringTie assembly for {sample}. Error: {e.stderr.decode()}", color="error")
        return None

import os
import re
import glob
import tempfile
import subprocess as sp
import pandas as pd
from datetime import datetime

from util.logger import *
from util.SQuIRE_toolkit import *
from util.toolkit import refine_result

def DE_func(args):
    # group1 = args.sample1_outdir
    # group2 = args.sample2_outdir
    outfolder=args.out_dir
    verbosity=True
    pthreads= args.threads
    table_only=args.table_only
    projectname=args.project
    verbosity = True
    label_no = args.label_no
    quant_dir = args.quant

    os.makedirs(outfolder, exist_ok=True)
    DE_out = os.path.join(outfolder, "DE")
    os.makedirs(DE_out, exist_ok=True)

    log_message("[INFO]", "Step 1: Undergoing differential exonTE expression analysis...", bold=True, color="step")

    if verbosity:
        CallTime = datetime.now()
        print("Script start time is:" + str(CallTime) + '\n', file = sys.stderr)# Prints Call time
        print("Script Arguments" + '\n' + "=================", file = sys.stderr)
        args_dict = vars(args)
        for option,arg in args_dict.items():
            print(str(option) + "=" + str(arg), file = sys.stderr) #prints all arguments to std err
        print("\n", file = sys.stderr)   

    # group1_exonTE=find_file(f'{group1}/quantification','exonTE_total.txt', False, 1, True)
    # group2_exonTE=find_file(f'{group2}/quantification','exonTE_total.txt', False, 1, True)
    
    project_TEexon = find_file(f'{quant_dir}/quantification','TEexons.txt', projectname, 1, True)
    group1_name, group2_name = args.samples.split(",")

    df = pd.read_csv(project_TEexon, sep="\t")

    # Expect first column to be metaexon
    all_samples = list(df.columns[1:])

    # Infer condition as everything before the last underscore
    # e.g. heart_postnatal_14d_rep1 -> heart_postnatal_14d
    def infer_condition(sample_name):
        parts = sample_name.split("_")
        if len(parts) < 2:
            raise ValueError(
                f"Sample name '{sample_name}' does not follow <condition>_<replicate> format"
            )
        return "_".join(parts[:-1])

    sample_to_condition = {
        s: infer_condition(s) for s in all_samples
    }

    # Build coldata
    coldata = pd.DataFrame({
        "sample": list(sample_to_condition.keys()),
        "condition": list(sample_to_condition.values())
    })

    # Select samples belonging to the two requested groups
    selected_coldata = coldata[
        coldata["condition"].isin([group1_name, group2_name])
    ].reset_index(drop=True)

    if selected_coldata.empty:
        raise RuntimeError(
            f"No samples matched requested groups: {group1_name}, {group2_name}\n"
            f"Available conditions: {sorted(coldata['condition'].unique())}"
        )

    # Ensure each group has at least one sample
    counts_per_group = selected_coldata["condition"].value_counts()
    for g in [group1_name, group2_name]:
        if g not in counts_per_group:
            raise RuntimeError(
                f"No samples found for condition '{g}'. "
                f"Check --samples argument and sample naming."
            )

    selected_samples = selected_coldata["sample"].tolist()

    # -------------------------------
    # Construct DESeq2 input matrix
    # -------------------------------
    deseq_df = df[["metaexon"] + selected_samples]

    # Sanity check: non-zero counts
    if deseq_df[selected_samples].sum().sum() == 0:
        raise RuntimeError(
            "All selected samples have zero counts. "
            "Check quantification results before running DESeq2."
        )

    # Write outputs
    deseq_input = os.path.join(DE_out, "DESeq2_input.txt")
    deseq_df.to_csv(deseq_input, sep="\t", index=False)

    col_output_file = os.path.join(DE_out, "coldata.txt")
    selected_coldata.to_csv(col_output_file, sep="\t", index=False)

    print(f"Saved DESeq2 input matrix: {deseq_input}")
    print(f"Saved sample metadata: {col_output_file}")  

    # df = pd.read_csv(project_TEexon, sep='\t')
    # sample_names = list(df.columns[1:])
    # conditions = [name.split('_')[0] for name in sample_names]

    # coldata = pd.DataFrame({
    #     'sample': sample_names,
    #     'condition': conditions
    # })

    # selected_coldata = coldata[coldata["condition"].isin([group1_name, group2_name])].reset_index(drop=True)
    # selected_samples = selected_coldata["sample"].tolist()

    # deseq_df = df[["metaexon"] + selected_samples]  
    # deseq_input = os.path.join(DE_out, "DESeq2_input.txt")
    # deseq_df.to_csv(deseq_input, sep="\t", index=False)
    # print(f"Saved merged file: {deseq_input}")

    # col_output_file = os.path.join(DE_out, "coldata.txt")
    # selected_coldata.to_csv(col_output_file, sep="\t", index=False)
    # print(f"Saved sample metadata: {col_output_file}")

    create_rscript(deseq_input,col_output_file,DE_out,projectname,verbosity,str(pthreads),group1_name,group2_name,label_no,table_only)

    deseq_out = find_file(DE_out,'DESeq2_all.txt', False, 1, True)
    metaexon = find_file(f'{quant_dir}/quantification','metaexon.bed', projectname, 1, True)
    refine_result(deseq_out, metaexon)

    ####### STOP TIMING SCRIPT #######################
    if verbosity:
        print("finished writing outputs at "+ str(datetime.now()) + "\n",file = sys.stderr)

        endTime = datetime.now()
        print('end time is: '+ str(endTime) + "\n", file = sys.stderr)
        print('it took: ' + str(endTime-CallTime) + "\n", file = sys.stderr)

    log_message("[SUCCESS]", f"Differential exonTE expression results were saved in {DE_out}", color="success")
import os
import re
import glob
import tempfile
import subprocess as sp
import pandas as pd
from datetime import datetime
from Bio import SeqIO

from util.logger import *
from util.SQuIRE_toolkit import *
from util.toolkit import generate_plek_annotation_result

def ncPred_func(args):
    pthreads= args.threads
    log2fc = args.log2fc
    padj = args.padj
    projectname = args.project
    quant_dir = args.quant
    prep_dir = args.prep
    genome_fasta = args.genome
    model = args.model

    ncPred_out = os.path.join(args.out_dir, "ncPred")
    os.makedirs(ncPred_out, exist_ok=True)

    DE_out = os.path.join(args.out_dir, "DE")
    DE_exonTE = find_file(DE_out,'DESeq2_all.txt', False, 1, True)

    total_exonTE = pd.read_csv(DE_exonTE, sep="\t")
    sig_exonTE = total_exonTE[(total_exonTE["padj"] < padj) & (abs(total_exonTE["log2FoldChange"]) > log2fc)]

    metaexon_bed_path = find_file(f'{quant_dir}/quantification','metaexon.bed', projectname, 1, True)
    meta = pd.read_csv(metaexon_bed_path, sep="\t", header=None)
    meta.columns = ["chr", "start", "end", "exon_info", "strand", "te_info"]
    meta["coord"] = meta["chr"] + ":" + (meta["start"]+1).astype(str) + "-" + meta["end"].astype(str) + ":" + meta["strand"]

    sig_ids = set(sig_exonTE["coord"]) 
    meta_sig = meta[meta["coord"].apply(lambda x: any(i in x for i in sig_ids))]

    def extract_transcripts(exon_info_str):
        items = exon_info_str.split(",")
        return [i.split(":")[1] for i in items]

    meta_sig["transcript_id"] = meta_sig["exon_info"].apply(extract_transcripts)
    meta_sig["gene_id"] = meta_sig["exon_info"].apply(lambda x: list(set(i.split(":")[0] for i in x.split(","))))

    all_transcripts = set(t for sublist in meta_sig["transcript_id"] for t in sublist)
    sig_transcripts = os.path.join(ncPred_out, f"sig_transcripts.txt")

    with open(sig_transcripts, "w") as f:
        f.write("\n".join(all_transcripts))

    gtf_input = find_file(os.path.join(prep_dir, 'assembly'),'transcripts.gtf', False, 1, True)
    gtf_output = os.path.join(ncPred_out, f"sig_transcripts.gtf")
    # fgrep_command = f"fgrep -f {sig_transcripts} {gtf_input} > {gtf_output}"
    awk_cmd = f"awk '{{print \"transcript_id \\\"\" $1 \"\\\"\"}}' {sig_transcripts}"
    grep_cmd = f"grep -F -f <({awk_cmd}) {gtf_input} > {gtf_output}"
    sp.run(grep_cmd, shell=True, check=True, executable="/bin/bash")

    # sp.run(fgrep_command, shell=True, check=True)
        
    fasta_output = os.path.join(ncPred_out, "sig_transcripts.fasta")
    gffread_command = f"gffread {gtf_output} -g {genome_fasta} -w {fasta_output}"
    sp.run(gffread_command, shell=True, check=True)

    plek_dirs = glob.glob(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), 'util', 'PLEK*'))
    plek_dirs = [d for d in plek_dirs if os.path.isdir(d)]
    
    if not plek_dirs:
        log_message("[ERROR]", f"PLEK tool not found in the expected location({os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), 'util')}).", color="error")
        return None
    
    # Assuming we want the first taco directory found
    plek_dir = plek_dirs[0]
    plek_path = os.path.join(plek_dir, 'PLEK2.py')

    if not os.path.exists(plek_path):
        log_message("[ERROR]", f"PLEK2 script not found in {plek_path}.", color="error")
        return None
    
    print(f"PLEK2 tool found at {plek_path}")

    fasta_output = os.path.abspath(fasta_output)
    plek_command = f"python PLEK2.py -i {fasta_output} -m {model}"
    sp.run(plek_command, shell=True, check=True, cwd=plek_dir)

    plek_result = os.path.join(plek_dir, 'results')
    plek_final = os.path.join(ncPred_out, 'plek_final_result.csv')
    te_bed_path = find_file(f'{quant_dir}/quantification', '_TEexon.bed', projectname, 1, True)

    result_df = generate_plek_annotation_result(plek_result,sig_transcripts,meta_sig,te_bed_path,plek_final)

    log_message("[SUCCESS]", f"PLEK2 result was saved at {plek_final}.", color="success")

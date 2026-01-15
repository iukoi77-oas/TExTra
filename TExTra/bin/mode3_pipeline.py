#!/usr/bin/env python
import sys
import os
import argparse
import pandas as pd

script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.insert(0, os.path.abspath(parent_dir))

from util.logger import *        
from src.mode3.DE import DE_func  
from TExTra.src.mode3.ncPred import ncPred_func 

def parse_arguments(args_list):
    parser = argparse.ArgumentParser(
        prog="TExTra diff",
        description="Differential Analysis Module",
        add_help=True,  
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Basic parameters.
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="Number of threads to use (default: 4)")
    parser.add_argument("-g", "--genome", type=str, required=False,
                        help="Path to the genome FASTA file")
    parser.add_argument("-o", "--out_dir", type=str, required=True,
                        help="Output directory for results (default: project_name)")
    # parser.add_argument("-s1", "--sample1_outdir", type=str, required=True,
    #                     help="Directory for sample 1 data")
    # parser.add_argument("-s2", "--sample2_outdir", type=str, required=True,
    #                     help="Directory for sample 2 data")
    parser.add_argument("--prep", required=True, help="Path to TExTra prep output")
    parser.add_argument("--quant", required=True, help="Path to TExTra quant output")
    
    # Analysis-specific parameters.
    parser.add_argument("--table-only", action="store_true", default=False,
                        help="Output only the differential analysis table (default: False)")
    parser.add_argument("--project", type=str, default="sample1_vs_sample2",
                        help="Project name (default: sample1_vs_sample2)")
    parser.add_argument("--DE", type=bool, default=True,
                        help="Run differential expression analysis (default: True)")
    parser.add_argument("--ncpred", action="store_true",
                        help="Run ncPred analysis (default: False)")
    parser.add_argument("--label_no", type=int, default=5,
                        help="Label number for volcano plot (default: 5)")
    parser.add_argument("--minlength", type=int, default=200,
                        help="Discard transcripts shorter than this length (default: 200)")
    parser.add_argument("--log2fc", type=float, default=1,
                        help="Log2 fold-change threshold (default: 1)")
    parser.add_argument("--padj", type=float, default=0.05,
                        help="Adjusted p-value threshold (default: 0.05)")
    
    parser.add_argument("-s", "--samples", required=True,
                        help="Group conditions for differential analysis, separated by commas (default: False)")
    parser.add_argument("-m", "--model", type=str, default="ve",
                        help="Model for coding prediction. Choose from the following:\n"
                            "'ve' (default): vertebrate model, optimized for predicting coding potential in vertebrate species.\n"
                            "'pl': plant model, optimized for predicting coding potential in plant species.\n"
                            "Select the appropriate model based on the organism type you are working with.")
    parser.add_argument("--keep-temp", action="store_true", help="Keep intermediate files (default: False).")
    
    return parser.parse_args(args_list)

def main(args_list=None):
    """
    Main function to run the TE_exon analysis tool.
    
    This function parses command-line arguments (with defaults defined above) and calls:
      - DE_func for differential expression analysis if enabled.
      - ncPred_func for ncPred analysis if both DE and ncPred analyses are enabled.
    """
    args_list = sys.argv[2:]
    args = parse_arguments(args_list)

    if args.out_dir is None:
        args.out_dir = args.project

    # Run differential expression analysis if enabled.
    if args.DE:
        DE_func(args)
    
    # Run ncPred analysis if both DE and ncPred analyses are enabled.
    if args.DE and args.ncpred:
        ncPred_func(args)

if __name__ == "__main__":
    main()

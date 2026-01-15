#!/usr/bin/env python
import sys
import os
import argparse

from TExTra.bin.mode0_pipeline import main as prep_main
from TExTra.bin.mode1_pipeline import main as qual_main
from TExTra.bin.mode2_pipeline import main as quant_main
from TExTra.bin.mode3_pipeline import main as diff_main

def main():
    # Create main parser for global options and subcommand routing
    main_parser = argparse.ArgumentParser(
        prog="TExTra",
        description="TExTra Analysis Pipeline",
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter
    )
    main_parser.add_argument('-v', '--version', action='version', version='TExTra 1.1.0')
    
    # Add subparsers for available commands
    subparsers = main_parser.add_subparsers(
        title="Available Commands",
        dest="command",
        metavar="<command>",
        required=True 
    )
    
    # Subcommand: prep
    prep_parser = subparsers.add_parser(
        "prep",
        add_help=False, 
        help="Read Mapping and Transcriptome Assembly"
    )

    # Subcommand: qual
    qual_parser = subparsers.add_parser(
        "qual",
        add_help=False, 
        help="TE-derived Exon Identification and Classification"
    )
    
    # Subcommand: quant
    quant_parser = subparsers.add_parser(
        "quant",
        add_help=False,
        help="Quantification of TE-derived Exons"
    )
    
    # Subcommand: diff
    diff_parser = subparsers.add_parser(
        "diff",
        add_help=False,
        help="Downstream Analysis of TE-derived Exons"
    )
    
    # Parse global arguments only (stop at subcommand)
    main_args, remaining_args = main_parser.parse_known_args()
    
    # Dispatch to the corresponding module based on the subcommand
    if main_args.command == "prep":
        prep_main()
    elif main_args.command == "qual":
        qual_main()
    elif main_args.command == "quant":
        quant_main()
    elif main_args.command == "diff":
        diff_main()

if __name__ == "__main__":
    main()

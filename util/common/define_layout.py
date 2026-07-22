"""Define project-wide output layout constants."""

import os

CONVERT_DIR = "00_convert"
ALIGNMENT_DIR = "01_alignment"
ASSEMBLY_DIR = "02_assembly"
CLASSIFY_DIR = "03_classify"
QUANT_DIR = "04_quantification"
DOWNSTREAM_DIR = "05_downstream"

CONSENSUS_GTF = "consensus_transcripts.gtf"
TRANSCRIPT_GENE_ASSIGNMENT_TSV = "transcript_gene_assignment.tsv"


def resolve_consensus_gtf(prep_dir):
    """Return the expected consensus transcript GTF path."""
    return os.path.join(prep_dir, ASSEMBLY_DIR, CONSENSUS_GTF)

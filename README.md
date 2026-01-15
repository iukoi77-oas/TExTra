<p align="center">
  <img src="image/logo.png" width="60%" alt="TExTra Logo">
</p>

# TExTra

**TExTra** is a computational framework that integrates splice-junction usage modelling with expectation–maximization–based multi-mapping read redistribution to enable accurate positional classification and locus-specific quantification of TE-derived exons. TExTra robustly distinguishes TE-initiated, internal, and TE-terminated exons, restores quantitative accuracy at highly repetitive loci, and consistently outperforms existing short-read approaches in both sensitivity and locus specificity.

## Overview

The workflow integrates four coordinated modules: 
1. **Consensus reconstruction** of TE-chimeric transcripts.
2. **High-confidence identification** and positional classification of TE-derived exons.
3. **Locus-specific quantification** through reallocation of multi-mapping reads.
4. **Downstream functional interpretation**. 

This design enables recovery of complete exon boundaries and accurate quantification of exon usage across samples, overcoming limitations of existing approaches that detect only TE-gene junction fragments.

![Workflow Overview](image/overview.png)

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/iukoi77-oas/TExTra.git
cd TExTra
```

### 2. Create and activate the conda environment

```bash
conda env create -n TExTra -f environment.yml
conda activate TExTra
```

Then install TExTra in editable mode:

```bash
pip install -e .
```

### 3. External Tools Setup (TACO and PLEK2)

TExTra expects external tools to be placed under the `util/` directory. These tools are provided as compressed packages in [the Zenodo archive](https://zenodo.org/record/18252614).

```bash
# Extract the archives
tar -xzvf taco-v0.7.3.Linux_x86_64.tar.gz
tar -xzvf PLEKv2_allfiles_240807.tar.gz

# Decompress the PLEK2 model files
bunzip2 Coding_Net_kmer6_orf_Arabidopsis.h5.bz2
bunzip2 Coding_Net_kmer6_orf.h5.bz2
```

TExTra expects external tools to be placed under the `util/` directory:

```text
TExTra/
├── util/
│   ├── taco*/          # TACO binary directory
│   └── PLEK*/            # PLEK2 directory
```

> **Note**: The TACO and PLEK2 tools are already provided as compressed packages under the util/ directory.
They can be used directly by extracting the corresponding archives in place, and no additional download or reinstallation is required.

* Ensure that **TACO** and **PLEK2** are installed and executable.
* The `taco_run` executable should be located under `util/taco*/`.
* PLEK2.py should be located under `util/PLEK*/`.

If installed elsewhere, either move them into the corresponding directories or update the paths in the TExTra configuration accordingly.

### 4. Verify installation

After installation, you should be able to run:

```bash
TExTra --help
```

and see the available subcommands:

```text
prep     Read mapping and transcriptome assembly
qual     TE-derived exon identification and classification
quant    Quantification of TE-derived exons
diff     Downstream analysis
```
### 5. Common Issues & Solutions

* **Conda Priority**: If environment creation fails, run `conda config --set channel_priority flexible`.
* **PLEK2 Dependencies**: If import errors occur, manually run:
`pip install keras==2.4.3 tensorflow==2.4.1 regex bio`


## Tutorials & Usage

### 1. Preparing Input Data

The `input.csv` (or `.tsv`) file (tab-separated) should contain sample names in the first column, followed by paths to BAM or FASTQ files for biological replicates.

* For paired-end FASTQ, use a comma to separate mates (e.g., `R1.fq,R2.fq`).
* Ensure the number of columns is consistent across all rows (use empty spaces for samples with fewer replicates).

**Example `input.csv`:**

```text
heart_14d    ENCFF014ZHD.bam    ENCFF890XGT.bam
heart_2m     ENCFF997DFX.bam    ENCFF100LXH.bam
heart_20m    ENCFF402AFJ.bam    ENCFF931AJZ.bam
```
### 2. Pipeline Commands

#### **(1) Module 1: prep (Mapping & Assembly)**

Processes raw reads and performs transcriptome assembly.

```bash
TExTra prep -i input.tsv -g genome.fa -G gene.gtf -r TE.gtf -o ./output --strand rf
```

* `-r / --te`: TE annotation (supports BED, TXT, GTF, or `.out`).
* `--taco-disable`: Use Stringtie instead of TACO for merging.
* `--best`: Use optimal parameters for merging assemblies.

#### **(2) Module 2: qual (Classification)**

Identifies and classifies TE-derived exons based on their relationship with genes.

```bash
TExTra qual -s sample1,sample2 --prep ./output/prep -o ./output/qual
```

* `--ss3buffer` / `--ss5buffer`: Splice site buffer sizes (Default: 20/50).
* `--threshold`: Threshold for candidate selection (Default: 0.75).

*Note: The `-s` sample names must match the first column of your input.csv.*

#### **(3) Module 3: quant (Quantification)**

Locus-specific quantification using EM-based read redistribution.

```bash
TExTra quant -s sample1,sample2 -r TE.gtf --prep ./output/prep --qual ./output/qual -o ./output/quant
```

* `-e / --EM`: Number of EM iterations (Default: auto).
* `--bw`: Generate BigWig files for visualization (Requires `--normlib RPM`).

#### **(4) Module 4: diff (Differential Analysis)**

Performs differential expression and functional analysis.

```bash
TExTra diff -s group1,group2 --prep ./output/prep --quant ./output/quant -o ./output/diff -m ve
```

* `-m / --model`: Coding prediction model. Use `ve` for vertebrates or `pl` for plants.
* `--ncpred`: Enable ncPred analysis to predict non-coding potential.
* `--log2fc` / `--padj`: Significance thresholds (Default: 1.0 / 0.05).

## Test

[Test datasets](https://zenodo.org/record/18252614) are available on **Zenodo**. After downloading, decompress to find:

* **BAMs**: `ENCFF*.bam` (Sample replicates)
* **Genome**: `GRCm38.primary_assembly.genome.fa`
* **Annotations**: `gencode.vM21.primary_assembly.annotation.gtf` and `GRCm38_GENCODE_rmsk_TE.gtf`

The test script is located in the `test/` directory. It will automatically process the provided sample BAMs and annotations.

## Troubleshooting

* **Chromosome Mismatch**: Ensure that the chromosome naming convention (e.g., "chr1" vs "1") is consistent between your Genome FASTA, Gene GTF, and TE GTF.
* **Input Formatting**: Ensure the `input.csv` has a consistent number of columns.
* **Sample Name Consistency**: Ensure sample names used in the `-s` flag match the metadata provided in the `input.csv`.

## Citation

TExTra is currently in preparation. If you use this software in your research, please cite the following:

> **TExTra: A computational framework for locus-specific quantification and classification of TE-derived exons.** > *Yanjie et al. (2024). In preparation.*

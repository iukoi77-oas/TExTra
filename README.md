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

TExTra is designed for Unix-like environments and has been tested on **CentOS Linux 7.9 (x86_64)**. All software dependencies are fully specified in the provided `environment.yml` file.

Installation takes approximately 10 minutes, although the time required for downloading dependencies may vary depending on internet connection speed.

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

Overall runtime is around 35 minutes, largely dominated by TExTra qual, which performs extensive resampling with 100,000 pseudo-replicates and 1,000 bootstrap iterations.

## Expected Output

TExTra generates structured, module-specific outputs. Below we describe the **core result files** produced by each module and their biological interpretation.

### 1. TExTra prep — Transcriptome Assembly

**Core output file:** `assembly/novel_transcripts.gtf`

```text
chr18	taco	transcript	4196962	4199260	.	+	.	transcript_id "TU6"; gene_id "ENSMUSG00000073647.4"; xloc "XLOC_000004"; cmp_ref "ENSMUST00000097690.4"; class_code "x"; cmp_ref_gene "Gm10557"; tss_id "TSS4";
chr18	taco	exon	4196962	4199260	.	+	.	transcript_id "TU6"; gene_id "ENSMUSG00000073647.4"; exon_number "1";
chr18	taco	transcript	4382398	4383185	.	+	.	transcript_id "TU7"; gene_id "ENSMUSG00000024234.7"; xloc "XLOC_000005"; cmp_ref "ENSMUST00000234186.1"; class_code "i"; cmp_ref_gene "Mtpap"; tss_id "TSS5";
chr18	taco	exon	4382398	4382424	.	+	.	transcript_id "TU7"; gene_id "ENSMUSG00000024234.7"; exon_number "1";
```

**Description:**

* This file contains the **consensus transcript annotation** reconstructed from **all samples and all biological replicates**.
* Transcripts are assembled by integrating splice junction evidence across samples.
* The resulting GTF represents a unified transcriptome reference used by all downstream modules (`qual`, `quant`, `diff`).

### 2. TExTra qual — Identification of TE-derived Exons

**Core output directory:**

```text
candidate/
├── sample_1/
│   └── exon_TE.bed
├── sample_2/
│   └── exon_TE.bed
```

**Key file:** `candidate/<sample_name>/exon_TE.bed`

```text
chr18	34631682	34634164	ENSMUSG00000003779.16:TU394:exon_1:first	1.0	+	chr18	34633521	34633608	chr18|34633521|34633608|PB1D10:Alu:SINE|307|+	307	+
chr18	36666680	36667659	ENSMUSG00000006050.12:TU468:exon_1:first	1.0	+	chr18	36666524	36666691	chr18|36666524|36666691|B2_Mm2:B2:SINE|835|+	835	+
chr18	68321317	68324845	ENSMUSG00000009535.13:TU950:exon_11:first	1.0	-	chr18	68324746	68324846	chr18|68324746|68324846|L1ME3B:L1:LINE|312|-	312	-
chr18	3511981	3516404	ENSMUSG00000024232.2:TU4:exon_4:first	1.0	-	chr18	3513877	3513957	chr18|3513877|3513957|B1_Mur2:Alu:SINE|257|-	257	-
```
**Column description:**

| Column | Name                    | Description                                                                              |
| -----: | ----------------------- | ---------------------------------------------------------------------------------------- |
|      1 | exon_chr                | Chromosome of the exon                                                                   |
|      2 | exon_start              | 0-based start coordinate of the exon                                                     |
|      3 | exon_end                | End coordinate of the exon                                                               |
|      4 | exon_id                 | Composite exon identifier                                                                |
|      5 | replicate_support_ratio | Fraction of biological replicates in which the TE-derived exon is detected               |
|      6 | exon_strand             | Strand of the exon                                                                       |
|      7 | TE_chr                  | Chromosome of the transposable element                                                   |
|      8 | TE_start                | Start coordinate of the TE insertion                                                     |
|      9 | TE_end                  | End coordinate of the TE insertion                                                       |
|     10 | TE_id                   | Composite TE identifier encoding genomic location, TE classification, length, and strand |
|     11 | TE_length               | Length of the TE insertion (bp)                                                          |
|     12 | TE_strand               | Strand of the TE                                                                         |

* **exon_id**: Composite exon identifier (Column 4) in the format `gene_id:transcript_id:exon_number:exon_type`
* **TE_id**: Composite TE identifier (Column 10) format: `chr|start|end|TE_subfamily:TE_family:TE_class|length|strand`
* **exon_type**: indicates the positional relationship between the TE-derived exon and the host transcript: `first` -> TE-initiated exon, `internal` -> TE-internal exon, `last` -> TE-terminated exon.

### 3. TExTra quant — Locus-specific Quantification

#### (1) TE–Transcript Mapping

**Core file:** `quantification/project_TEexon.bed`

```text
chr18	3511981	3516404	ENSMUSG00000024232.2:TU4:exon_4:first	1.0	-	chr18|3513877|3513957|B1_Mur2:Alu:SINE|257|-
chr18	6788479	6794429	ENSMUSG00000073639.6:TU56:exon_7:first	1.0	-	chr18|6792451|6792700|B4A:B4:SINE|731|-
chr18	6788479	6794429	ENSMUSG00000073639.6:TU56:exon_7:first	1.0	-	chr18|6790590|6790661|RLTR20A4:ERVK:LTR|267|-
chr18	6788479	6794429	ENSMUSG00000073639.6:TU56:exon_7:first	1.0	-	chr18|6793898|6793973|L3:CR1:LINE|182|-
```

**Column description:**

| Column | Name                 | Description                                                                              |
| -----: | -------------------- | -----------------------------------------------------------------------------------------|
|      1 | metaexon_chr         | Chromosome of the metaexon                                                               |
|      2 | metaexon_start       | 0-based start coordinate of the metaexon                                                 |
|      3 | metaexon_end         | End coordinate of the metaexon                                                           |
|      4 | exon_id              | Composite exon identifier linking gene, transcript, exon number, and exon type           |
|      5 | sample_support_ratio | Proportion of samples in which the TE-derived metaexon is detected                       |
|      6 | metaexon_strand      | Strand of the metaexon                                                                   |
|      7 | TE_id                | Composite TE identifier encoding genomic location, TE classification, length, and strand |

#### (2) Quantitative Expression Matrix

**Core file:** `quantification/project_TEexons.txt`

```text
metaexon	heart_postnatal_14d_rep1	heart_postnatal_14d_rep2	heart_adult_2m_rep2	heart_adult_2m_rep1
chr18:76304607-76311681:-	1217.0	1685.0	955.0	1460.0
chr18:76944092-76944531:-	102.0	169.0	87.5	99.5
chr18:76954978-76955186:-	741.0	919.0	719.5	862.5
chr18:76970611-76975307:-	10067.0	10088.0	7782.0	9737.0
```

**Description:**

* A matrix of **TE-derived exon expression values**.
* Rows: metaexons (genomic coordinates)
* Columns: each biological replicate of each sample
* Values: EM-corrected expression estimates


### 4. TExTra diff — Differential Analysis & Functional Annotation

#### (1) Differential Expression of TE-derived Exons

**Core file:** `DE/DESeq2_all.txt`

```text
coord	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	genes	TE_info
chr18:32528322-32530201:+	5639.7595556878	1.25648720324594	0.222127995490816	5.65659092393822	1.5440919540921e-08	7.25723218423285e-07	ENSMUSG00000090523.2	chr18|32528245|32528502|MER44A:TcMar-Tigger:DNA|692|+
chr18:37640491-37643532:+	2118.36948661845	1.15307733133441	0.241734237426584	4.77002076168298	1.84206932438774e-06	4.3288629123112e-05	ENSMUSG00000051316.8	chr18|37641009|37641107|PB1D7:Alu:SINE|272|+
chr18:12990493-12992946:-	980.914529681156	-1.16959737465939	0.263338124634987	-4.44142820672309	8.93637627899472e-06	0.0001400032283709	ENSMUSG00000024423.6	chr18|12992799|12992994|L1ME4a:L1:LINE|449|-
chr18:9876988-9882643:-	1727.55716758539	0.943080751628008	0.245639949098282	3.83928084617326	0.0001233952173325	0.0014498938036579	ENSMUSG00000036103.9	chr18|9877336|9877731|MTA_Mm:ERVL-MaLR:LTR|3728|-,chr18|9878899|9879302|ERVB7_3-LTR_MM:ERVK:LTR|3114|-
```

**Description:**

* Differential expression results for TE-derived exons.
* Includes:

  * log2 fold change
  * statistical significance (p-value, adjusted p-value)
  * associated genes
  * TE annotation information

#### (2) Coding Potential Prediction of Differential Transcripts

**Core output files:** `ncPred/plek_final_result.csv` and `ncPred/sig_transcripts.gtf`

```text
Transcript,Prediction,Gene,DiffExon,MetaExon,TE
TU192,0,ENSMUSG00000024423.6,exon_11:first,chr18:12990493-12992946:-,chr18|12992799|12992994|L1ME4a:L1:LINE|449|-
TU350,0,ENSMUSG00000090523.2,exon_1:first,chr18:32528322-32530201:+,chr18|32528245|32528502|MER44A:TcMar-Tigger:DNA|692|+
TU492,0,ENSMUSG00000051316.8,exon_1:first,chr18:37640491-37643532:+,chr18|37641009|37641107|PB1D7:Alu:SINE|272|+
```

**Descriptions:**

* `plek_final_result.csv`:

  * Coding potential prediction per transcript
  * `0` = non-coding
  * `1` = coding potential detected
* `sig_transcripts.gtf`:

  * GTF annotation of transcripts containing **significant differential TE-derived exons**
  * Can be directly visualized in genome browsers

## Troubleshooting

* **Chromosome Mismatch**: Ensure that the chromosome naming convention (e.g., "chr1" vs "1") is consistent between your Genome FASTA, Gene GTF, and TE GTF.
* **Input Formatting**: Ensure the `input.csv` has a consistent number of columns.
* **Sample Name Consistency**: Ensure sample names used in the `-s` flag match the metadata provided in the `input.csv`.

## Citation

TExTra is currently in preparation. If you use this software in your research, please cite the following:

> **TExTra enables locus-specific quantification of transposable element–derived exonization from short-read RNA sequencing.**
> *Yan et al., under review.*

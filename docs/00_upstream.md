# 00 - TExTra upstream

`TExTra upstream` runs `prep + qual + quant` in one command. It is the recommended entry point for regular analysis when differential testing is not yet needed.

## Minimal Command

```bash
TExTra upstream \
  --input input.tsv \
  --out_dir result \
  --genome genome.fa \
  --gene gene.gtf \
  --te TE_annotation
```

## Required Arguments

| Argument | Meaning |
| --- | --- |
| `-i, --input` | Input TSV used by prep and quant read resolution. |
| `-o, --out_dir` | Output root directory shared by prep, qual, and quant. |
| `-g, --genome` | Genome FASTA used by prep and quant. |
| `-G, --gene` | Reference gene annotation GTF used by prep. |
| `-r, --te` | TE annotation file used by prep. Supported formats: `.gtf`, `.gff`, `.bed`, RepeatMasker `.out`, and RepeatMasker `.txt`. |

## Key Options

| Argument | Meaning |
| --- | --- |
| `--samples` | Sample/condition names, comma-separated. If omitted, infer from the input TSV first column. |
| `--threads` | Total thread budget passed to child modules. |
| `--njobs` | Maximum number of parallel jobs. |
| `--project` | Project prefix used by qual and quant result tables. |
| `--strand` | Library strand mode. |
| `--readtype` | `paired` or `single`. |
| `--detail` | Propagate detail mode to prep, qual, and quant. |
| `--debug` | Propagate debug mode to prep, qual, and quant. |
| `--calculate-afe-ale` | Enable AFE/ALE usage outputs in qual. |
| `--skip-hitindex` | Skip HITindex positional classification in qual. |
| `--ignore-junction` | Ignore TE-overlap junction support/degradation checks in qual. |
| `--quantifier` | Quantification backend for quant: `rsem` or `salmon`. |
| `--quant-result-dir` | Reuse existing RSEM/Salmon backend outputs. |
| `--compute-gene-abundance` | Export gene-level abundance in quant. |

## Terminal Output

`upstream` separates child module output with numbered banners:

```text
==================== TExTra 1: prep ====================
==================== TExTra 2: qual ====================
==================== TExTra 3: quant ====================
```

Each child module writes the same outputs as its standalone command. See:

* [01 - TExTra prep](01_prep.md)
* [02 - TExTra qual](02_qual.md)
* [03 - TExTra quant](03_quant.md)

## Main Outputs

Default upstream results include:

```text
00_convert/TE_anno.bed
00_convert/gene_anno.bed
02_assembly/consensus_transcripts.gtf
02_assembly/transcript_gene_assignment.tsv
03_classify/transcript_exon_te_annotation.tsv
03_classify/<project>.TE_overlap.exon.tsv
04_quantification/<project>.TE_overlap.exon_usage.tsv
04_quantification/transcript_abundance.tsv
TExTra.config.json
```

`TExTra.config.json` records the shared run context and is used by later modules with `--reuse`.

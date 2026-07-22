# 01 - TExTra prep

`TExTra prep` maps reads, assembles per-replicate transcripts, merges assemblies, assigns transcript gene IDs, and converts TE/gene annotations.

## Minimal Command

```bash
TExTra prep \
  --input input.tsv \
  --out_dir result \
  --genome genome.fa \
  --gene gene.gtf \
  --te TE_annotation
```

## Required Arguments

| Argument | Meaning |
| --- | --- |
| `-i, --input` | Input TSV. First column is sample/condition; following columns are replicate inputs. |
| `-o, --out_dir` | Prep output root directory. |
| `-g, --genome` | Genome FASTA used for STAR index generation and annotation conversion. |
| `-G, --gene` | Reference gene annotation GTF. |
| `-r, --te` | TE annotation file converted to TE BED. Supported formats: `.gtf`, `.gff`, `.bed`, RepeatMasker `.out`, and RepeatMasker `.txt`. |

## Key Options

| Argument | Meaning |
| --- | --- |
| `--assembly-mode {de-novo,reference-guided}` | Per-replicate StringTie assembly mode. |
| `--merge-method {taco,stringtie}` | Transcript merge backend. |
| `--optimized` | Use optimized merge filtering defaults unless explicit `--min-*` values are provided. |
| `--junction` | StringTie first-pass minimum junction coverage. The same value is stored in `TExTra.config.json` and used by downstream qual as the TE-side junction degrade read threshold. |
| `--min-expr` | Merge minimum expression filter. |
| `--min-length` | Merge minimum transcript length filter. |
| `--min-frac` | Merge minimum isoform-fraction filter. |
| `--index` | Existing STAR index directory. |
| `--taco-path` | Path to `taco_run` or a TACO installation directory. |
| `--detail` | Add result-checking detail tables and summaries. |
| `--debug` | Keep selected intermediate files, debug config, and extended command logs. |

Advanced STAR-native options:

* `--sjdbOverhang`
* `--seedMultimapNmax`
* `--winAnchorMultimapNmax`
* `--outFilterMultimapNmax`

## Default Outputs

```text
00_convert/TE_anno.bed
00_convert/gene_anno.bed
02_assembly/consensus_transcripts.gtf
02_assembly/transcript_gene_assignment.tsv
TExTra.config.json
```

### transcript_gene_assignment.tsv

| Column | Meaning |
| --- | --- |
| `transcript_id` | Assembled transcript ID. |
| `original_gene_id` | Gene ID assigned by the assembly/merge step before TExTra reference-gene reassignment. |
| `assigned_gene_id` | Final gene ID used by downstream modules. |
| `assigned_gene_name` | Reference gene name when the transcript is assigned to a reference gene with a known name. |
| `unassigned_gene_type` | Retained unassigned category when no confident reference-gene assignment is made. |

## Detail Outputs

Detail mode keeps default outputs and adds:

```text
02_assembly/transcript_gene_assignment.detail.tsv
```

### transcript_gene_assignment.detail.tsv

| Column | Meaning |
| --- | --- |
| `transcript_id` | Assembled transcript ID. |
| `original_gene_id` | Gene ID before TExTra reassignment. |
| `assigned_gene_id` | Final downstream gene ID. |
| `assigned_gene_name` | Reference gene name for clear reference-gene assignments. |
| `reference_gene_id` | Reference gene ID reported by the reference comparison step. |
| `reference_gene_name` | Reference gene name associated with `reference_gene_id`. |
| `reference_transcript_id` | Reference transcript ID reported by the reference comparison step. |
| `class_code` | gffcompare class code used as assignment/filtering evidence. |
| `strand_consistent` | Whether transcript/reference overlap is strand-consistent. |
| `exon_count` | Number of exons in the assembled transcript. |
| `unassigned_gene_type` | Retained unassigned category when no confident reference-gene assignment is made. |

## Debug Outputs

Debug mode keeps detail outputs and additionally preserves:

```text
01_alignment/ STAR original outputs, samtools intermediates, and STAR index files
02_assembly/tmp/
logs/prep_config.json
logs/prep_extend.log
```

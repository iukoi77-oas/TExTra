# 03 - TExTra quant

`TExTra quant` estimates transcript abundance with RSEM or Salmon and projects transcript abundance to TE-overlap exon usage.

## Minimal Command

```bash
TExTra quant \
  --prep result \
  --qual result \
  --out_dir result \
  --samples heart_14d,heart_2m \
  --genome genome.fa
```

If `TExTra.config.json` from previous modules is available, the shortest form is:

```bash
TExTra quant --reuse --out_dir result --genome genome.fa
```

## Key Arguments

| Argument | Meaning |
| --- | --- |
| `--prep` | Path to prep output; required unless inferred by `--reuse`. |
| `--qual` | Path to qual output; required unless inferred by `--reuse`. |
| `-o, --out_dir` | Quant output root directory. |
| `-s, --samples` | Sample/condition names, comma-separated; required unless inferred by `--reuse`. |
| `-i, --input` | Original prep input TSV used to resolve quant reads. If omitted, quant may use prep BAM fallback. |
| `-g, --genome` | Genome FASTA used to build RSEM reference or Salmon transcript FASTA when needed. |
| `--reuse` | Infer omitted prep/qual context. Does not reuse quant result files. |

## Key Options

| Argument | Meaning |
| --- | --- |
| `--quantifier {rsem,salmon}` | Quantification backend. |
| `--quant-result-dir` | Reuse RSEM/Salmon backend outputs, typically from a debug quant run. Final usage tables are regenerated. |
| `--compute-gene-abundance` | Also write gene-level abundance. |
| `--detail` | Add result-checking detail tables and summaries. |
| `--debug` | Keep organized backend artifacts, debug config, and extended command logs. |

## Default Outputs

```text
04_quantification/<project>.TE_overlap.exon_usage.tsv
04_quantification/transcript_abundance.tsv
04_quantification/gene_abundance.tsv    # only with --compute-gene-abundance
```

### <project>.TE_overlap.exon_usage.tsv

| Column | Meaning |
| --- | --- |
| `exon_id` | Coordinate-style exon ID: `chr:start-end:strand`. |
| `<sample usage columns...>` | TE-overlap exon usage value for each quantified sample/replicate. |
| `metaexon_id` | HITindex/ID-position mapping source from qual. |
| `gene_id` | Gene ID associated with the exon event. |
| `transcript_id` | Transcript IDs supporting the exon event. |
| `te_overlap_label` | Final TE-overlap status from qual. |
| `te_boundary_side` | TE boundary side in transcript/exon direction: `5'`, `3'`, `both`, or `none`. |
| `te_splice_site_repeat_TE` | TE annotation supporting splice-site repeat/TE-overlap interpretation. |
| `ID_position_summary` | Summary of inferred TE-exon position class. |
| `ID_position_confidence` | Confidence label for `ID_position_summary`. |
| `candidate_TE_event` | Candidate TE event label derived from position and TE-overlap evidence. |
| `candidate_TE_confidence` | Confidence label for the candidate TE event. |

Detail/debug additional columns:

| Column | Meaning |
| --- | --- |
| `te_overlap_bp_max` | Maximum base-pair overlap between exon and TE annotations. |
| `te_overlap_frac_max` | Maximum fractional overlap between exon and TE annotations. |
| `junction_te_side_reads_max` | Maximum TE-side junction read count. |
| `support_sample_n` | Number of samples with usage `>= 0.1`. |
| `support_sample_ratio` | `support_sample_n` divided by total sample count. |
| `ID_position_source` | Evidence source used for ID-position assignment. |
| `transcript_structure_roles` | Position roles inferred from transcript structure. |
| `HITindex_structure_roles` | Position roles inferred from HITindex results. |
| `HITindex_evaluable_sample_n` | Number of samples evaluable by HITindex. |
| `HITindex_evaluable_replicate_n` | Number of replicates evaluable by HITindex. |

### transcript_abundance.tsv

| Column | Meaning |
| --- | --- |
| `sample` | Quantified sample or replicate ID. |
| `transcript_id` | Transcript ID from the consensus transcriptome. |
| `gene_id` | Gene ID associated with the transcript. |
| `te_overlap_exon_transcript` | `yes` when the transcript supports at least one TE-overlap exon event; otherwise `no`. This is not a transcript biotype. |
| `estimated_count` | Backend-normalized estimated read/count value. |
| `tpm` | Transcript TPM. |

### gene_abundance.tsv

Written only with `--compute-gene-abundance`.

| Column | Meaning |
| --- | --- |
| `sample` | Quantified sample or replicate ID. |
| `gene_id` | Gene ID. |
| `estimated_count` | Backend-normalized estimated gene count. |
| `tpm` | Gene TPM. |

## Debug Outputs

Debug mode may additionally retain organized backend artifacts:

```text
04_quantification/RSEM/index/
04_quantification/RSEM/result/
04_quantification/RSEM/bam2fastq/       # only when BAM fallback conversion is used

04_quantification/salmon/index/
04_quantification/salmon/result/
04_quantification/salmon/bam2fastq/     # only when BAM fallback conversion is used

logs/quant_config.json
logs/quant_extend.log
```

Backend command logs are aggregated into `logs/quant_extend.log`, not duplicated under backend result directories.

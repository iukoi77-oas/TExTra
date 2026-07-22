# 04 - TExTra diff

`TExTra diff` tests differential TE-overlap exon usage and can optionally run ncPred/PLEK2 coding-potential prediction.

## Minimal Command

```bash
TExTra diff \
  --prep result \
  --quant result \
  --out_dir result \
  --groups heart_14d,heart_2m
```

If `TExTra.config.json` from previous modules is available and exactly two groups can be inferred, the shortest form is:

```bash
TExTra diff --reuse --out_dir result
```

## Key Arguments

| Argument | Meaning |
| --- | --- |
| `--prep` | Path to prep output; required unless inferred by `--reuse`. |
| `--quant` | Path to quant output; required unless inferred by `--reuse`. |
| `-o, --out_dir` | Diff output root directory. |
| `--groups` | Two condition names to compare, comma-separated. |
| `--reuse` | Infer omitted prep/quant/project/groups/genome context. Does not reuse diff result files. |

## Key Options

| Argument | Meaning |
| --- | --- |
| `--test-method {classical,empirical}` | Differential usage test method. |
| `--delta-exon-usage` | Minimum absolute exon-usage difference. Default: `0.1`. |
| `--padj` | Adjusted p-value threshold. Default: `0.05`. |
| `--pvalue` | Optional raw p-value threshold. |
| `--paired` | Use paired replicate testing for classical mode. |
| `--empirical-background-n` | Local background size for empirical mode. |
| `--ncpred` | Run ncPred coding-potential prediction. |
| `-g, --genome` | Genome FASTA; required when `--ncpred` is enabled. |
| `--plek-path` | Path to `PLEK2.py` or a PLEK installation directory. If omitted, TExTra searches `TEXTRA_EXTERNAL_DIR/PLEK*/PLEK2.py` and `util/external/PLEK*/PLEK2.py`. |
| `--ncpred-model {ve,pl}` | PLEK2 model, vertebrate or plant. |
| `--min-length` | Minimum transcript length retained for ncPred candidates. |
| `--detail` | Add all-event differential usage table and result-checking summaries. |
| `--debug` | Keep debug config, extended logs, and selected ncPred FASTA when applicable. |

By default, significant events satisfy:

```text
abs(delta_usage) > --delta-exon-usage and padj < --padj
```

## Default Outputs

```text
05_downstream/DE/differential_significant_usage.tsv
05_downstream/ncPred/plek_result.csv             # only with --ncpred
05_downstream/ncPred/selected_transcripts.gtf    # only with --ncpred
```

### differential_significant_usage.tsv

| Column | Meaning |
| --- | --- |
| `exon_id` | Coordinate-style exon ID tested for differential usage. |
| `group1` | First condition in `--groups`. |
| `group2` | Second condition in `--groups`. |
| `mean_usage_group1` | Mean exon usage in `group1`. |
| `mean_usage_group2` | Mean exon usage in `group2`. |
| `delta_usage` | `mean_usage_group1 - mean_usage_group2`. |
| `higher_usage_group` | Group with higher mean exon usage. |
| `pvalue` | Raw p-value from the selected differential usage test. |
| `padj` | Multiple-testing adjusted p-value. |
| `metaexon_id` | HITindex/ID-position mapping source from upstream output, when present. |
| `gene_id` | Gene ID associated with the exon event. |
| `transcript_id` | Transcript IDs supporting the exon event. |
| `te_overlap_label` | Final TE-overlap status from qual. |
| `te_boundary_side` | TE boundary side in transcript/exon direction. |
| `te_splice_site_repeat_TE` | TE annotation supporting splice-site repeat/TE-overlap interpretation. |

Because this file contains only significant events, it does not include a redundant `significant_usage` column.

### ncPred default outputs

`05_downstream/ncPred/plek_result.csv`:

| Column | Meaning |
| --- | --- |
| `Transcript` | Transcript ID submitted to PLEK2. |
| `Prediction` | PLEK2 coding-potential label, where `0` indicates non-coding and `1` indicates coding potential. |
| `Gene` | Gene ID associated with the transcript. |
| `DiffExon` | Differential TE-overlap exon event associated with the transcript. |
| `MetaExon` | Metaexon/event coordinate associated with the differential event. |
| `TE` | TE annotation associated with the differential event. |

`05_downstream/ncPred/selected_transcripts.gtf` contains the GTF records for transcripts carrying significant TE-overlap exon events selected for ncPred.

## Detail/Debug Outputs

Detail/debug add:

```text
05_downstream/DE/differential_usage.tsv
05_downstream/ncPred/selected_transcripts.detail.tsv    # only with --ncpred
```

### differential_usage.tsv

| Column | Meaning |
| --- | --- |
| `exon_id` | Coordinate-style exon ID tested for differential usage. |
| `<original sample usage columns...>` | Original exon usage values from quant for each sample/replicate. |
| `group1` | First condition in `--groups`. |
| `group2` | Second condition in `--groups`. |
| `mean_usage_group1` | Mean exon usage in `group1`. |
| `mean_usage_group2` | Mean exon usage in `group2`. |
| `delta_usage` | `mean_usage_group1 - mean_usage_group2`. |
| `higher_usage_group` | Group with higher mean exon usage. |
| `pvalue` | Raw p-value from the selected test. |
| `padj` | Multiple-testing adjusted p-value. |
| `test_method` | Differential usage test method used for the row. |
| `n_non_na_group1` | Number of non-missing usage values in `group1`. |
| `n_non_na_group2` | Number of non-missing usage values in `group2`. |
| `usage_pattern_flag` | Result-checking flag describing the usage-data pattern used in testing/filtering. |
| `metaexon_id` | HITindex/ID-position mapping source from upstream output, when present. |
| `gene_id` | Gene ID associated with the exon event. |
| `transcript_id` | Transcript IDs supporting the exon event. |
| `te_overlap_label` | Final TE-overlap status from qual. |
| `te_boundary_side` | TE boundary side in transcript/exon direction. |
| `te_splice_site_repeat_TE` | TE annotation supporting splice-site repeat/TE-overlap interpretation. |

Debug mode additionally retains:

```text
logs/diff_config.json
logs/diff_extend.log
05_downstream/ncPred/selected_transcripts.fasta    # only with --ncpred
```

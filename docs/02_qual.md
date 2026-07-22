# 02 - TExTra qual

`TExTra qual` identifies TE-overlap exons, optionally fits HITindex models, and writes exon-level positional evidence used by quantification.

## Minimal Command

```bash
TExTra qual \
  --prep result \
  --out_dir result \
  --samples heart_14d,heart_2m
```

If `TExTra.config.json` from `prep` is available, the shortest form is:

```bash
TExTra qual --reuse --out_dir result
```

## Key Arguments

| Argument | Meaning |
| --- | --- |
| `--prep` | Path to TExTra prep output; required unless inferred by `--reuse`. |
| `-o, --out_dir` | Qual output root directory. |
| `-s, --samples` | Sample/condition names, comma-separated; required unless inferred by `--reuse`. |
| `--reuse` | Infer omitted prep context. Does not reuse qual result files. |

## Key Options

| Argument | Meaning |
| --- | --- |
| `--calculate-afe-ale` | Calculate AFE/ALE outputs. Requires HITindex. |
| `--skip-hitindex` | Skip HITindex positional classification. Junction evidence is disabled. |
| `--ignore-junction` | Ignore junction support/degradation checks while still generating TE-overlap annotation. Cannot be combined with `--junction`. |
| `-j, --junction` | Minimum TE-side junction reads required to retain a TE-overlap call during junction evidence degradation. Default: `2.0`. If omitted and prep used a non-default `--junction`, inherit prep `--junction` and print an INFO message. Cannot be combined with `--ignore-junction` or `--skip-hitindex`. |
| `--hitindex-dir` | Reuse completed per-replicate HITindex outputs. |
| `--ss3-buffer` | 3' splice-site buffer size. |
| `--ss5-buffer` | 5' splice-site buffer size. |
| `--te-overlap-min-bp` | Minimum overlap bp for TE-overlap metrics. |
| `--te-overlap-min-frac` | Minimum overlap fraction for TE-overlap metrics. |
| `--splice-site-flank-bp` | Flank window around splice sites for boundary/anchor/var-site checks. |
| `--genmodel-iters` | ADVI iterations in HITindex generative model. |
| `--bootstrap-n` | Bootstrap iterations for HITindex significance. |
| `--seed` | Optional random seed. |
| `--detail` | Add result-checking detail tables and summaries. |
| `--debug` | Keep selected reusable HITindex outputs, debug config, and extended logs. |

Conflict rules:

* `--calculate-afe-ale` cannot be combined with `--skip-hitindex`.
* `--hitindex-dir` cannot be combined with `--skip-hitindex`.
* `--skip-hitindex` disables junction evidence.
* Explicit `--junction` cannot be combined with `--ignore-junction` or `--skip-hitindex`.
* If `--hitindex-dir` lacks valid junction support, qual warns and continues as if `--ignore-junction` was used.
* If `--hitindex-dir` contains complete AFE/ALE outputs, qual enables AFE/ALE output generation.

Junction evidence note:

* qual `--junction` defaults to `2.0`.
* If `--junction` is not provided and prep used a non-default `--junction`, qual reads prep's `junction` value from `<prep>/TExTra.config.json`, applies that value, and prints an INFO message.
* This value is reported as `junction_degrade_min_reads` and is used as the TE-side junction read threshold for TE-overlap degradation.
* If the prep config is unavailable or prep used the default value, qual uses `2.0`.

## Default Outputs

```text
03_classify/transcript_exon_te_annotation.tsv
03_classify/<project>.TE_overlap.exon.tsv
03_classify/<project>.TE_overlap.AFEPSI.tsv    # only with AFE/ALE
03_classify/<project>.TE_overlap.ALEPSI.tsv    # only with AFE/ALE
```

### transcript_exon_te_annotation.tsv

This is the transcript-exon interface table used by downstream modules. It is written directly under `03_classify/`.

| Column | Meaning |
| --- | --- |
| `exon_id` | Coordinate-style exon ID: `chr:start-end:strand`. |
| `metaexon_id` | HITindex metaexon/source event ID. Included when HITindex is not skipped. |
| `gene_id` | Gene ID associated with the transcript exon. |
| `transcript_id` | Transcript ID containing the exon. |
| `exon_number` | Exon order within the transcript. |
| `transcript_exon_count` | Total number of exons in the transcript. |
| `transcript_strand` | Strand of the parent transcript. |
| `exon_chrom` | Exon chromosome. |
| `exon_start` | Exon start coordinate. |
| `exon_end` | Exon end coordinate. |
| `exon_strand` | Exon strand. |
| `te_overlap_label` | Final TE-overlap status. If junction evidence is effective, this reflects junction-based degradation checks. |
| `te_splice_site_repeat_TE` | TE annotation overlapping or nearest to the splice-site repeat evidence used for TE-overlap interpretation. |

Detail/debug additional columns:

| Column | Meaning |
| --- | --- |
| `te_overlap_bp_max` | Maximum base-pair overlap between the exon and TE annotations. |
| `te_overlap_frac_max` | Maximum fractional overlap between the exon and TE annotations. |
| `te_boundary_hit_any` | Whether a TE boundary hit was detected at either splice-site side. |
| `te_boundary_side` | Transcript/exon direction of the TE boundary hit: `5'`, `3'`, `both`, or `none`. |
| `junction_support_sample_n` | Number of samples with junction support. Present only when junction evidence is effective. |
| `junction_support_sample_ids` | Sample IDs with junction support. Present only when junction evidence is effective. |
| `junction_te_side_reads_max` | Maximum TE-side junction read count. Present only when junction evidence is effective. |
| `junction_te_side_reads_mean_supported_samples` | Mean TE-side junction read count among supported samples. Present only when junction evidence is effective. |

### <project>.TE_overlap.exon.tsv

This is the default exon-level TE-overlap result table and the main qual-to-quant contract. It is one row per deduplicated `exon_id`.

| Column | Meaning |
| --- | --- |
| `exon_id` | Coordinate-style exon ID: `chr:start-end:strand`. |
| `metaexon_id` | HITindex/ID-position mapping source. This is not the TE-overlap decision unit. |
| `gene_id` | Gene ID associated with the exon event. |
| `transcript_id` | Transcript IDs supporting the exon event; comma-merged when multiple transcripts share the exon. |
| `transcript_exon_id` | Transcript-specific exon IDs; comma-merged when multiple transcripts share the exon. |
| `te_overlap_label` | Final TE-overlap status. |
| `te_splice_site_repeat_TE` | TE annotation supporting the splice-site repeat/TE-overlap interpretation. |
| `te_boundary_side` | TE boundary side in transcript/exon direction: `5'`, `3'`, `both`, or `none`. |
| `ID_position_summary` | Summary of inferred TE-exon position class from HITindex and/or transcript structure evidence. |
| `ID_position_confidence` | Confidence label for `ID_position_summary`. |
| `candidate_TE_event` | Candidate TE event label derived from position and TE-overlap evidence. |
| `candidate_TE_confidence` | Confidence label for the candidate TE event. |

Detail/debug additional columns in `03_classify/<project>.TE_overlap.exon.detail.tsv`:

| Column | Meaning |
| --- | --- |
| `te_overlap_bp_max` | Maximum base-pair overlap between exon and TE annotations. |
| `te_overlap_frac_max` | Maximum fractional overlap between exon and TE annotations. |
| `junction_te_side_reads_max` | Maximum TE-side junction read count. Present only when junction evidence is effective. |
| `ID_position_source` | Evidence source used for ID-position assignment. |
| `transcript_structure_roles` | Position roles inferred from transcript structure. |
| `HITindex_structure_roles` | Position roles inferred from HITindex results. |
| `HITindex_evaluable_sample_n` | Number of samples evaluable by HITindex. |
| `HITindex_evaluable_replicate_n` | Number of replicates evaluable by HITindex. |

`event_id` and `te_exon_id` are not output.

### AFE/ALE Outputs

Default AFE/ALE outputs are TE-overlap-only exon-level PSI tables. They are one row per deduplicated `exon_id`.

| Column | Meaning |
| --- | --- |
| `exon_id` | Coordinate-style exon ID: `chr:start-end:strand`. |
| `metaexon_id` | Original AFE/ALE metaexon event coordinate. |
| `gene_id` | Gene ID associated with the AFE/ALE event. |
| `transcript_id` | Transcript IDs supporting the event. |
| `transcript_exon_id` | Transcript-specific exon IDs supporting the event. |
| `te_overlap_label` | Final TE-overlap status. |
| `te_splice_site_repeat_TE` | TE annotation supporting the splice-site repeat/TE-overlap interpretation. |
| `<sample/replicate PSI columns...>` | AFE/ALE PSI values for each sample or replicate. |

Detail/debug additionally retain:

```text
03_classify/<project>.AFEPSI.detail.tsv
03_classify/<project>.ALEPSI.detail.tsv
03_classify/<project>.TE_overlap.AFEPSI.detail.tsv
03_classify/<project>.TE_overlap.ALEPSI.detail.tsv
```

AFE/ALE detail TE-overlap tables add:

| Column | Meaning |
| --- | --- |
| `te_overlap_bp_max` | Maximum base-pair overlap between event exon and TE annotations. |
| `te_overlap_frac_max` | Maximum fractional overlap between event exon and TE annotations. |
| `te_boundary_side` | TE boundary side in transcript/exon direction. |
| `junction_te_side_reads_max` | Maximum TE-side junction read count. Present only when junction evidence is effective. |

## Debug Outputs

Debug mode additionally retains selected reusable HITindex outputs:

```text
03_classify/HITindex/<replicate>.exon
03_classify/HITindex/<replicate>.AFEPSI        # only when AFE/ALE effective
03_classify/HITindex/<replicate>.ALEPSI        # only when AFE/ALE effective
03_classify/HITindex/<project>.te_overlap_transcript_exon_junction_support.tsv
logs/qual_config.json
logs/qual_extend.log
```

### te_overlap_transcript_exon_junction_support.tsv

| Column | Meaning |
| --- | --- |
| `sample` | Sample or replicate ID. |
| `gene_id` | Gene ID associated with the exon. |
| `transcript_id` | Transcript ID containing the exon. |
| `exon` | Exon coordinate used by the junction-support table. |
| `chrom` | Exon chromosome. |
| `start` | Exon start coordinate. |
| `end` | Exon end coordinate. |
| `strand` | Exon strand. |
| `nleft` | Junction read count on the left genomic side. |
| `nright` | Junction read count on the right genomic side. |
| `total_junction_reads` | Total junction reads across both sides. |
| `junction_supported` | Whether junction support passes the support rule. |
| `te_overlap_label` | TE-overlap label after junction evidence evaluation. |
| `te_overlap_bp_max` | Maximum TE-overlap base pairs. |
| `te_overlap_frac_max` | Maximum TE-overlap fraction. |
| `te_boundary_hit_any` | Whether any TE boundary hit was found. |
| `te_splice_site_repeat_TE` | TE annotation supporting splice-site repeat evidence. |

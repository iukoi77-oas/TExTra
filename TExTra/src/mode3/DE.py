"""Differential TE-derived exon usage tests and output formatting."""

import os
import glob
import pandas as pd

from util.common.write_logs import log_message
from util.common.define_layout import QUANT_DIR

ZERO_TOL = 1e-6
SUMMARY_ROUND_DECIMALS = 6
DEFAULT_NONZERO_USAGE_CUTOFF = 0.05
DEFAULT_SWITCH_LOW_USAGE = 0.05
DEFAULT_SWITCH_HIGH_USAGE = 0.95
DEFAULT_EMPIRICAL_AREA = 1000
DIFF_META_OUTPUT_COLUMNS = [
    "metaexon_id",
    "gene_id",
    "transcript_id",
    "te_overlap_label",
    "te_boundary_side",
    "te_splice_site_repeat_TE",
]


def _fmt_int(value):
    return f"{int(value):,}"


def _fmt_pct(count, total):
    total = int(total)
    if total <= 0:
        return "0.0%"
    return f"{float(count) / float(total) * 100.0:.1f}%"


def _fmt_count_pct(count, total):
    return f"{_fmt_int(count)} ({_fmt_pct(count, total)})"


def _format_counts(series, denominator=None, order=None):
    counts = series.fillna("NA").astype(str).value_counts().to_dict()
    keys = list(order or [])
    keys.extend(sorted(k for k in counts if k not in keys))
    total = int(denominator if denominator is not None else sum(counts.values()))
    items = []
    for key in keys:
        count = int(counts.get(key, 0))
        if count <= 0:
            continue
        items.append(f"{key}={_fmt_count_pct(count, total)}")
    return ", ".join(items) if items else "none"


def _format_bool(value):
    return "on" if bool(value) else "off"


def _higher_usage_group(delta_usage, group1, group2):
    delta = pd.to_numeric(pd.Series([delta_usage]), errors="coerce").iloc[0]
    if pd.isna(delta) or abs(float(delta)) <= ZERO_TOL:
        return "tie"
    return group1 if float(delta) > 0 else group2


def _add_missing_columns(df, columns):
    out = df.copy()
    for col in columns:
        if col not in out.columns:
            out[col] = ""
    return out


def _build_diff_output_table(
    result_df,
    meta_df,
    usage_df,
    selected_samples,
    group1,
    group2,
    include_samples=False,
    include_detail=False,
):
    """Normalize public diff result columns for default/detail/debug contracts."""
    rename_map = {
        "mean_A": "mean_usage_group1",
        "mean_B": "mean_usage_group2",
        "n_non_na_A": "n_non_na_group1",
        "n_non_na_B": "n_non_na_group2",
    }
    out = result_df.copy().rename(columns=rename_map)
    if include_samples:
        sample_df = usage_df[["exon_id"] + selected_samples].copy()
        sample_df["exon_id"] = sample_df["exon_id"].astype(str)
        out["exon_id"] = out["exon_id"].astype(str) if "exon_id" in out.columns else ""
        out = pd.merge(out, sample_df.drop_duplicates(subset=["exon_id"]), on="exon_id", how="left")

    out["group1"] = group1
    out["group2"] = group2
    if "delta_usage" in out.columns:
        out["higher_usage_group"] = out["delta_usage"].apply(lambda x: _higher_usage_group(x, group1, group2))
    else:
        out["higher_usage_group"] = ""

    meta_cols = [col for col in DIFF_META_OUTPUT_COLUMNS if col in meta_df.columns]
    if not meta_df.empty and meta_cols:
        meta_merge = meta_df[["exon_id"] + meta_cols].drop_duplicates(subset=["exon_id"]).copy()
        meta_merge["exon_id"] = meta_merge["exon_id"].astype(str)
        out["exon_id"] = out["exon_id"].astype(str) if "exon_id" in out.columns else ""
        out = pd.merge(out, meta_merge, on="exon_id", how="left")

    default_cols = [
        "exon_id",
        "group1",
        "group2",
        "mean_usage_group1",
        "mean_usage_group2",
        "delta_usage",
        "higher_usage_group",
        "pvalue",
        "padj",
    ] + meta_cols
    detail_cols = (
        ["exon_id"]
        + list(selected_samples)
        + [
            "group1",
            "group2",
            "mean_usage_group1",
            "mean_usage_group2",
            "delta_usage",
            "higher_usage_group",
            "pvalue",
            "padj",
            "test_method",
            "n_non_na_group1",
            "n_non_na_group2",
            "usage_pattern_flag",
        ]
        + meta_cols
    )
    columns = detail_cols if include_detail else default_cols
    out = _add_missing_columns(out, columns)
    return out[columns]


DIFF_METADATA_COLUMNS = {
    "metaexon",
    "metaexon_id",
    "exon",
    "exon_id",
    "event_id",
    "gene",
    "gene_id",
    "transcript_id",
    "chrom",
    "start",
    "end",
    "strand",
    "source",
    "te_overlap_label",
    "te_overlap_label_raw",
    "te_overlap_label_final",
    "te_overlap_n",
    "te_overlap_bp_max",
    "te_overlap_frac_max",
    "te_boundary_side",
    "te_boundary_hit_any",
    "te_overlap_pass_any",
    "te_splice_site_repeat_TE",
    "te_other_overlap_TE",
    "usage_definition",
    "expressed_any",
    "ID_position",
    "ID_position_summary",
    "ID_position_detail",
    "ID_position_source",
    "ID_position_roles",
    "ID_position_confidence",
    "transcript_structure_roles",
    "HITindex_structure_roles",
    "position_evidence_relation",
    "candidate_TE_event",
    "candidate_TE_confidence",
    "ID_position_evaluable_sample_n",
    "ID_position_evaluable_replicate_n",
    "HITindex_evaluable_sample_n",
    "HITindex_evaluable_replicate_n",
    "support_sample_n",
    "support_sample_ratio",
    "junction_evidence_enabled",
    "junction_support_sample_n",
    "junction_support_sample_ids",
    "junction_te_side_reads_max",
    "junction_te_side_reads_mean_supported_samples",
    "junction_supported",
    "junction_pass",
    "te_overlap_degraded",
    "te_overlap_degrade_reason",
}


def _apply_zero_tolerance(df):
    """Set near-zero numeric values to exact zero for stable summaries/labels."""
    out = df.copy()
    num_cols = out.select_dtypes(include=["floating"]).columns
    if len(num_cols) > 0:
        out[num_cols] = out[num_cols].mask(out[num_cols].abs() < ZERO_TOL, 0.0)
    return out


def _resolve_te_exon_table(quant_dir, projectname):
    quant_folder = os.path.join(quant_dir, QUANT_DIR)
    expected = os.path.join(quant_folder, f"{projectname}.TE_overlap.exon_usage.tsv")
    if os.path.isfile(expected):
        return expected

    all_tables = sorted(glob.glob(os.path.join(quant_folder, "*.TE_overlap.exon_usage.tsv")))
    if len(all_tables) == 1:
        return all_tables[0]
    if len(all_tables) > 1:
        raise RuntimeError(
            f"Multiple TE-exon usage tables found in {quant_folder}; "
            "please pass --project matching quant output."
        )
    return None


def _compose_te_info(row):
    tokens = []
    for col in ["te_splice_site_repeat_TE", "te_other_overlap_TE", "TE_info"]:
        if col not in row:
            continue
        for te in str(row[col]).split(","):
            te = te.strip()
            if te and te.lower() != "nan":
                tokens.append(te)
    return ",".join(sorted(set(tokens)))


def _infer_condition(sample_name):
    parts = str(sample_name).split("_")
    if len(parts) < 2:
        raise ValueError(
            f"Sample name '{sample_name}' does not follow <condition>_<replicate> format"
        )
    return "_".join(parts[:-1])


def _parse_group_pair(samples_arg):
    groups = [x.strip() for x in str(samples_arg).split(",") if x.strip()]
    if len(groups) != 2:
        raise RuntimeError(
            f"--groups must contain exactly two condition names separated by comma; got: {samples_arg}"
        )
    if groups[0] == groups[1]:
        raise RuntimeError(f"--groups conditions must be different; got: {samples_arg}")
    return groups[0], groups[1]


def _bh_adjust(pvalues):
    """Benjamini-Hochberg adjusted p-values without requiring statsmodels."""
    p = pd.to_numeric(pd.Series(pvalues), errors="coerce")
    out = pd.Series([float("nan")] * len(p), index=p.index, dtype="float64")
    valid = p.notna()
    if valid.sum() == 0:
        return out.tolist()

    valid_p = p.loc[valid].astype(float)
    order = valid_p.sort_values(kind="mergesort").index.tolist()
    m = float(len(order))
    prev = 1.0
    adjusted = {}
    for rank_from_end, idx in enumerate(reversed(order), start=1):
        rank = m - rank_from_end + 1.0
        val = min(prev, float(valid_p.loc[idx]) * m / rank)
        prev = val
        adjusted[idx] = min(max(val, 0.0), 1.0)
    for idx, val in adjusted.items():
        out.loc[idx] = val
    return out.tolist()


def _resolve_sample_columns(df, requested_groups):
    """Select only replicate columns whose inferred condition is requested."""
    sample_to_condition = {}
    available_conditions = set()
    skipped_malformed = []
    for col in df.columns:
        if col in DIFF_METADATA_COLUMNS:
            continue
        vals = pd.to_numeric(df[col], errors="coerce")
        if not vals.notna().any():
            continue
        try:
            condition = _infer_condition(col)
        except ValueError:
            skipped_malformed.append(str(col))
            continue
        available_conditions.add(condition)
        if condition in requested_groups:
            sample_to_condition[str(col)] = condition

    if not sample_to_condition:
        raise RuntimeError(
            "No sample columns matched requested groups. "
            f"Requested groups: {', '.join(requested_groups)}. "
            f"Available inferred conditions: {', '.join(sorted(available_conditions)) or 'none'}."
        )
    selected_coldata = pd.DataFrame(
        {
            "sample": list(sample_to_condition.keys()),
            "condition": list(sample_to_condition.values()),
        }
    )
    counts_per_group = selected_coldata["condition"].value_counts()
    missing_groups = [g for g in requested_groups if g not in counts_per_group]
    if missing_groups:
        raise RuntimeError(
            f"No samples found for condition(s): {', '.join(missing_groups)}. "
            f"Available inferred conditions: {', '.join(sorted(available_conditions)) or 'none'}."
        )
    if skipped_malformed:
        log_message(
            "[WARNING]",
            "Ignored numeric column(s) that do not follow <condition>_<replicate>: "
            + ", ".join(skipped_malformed[:10])
            + (" ..." if len(skipped_malformed) > 10 else ""),
            color="warning",
        )
    return selected_coldata.reset_index(drop=True), sorted(available_conditions)


def _classical_pvalue(g1_values, g2_values, paired=False):
    """Classical replicate-level PSI/usage test, following SUPPA's classical option."""
    try:
        if paired:
            from scipy.stats import wilcoxon

            paired_vals = [
                (float(a), float(b))
                for a, b in zip(g1_values, g2_values)
                if pd.notna(a) and pd.notna(b)
            ]
            if len(paired_vals) < 2:
                return float("nan")
            a_vals, b_vals = zip(*paired_vals)
            if all(abs(float(a) - float(b)) <= ZERO_TOL for a, b in paired_vals):
                return 1.0
            return float(wilcoxon(a_vals, b_vals).pvalue)

        from scipy.stats import mannwhitneyu

        a_vals = [float(x) for x in g1_values if pd.notna(x)]
        b_vals = [float(x) for x in g2_values if pd.notna(x)]
        if len(a_vals) < 1 or len(b_vals) < 1:
            return float("nan")
        if len(a_vals) + len(b_vals) < 3:
            return float("nan")
        if all(abs(x - a_vals[0]) <= ZERO_TOL for x in a_vals + b_vals):
            return 1.0
        return float(mannwhitneyu(a_vals, b_vals, alternative="two-sided").pvalue)
    except ImportError as exc:
        raise RuntimeError(
            "Python package 'scipy' is required for --test-method classical."
        ) from exc
    except ValueError:
        return 1.0


def _run_classical(usage_df, selected_coldata, group1, group2, paired=False):
    g1_samples = selected_coldata.loc[selected_coldata["condition"] == group1, "sample"].astype(str).tolist()
    g2_samples = selected_coldata.loc[selected_coldata["condition"] == group2, "sample"].astype(str).tolist()
    rows = []
    for _, row in usage_df.iterrows():
        g1_vals = pd.to_numeric(row[g1_samples], errors="coerce")
        g2_vals = pd.to_numeric(row[g2_samples], errors="coerce")
        mean_a = float(g1_vals.mean(skipna=True)) if g1_vals.notna().any() else 0.0
        mean_b = float(g2_vals.mean(skipna=True)) if g2_vals.notna().any() else 0.0
        pvalue = _classical_pvalue(g1_vals.tolist(), g2_vals.tolist(), paired=paired)
        rows.append(
            {
                "exon_id": str(row["exon_id"]),
                "mean_A": mean_a,
                "mean_B": mean_b,
                "delta_usage": mean_a - mean_b,
                "pvalue": pvalue,
                "n_non_na_A": int(g1_vals.notna().sum()),
                "n_non_na_B": int(g2_vals.notna().sum()),
                "test_method": "classical_paired" if paired else "classical_unpaired",
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        out["padj"] = []
    else:
        out["padj"] = _bh_adjust(out["pvalue"])
        for c in ["mean_A", "mean_B", "delta_usage"]:
            out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0.0).round(SUMMARY_ROUND_DECIMALS)
    return _apply_zero_tolerance(out)


def _load_event_transcript_map(quant_dir, projectname, event_ids):
    usage_path = _resolve_te_exon_table(quant_dir, projectname)
    if usage_path is None:
        raise FileNotFoundError(
            f"TE-overlap exon usage table not found under {os.path.join(quant_dir, QUANT_DIR)}."
        )
    df = pd.read_csv(usage_path, sep="\t", dtype=str).fillna("")
    if "exon_id" not in df.columns or "transcript_id" not in df.columns:
        raise RuntimeError(
            f"Invalid empirical event-transcript table: {usage_path} must include exon_id and transcript_id."
        )
    wanted = set(str(x) for x in event_ids)
    event_to_txs = {}
    for _, row in df.iterrows():
        event_id = str(row["exon_id"])
        if event_id not in wanted or event_id in event_to_txs:
            continue
        txs = [x.strip() for x in str(row["transcript_id"]).split(",") if x.strip()]
        event_to_txs[event_id] = txs
    missing = sorted(wanted - set(event_to_txs))
    if missing:
        raise RuntimeError(
            f"empirical mode could not find supporting transcripts for {len(missing)} event(s) in {usage_path}."
        )
    return event_to_txs, usage_path


def _load_transcript_tpm_by_sample(quant_dir, sample_names):
    quant_folder = os.path.join(quant_dir, QUANT_DIR)
    abundance_path = os.path.join(quant_folder, "transcript_abundance.tsv")
    if not os.path.isfile(abundance_path):
        raise FileNotFoundError(
            f"empirical mode requires transcript abundance table from quant output: {abundance_path}"
        )
    df = pd.read_csv(abundance_path, sep="\t")
    required = {"sample", "transcript_id", "tpm"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise RuntimeError(f"Invalid transcript abundance table {abundance_path}; missing column(s): {', '.join(missing)}")
    out = {}
    for sample in sample_names:
        sub = df[df["sample"].astype(str) == str(sample)].copy()
        if sub.empty:
            raise FileNotFoundError(
                f"empirical mode could not find sample '{sample}' in transcript abundance table: {abundance_path}"
            )
        out[sample] = pd.to_numeric(sub["tpm"], errors="coerce").fillna(0.0).groupby(sub["transcript_id"].astype(str)).sum().to_dict()
    return out


def _event_log_tpm(event_txs, sample_tpm):
    import math

    total = float(sum(float(sample_tpm.get(tx, 0.0)) for tx in event_txs))
    if total <= 0:
        return float("nan")
    return math.log10(total)


def _empirical_pvalue(local_abs_deltas, observed_abs_delta):
    vals = [float(x) for x in local_abs_deltas if pd.notna(x)]
    if not vals:
        return float("nan")
    ge_n = sum(1 for x in vals if x >= float(observed_abs_delta))
    return min(max((ge_n + 1.0) / (len(vals) + 1.0), 0.0), 1.0)


def _run_empirical(usage_df, selected_coldata, group1, group2, quant_dir, projectname, area=DEFAULT_EMPIRICAL_AREA):
    sample_names = selected_coldata["sample"].astype(str).tolist()
    event_ids = usage_df["exon_id"].astype(str).tolist()
    event_to_txs, usage_tx_path = _load_event_transcript_map(quant_dir, projectname, event_ids)
    sample_tpm = _load_transcript_tpm_by_sample(quant_dir, sample_names)

    tpm_log = {}
    for event_id, txs in event_to_txs.items():
        tpm_log[event_id] = {sample: _event_log_tpm(txs, sample_tpm[sample]) for sample in sample_names}

    g1_samples = selected_coldata.loc[selected_coldata["condition"] == group1, "sample"].astype(str).tolist()
    g2_samples = selected_coldata.loc[selected_coldata["condition"] == group2, "sample"].astype(str).tolist()

    background = []
    for _, row in usage_df.iterrows():
        event_id = str(row["exon_id"])
        for group_samples in (g1_samples, g2_samples):
            for i in range(len(group_samples)):
                for j in range(i + 1, len(group_samples)):
                    s1, s2 = group_samples[i], group_samples[j]
                    u1 = pd.to_numeric(pd.Series([row[s1]]), errors="coerce").iloc[0]
                    u2 = pd.to_numeric(pd.Series([row[s2]]), errors="coerce").iloc[0]
                    l1 = tpm_log.get(event_id, {}).get(s1, float("nan"))
                    l2 = tpm_log.get(event_id, {}).get(s2, float("nan"))
                    if pd.notna(u1) and pd.notna(u2) and pd.notna(l1) and pd.notna(l2):
                        background.append((abs(float(u2) - float(u1)), (float(l1) + float(l2)) / 2.0))
    if not background:
        raise RuntimeError(
            "empirical mode could not build a replicate background distribution. "
            "At least two replicates per condition with positive supporting transcript TPM are required."
        )
    background = sorted(background, key=lambda x: x[1])

    rows = []
    window_n = max(1, int(area))
    for _, row in usage_df.iterrows():
        event_id = str(row["exon_id"])
        g1_vals = pd.to_numeric(row[g1_samples], errors="coerce")
        g2_vals = pd.to_numeric(row[g2_samples], errors="coerce")
        mean_a = float(g1_vals.mean(skipna=True)) if g1_vals.notna().any() else 0.0
        mean_b = float(g2_vals.mean(skipna=True)) if g2_vals.notna().any() else 0.0
        delta = mean_a - mean_b
        event_logs = [v for v in tpm_log.get(event_id, {}).values() if pd.notna(v)]
        if event_logs:
            event_log_tpm = sum(event_logs) / float(len(event_logs))
            ranked = sorted(background, key=lambda x: abs(x[1] - event_log_tpm))
            local = [x[0] for x in ranked[: min(window_n, len(ranked))]]
            pvalue = _empirical_pvalue(local, abs(delta))
            local_n = len(local)
        else:
            event_log_tpm = float("nan")
            pvalue = float("nan")
            local_n = 0
        rows.append(
            {
                "exon_id": event_id,
                "mean_A": mean_a,
                "mean_B": mean_b,
                "delta_usage": delta,
                "pvalue": pvalue,
                "padj": float("nan"),
                "n_non_na_A": int(g1_vals.notna().sum()),
                "n_non_na_B": int(g2_vals.notna().sum()),
                "event_avg_log10_tpm": event_log_tpm,
                "empirical_local_background_n": int(local_n),
                "empirical_background_n": int(window_n),
                "event_transcript_map": usage_tx_path,
                "test_method": "empirical",
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["padj"] = _bh_adjust(out["pvalue"])
        for c in ["mean_A", "mean_B", "delta_usage", "event_avg_log10_tpm"]:
            out[c] = pd.to_numeric(out[c], errors="coerce").round(SUMMARY_ROUND_DECIMALS)
    return _apply_zero_tolerance(out)


def _classify_usage_events(usage_df, g1_samples, g2_samples):
    """
    Add usage-pattern annotations while only removing rows with no usable signal.

    near_binary_switch and sparse_replicate_support are descriptive flags only;
    they do not exclude events from classical/empirical statistical testing.
    """
    sample_cols = [c for c in usage_df.columns if c != "exon_id"]
    mat = usage_df[sample_cols].apply(pd.to_numeric, errors="coerce").astype("float64")
    mat = mat.mask(mat.abs() < ZERO_TOL, 0.0)
    g1 = mat[g1_samples]
    g2 = mat[g2_samples]

    summary = pd.DataFrame({"exon_id": usage_df["exon_id"].astype(str)})
    summary["mean_A"] = g1.mean(axis=1, skipna=True).fillna(0.0)
    summary["mean_B"] = g2.mean(axis=1, skipna=True).fillna(0.0)
    summary["delta_usage"] = summary["mean_A"] - summary["mean_B"]
    summary["n_non_na_A"] = g1.notna().sum(axis=1).astype(int)
    summary["n_non_na_B"] = g2.notna().sum(axis=1).astype(int)
    summary["n_non_na_total"] = mat.notna().sum(axis=1).astype(int)
    summary["n_nonzero_A"] = (g1.fillna(0.0) > DEFAULT_NONZERO_USAGE_CUTOFF).sum(axis=1).astype(int)
    summary["n_nonzero_B"] = (g2.fillna(0.0) > DEFAULT_NONZERO_USAGE_CUTOFF).sum(axis=1).astype(int)
    summary["n_nonzero_total"] = summary["n_nonzero_A"] + summary["n_nonzero_B"]
    row_max = mat.max(axis=1, skipna=True)
    row_min = mat.min(axis=1, skipna=True)
    summary["row_range"] = (row_max - row_min).fillna(0.0)

    all_na = summary["n_non_na_total"].eq(0)
    insufficient = (
        summary["n_non_na_total"].lt(2)
        | summary["n_non_na_A"].lt(1)
        | summary["n_non_na_B"].lt(1)
    )
    constant = (~all_na) & summary["row_range"].abs().le(1e-12)
    switch_like = (
        (~all_na)
        & (~insufficient)
        & (~constant)
        & (
            ((summary["mean_A"] <= DEFAULT_SWITCH_LOW_USAGE) & (summary["mean_B"] >= DEFAULT_SWITCH_HIGH_USAGE))
            | ((summary["mean_B"] <= DEFAULT_SWITCH_LOW_USAGE) & (summary["mean_A"] >= DEFAULT_SWITCH_HIGH_USAGE))
        )
    )
    sparse = (
        (~all_na)
        & (~insufficient)
        & (~constant)
        & (~switch_like)
        & (
            summary["n_nonzero_total"].le(1)
            | ((summary["n_nonzero_A"] <= 1) & (summary["mean_B"] <= DEFAULT_SWITCH_LOW_USAGE))
            | ((summary["n_nonzero_B"] <= 1) & (summary["mean_A"] <= DEFAULT_SWITCH_LOW_USAGE))
        )
    )

    summary["filter_status"] = "tested"
    summary.loc[all_na, "filter_status"] = "removed_all_na"
    summary.loc[insufficient, "filter_status"] = "removed_insufficient_samples"
    summary.loc[constant, "filter_status"] = "removed_constant_all_samples"
    summary["usage_pattern_flag"] = "none"
    summary.loc[switch_like, "usage_pattern_flag"] = "near_binary_switch"
    summary.loc[sparse, "usage_pattern_flag"] = "sparse_replicate_support"

    summary["support_note"] = ""
    summary.loc[summary["usage_pattern_flag"] == "near_binary_switch", "support_note"] = "near-binary group switch"
    summary.loc[summary["usage_pattern_flag"] == "sparse_replicate_support", "support_note"] = (
        "single-sample-driven or sparse non-zero usage"
    )
    summary.loc[summary["filter_status"] == "removed_constant_all_samples", "support_note"] = (
        "all selected samples have identical usage"
    )
    summary.loc[summary["filter_status"] == "removed_insufficient_samples", "support_note"] = (
        "insufficient non-NA samples for group comparison"
    )
    summary.loc[summary["filter_status"] == "removed_all_na", "support_note"] = "all values are NA"

    summary_num_cols = [
        "mean_A",
        "mean_B",
        "delta_usage",
        "row_range",
    ]
    summary[summary_num_cols] = summary[summary_num_cols].round(SUMMARY_ROUND_DECIMALS)

    tested_ids = set(summary.loc[summary["filter_status"] == "tested", "exon_id"].astype(str))
    flag_ids = set(
        summary.loc[
            (summary["filter_status"] == "tested") & (summary["usage_pattern_flag"] != "none"),
            "exon_id",
        ].astype(str)
    )

    tested_df = usage_df[usage_df["exon_id"].astype(str).isin(tested_ids)].copy()
    pattern_df = summary[summary["exon_id"].astype(str).isin(flag_ids)].copy()

    return tested_df, pattern_df, summary


def DE_func(args):
    """Run differential exon-usage analysis on exon usage matrix."""
    outfolder = args.out_dir
    method = str(getattr(args, "test_method", "classical")).strip().lower()
    projectname = args.project
    quant_dir = args.quant
    usage_delta_cutoff = float(getattr(args, "delta_exon_usage", 0.1))
    padj_cutoff = float(args.padj)
    raw_pvalue_cutoff = getattr(args, "pvalue", None)
    raw_pvalue_cutoff = None if raw_pvalue_cutoff is None else float(raw_pvalue_cutoff)
    paired = bool(getattr(args, "paired", False))
    empirical_background_n = int(getattr(args, "empirical_background_n", DEFAULT_EMPIRICAL_AREA))

    os.makedirs(outfolder, exist_ok=True)
    de_out = os.path.join(outfolder, "DE")
    os.makedirs(de_out, exist_ok=True)

    if method not in {"classical", "empirical"}:
        raise RuntimeError(f"Unsupported diff method: {method}")
    verbosity = bool(getattr(args, "detail", False) or getattr(args, "debug", False))
    total_steps = int(getattr(args, "diff_total_steps", 1))

    log_message(
        "[INFO]",
        f"Step 1/{total_steps}: Differential usage analysis",
        bold=True,
        color="step",
    )
    pvalue_text = "off" if raw_pvalue_cutoff is None else f"<={raw_pvalue_cutoff:g}"
    log_message(
        "[INFO]",
        "Diff parameters: "
        f"groups={args.groups}, method={method}, padj<{padj_cutoff:g}, "
        f"|delta_exon_usage|>{usage_delta_cutoff:g}, pvalue={pvalue_text}, paired={_format_bool(paired)}",
    )
    if getattr(args, "ncpred", False):
        log_message(
            "[INFO]",
            "ncPred parameters: "
            f"model={getattr(args, 'ncpred_model', 've')}, "
            f"min_length={int(getattr(args, 'min_length', 200))}, "
            f"PLEK2={os.path.abspath(getattr(args, 'plek2_path', ''))}",
        )

    quant_table = _resolve_te_exon_table(quant_dir, projectname)
    if quant_table is None:
        raise FileNotFoundError(
            f"TE-exon usage table not found under {os.path.join(quant_dir, QUANT_DIR)}. "
            f"Expected the current quant contract output {projectname}.TE_overlap.exon_usage.tsv."
        )
    group1_name, group2_name = _parse_group_pair(args.groups)
    df = pd.read_csv(quant_table, sep="\t")

    id_col = "exon_id"
    if id_col not in df.columns:
        raise RuntimeError(f"{quant_table} must include column: exon_id.")
    te_label_col = "te_overlap_label"
    if te_label_col not in df.columns:
        raise RuntimeError(f"{quant_table} must include column: te_overlap_label.")
    df = df[df[te_label_col].fillna("no_overlap") == "TE_overlap"].copy()
    if df.empty:
        raise RuntimeError(
            f"No TE-overlap exons found in {quant_table}; cannot run diff."
        )

    selected_coldata, _available_conditions = _resolve_sample_columns(df, {group1_name, group2_name})
    selected_samples = selected_coldata["sample"].astype(str).tolist()
    group_sample_counts = selected_coldata["condition"].value_counts().to_dict()
    log_message(
        "[INFO]",
        f"Input TE-overlap usage: events={_fmt_int(len(df))}, samples={_fmt_int(len(selected_samples))}, source={os.path.abspath(quant_table)}",
    )
    log_message(
        "[INFO]",
        "Group comparison: "
        f"{group1_name}={_fmt_int(group_sample_counts.get(group1_name, 0))}, "
        f"{group2_name}={_fmt_int(group_sample_counts.get(group2_name, 0))}",
    )
    usage_df = df[[id_col] + selected_samples].rename(columns={id_col: "exon_id"}).copy()
    meta_cols_present = [col for col in DIFF_META_OUTPUT_COLUMNS if col in df.columns]
    meta_df = df[[id_col] + meta_cols_present].rename(columns={id_col: "exon_id"}).copy()
    meta_df["exon_id"] = meta_df["exon_id"].astype(str)
    meta_df = meta_df.drop_duplicates(subset=["exon_id"])

    for sample in selected_samples:
        usage_df[sample] = pd.to_numeric(usage_df[sample], errors="coerce")
    if usage_df[selected_samples].apply(pd.to_numeric, errors="coerce").fillna(0.0).sum().sum() == 0:
        log_message(
            "[WARNING]",
            "All selected final TE-overlap exons have zero usage in the requested samples; "
            "diff will write preprocessing summaries and empty statistical result tables.",
            color="warning",
        )

    g1_samples = selected_coldata.loc[selected_coldata["condition"] == group1_name, "sample"].astype(str).tolist()
    g2_samples = selected_coldata.loc[selected_coldata["condition"] == group2_name, "sample"].astype(str).tolist()

    tested_df, pattern_df, preprocess_df = _classify_usage_events(
        usage_df=usage_df[["exon_id"] + selected_samples].copy(),
        g1_samples=g1_samples,
        g2_samples=g2_samples,
    )
    if method == "empirical" and not tested_df.empty:
        _load_event_transcript_map(quant_dir, projectname, tested_df["exon_id"].astype(str).tolist())
        _load_transcript_tpm_by_sample(quant_dir, selected_samples)
    preprocess_df = _apply_zero_tolerance(preprocess_df)
    total_events = len(preprocess_df)
    tested_events = int((preprocess_df["filter_status"] == "tested").sum()) if "filter_status" in preprocess_df else len(tested_df)
    filtered_events = max(0, int(total_events) - int(tested_events))
    log_message(
        "[INFO]",
        "Usage preprocessing: "
        f"total={_fmt_int(total_events)}, tested={_fmt_count_pct(tested_events, total_events)}, "
        f"filtered={_fmt_count_pct(filtered_events, total_events)}",
    )

    full_out = os.path.join(de_out, "differential_usage.tsv")
    preprocess_out = os.path.join(de_out, "usage_preprocess_summary.tsv")
    canonical_out = os.path.join(de_out, "differential_usage.tsv")
    canonical_sig_out = os.path.join(de_out, "differential_significant_usage.tsv")
    internal_cols = [
        "exon_id",
        "mean_A",
        "mean_B",
        "delta_usage",
        "pvalue",
        "padj",
        "n_non_na_A",
        "n_non_na_B",
        "test_method",
        "usage_pattern_flag",
        "filter_status",
        "support_note",
    ]

    if tested_df.empty:
        result_df = pd.DataFrame(
            columns=internal_cols
            + [
                "usage_delta_cutoff",
                "padj_cutoff",
                "pvalue_cutoff",
                "significant_usage",
            ]
        )
    elif method == "classical":
        result_df = _run_classical(
            usage_df=tested_df[["exon_id"] + selected_samples].copy(),
            selected_coldata=selected_coldata,
            group1=group1_name,
            group2=group2_name,
            paired=paired,
        )
    else:
        result_df = _run_empirical(
            usage_df=tested_df[["exon_id"] + selected_samples].copy(),
            selected_coldata=selected_coldata,
            group1=group1_name,
            group2=group2_name,
            quant_dir=quant_dir,
            projectname=projectname,
            area=empirical_background_n,
        )

    if not result_df.empty:
        preprocess_anno = preprocess_df[
            ["exon_id", "usage_pattern_flag", "filter_status", "support_note"]
        ].drop_duplicates(subset=["exon_id"]).copy()
        result_df["exon_id"] = result_df["exon_id"].astype(str)
        preprocess_anno["exon_id"] = preprocess_anno["exon_id"].astype(str)
        result_df = pd.merge(result_df, preprocess_anno, on="exon_id", how="left")
        result_df["usage_pattern_flag"] = result_df["usage_pattern_flag"].fillna("none")
        result_df["filter_status"] = result_df["filter_status"].fillna("tested")
        result_df["support_note"] = result_df["support_note"].fillna("")

    if not result_df.empty:
        result_df["significant_usage"] = (
            pd.to_numeric(result_df["delta_usage"], errors="coerce").abs().gt(usage_delta_cutoff)
            & pd.to_numeric(result_df["padj"], errors="coerce").lt(padj_cutoff)
        )
        if raw_pvalue_cutoff is not None:
            result_df["significant_usage"] = result_df["significant_usage"] & pd.to_numeric(
                result_df["pvalue"], errors="coerce"
            ).lt(raw_pvalue_cutoff)
        result_df["significant_usage"] = result_df["significant_usage"].astype(int)
        result_df["usage_delta_cutoff"] = usage_delta_cutoff
        result_df["padj_cutoff"] = padj_cutoff
        result_df["pvalue_cutoff"] = raw_pvalue_cutoff if raw_pvalue_cutoff is not None else ""
        result_df = result_df[["exon_id"] + [c for c in result_df.columns if c != "exon_id"]]

    if getattr(args, "detail", False) or getattr(args, "debug", False):
        detail_out_df = _build_diff_output_table(
            result_df=result_df,
            meta_df=meta_df,
            usage_df=usage_df,
            selected_samples=selected_samples,
            group1=group1_name,
            group2=group2_name,
            include_samples=True,
            include_detail=True,
        )
        detail_out_df.to_csv(full_out, sep="\t", index=False)
        if "filter_status" in preprocess_df:
            log_message(
                "[INFO]",
                "Detail filter_status: "
                + _format_counts(
                    preprocess_df["filter_status"],
                    denominator=len(preprocess_df),
                    order=[
                        "tested",
                        "removed_constant_all_samples",
                        "removed_insufficient_samples",
                        "removed_all_na",
                    ],
                ),
            )
        if "usage_pattern_flag" in preprocess_df:
            log_message(
                "[INFO]",
                "Detail usage_pattern_flag: "
                + _format_counts(
                    preprocess_df["usage_pattern_flag"],
                    denominator=len(preprocess_df),
                    order=["none", "near_binary_switch", "sparse_replicate_support"],
                ),
            )
    else:
        for stale in [canonical_out]:
            if os.path.isfile(stale):
                os.remove(stale)
    if os.path.isfile(preprocess_out):
        os.remove(preprocess_out)

    if result_df.empty:
        sig_internal_df = pd.DataFrame(columns=result_df.columns)
    else:
        sig_internal_df = result_df.loc[result_df["significant_usage"] == 1].copy()
    sig_df = _build_diff_output_table(
        result_df=sig_internal_df,
        meta_df=meta_df,
        usage_df=usage_df,
        selected_samples=selected_samples,
        group1=group1_name,
        group2=group2_name,
        include_samples=False,
        include_detail=False,
    )
    sig_df.to_csv(canonical_sig_out, sep="\t", index=False)
    significant_n = len(sig_df)
    if not result_df.empty:
        log_message(
            "[INFO]",
            f"Differential usage summary: tested={_fmt_int(len(result_df))}, significant={_fmt_count_pct(significant_n, len(result_df))}",
        )
    else:
        log_message("[INFO]", "Differential usage summary: tested=0, significant=0 (0.0%)")
    if verbosity and not sig_df.empty and "delta_usage" in sig_df.columns:
        group1_higher = int(pd.to_numeric(sig_df["delta_usage"], errors="coerce").gt(0).sum())
        group2_higher = int(pd.to_numeric(sig_df["delta_usage"], errors="coerce").lt(0).sum())
        log_message(
            "[INFO]",
            "Detail significant direction: "
            f"{group1_name}_higher={_fmt_count_pct(group1_higher, significant_n)}, "
            f"{group2_name}_higher={_fmt_count_pct(group2_higher, significant_n)}",
        )
    if verbosity:
        log_message("[SUCCESS]", f"Differential usage detail: {os.path.abspath(full_out)}", color="success")
    stale_outputs = [
        f"{method}_all.tsv",
        f"{method}_significant_usage.tsv",
        f"{method}_input.tsv",
        "classical_all.tsv",
        "classical_significant_usage.tsv",
        "classical_input.tsv",
        "empirical_all.tsv",
        "empirical_significant_usage.tsv",
        "empirical_input.tsv",
        "coldata.txt",
        "exon_meta.from_exon_usage.tsv",
        "usage_pattern_flags.tsv",
        "usage_preprocess_summary.tsv",
    ]
    for name in sorted(set(stale_outputs)):
        path = os.path.join(de_out, name)
        if os.path.isfile(path):
            os.remove(path)

    log_message("[SUCCESS]", f"Significant differential usage: {os.path.abspath(canonical_sig_out)}", color="success")

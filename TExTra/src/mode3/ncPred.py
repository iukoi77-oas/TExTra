"""ncPred transcript selection and PLEK2 coding-potential prediction workflow."""

import os
import glob
import subprocess as sp
import shutil
import pandas as pd
from collections import defaultdict

from util.common.write_logs import log_message
from util.common.define_layout import QUANT_DIR, resolve_consensus_gtf
from util.common.external_tools import resolve_external_file


def _parse_attr(attr_text):
    attr = {}
    for part in attr_text.strip().split(";"):
        token = part.strip()
        if not token:
            continue
        if " " not in token:
            continue
        key, value = token.split(" ", 1)
        attr[key] = value.strip().strip('"')
    return attr


def _parse_coord(coord):
    left, strand = str(coord).rsplit(":", 1)
    chrom, span = left.split(":", 1)
    start_s, end_s = span.split("-", 1)
    return chrom, int(start_s), int(end_s), strand


def _normalize_sig_schema(df):
    if df.empty:
        return pd.DataFrame(columns=["exon_id", "branch_type", "coord", "gene", "TE_info"])
    if "exon_id" not in df.columns:
        if "event_id" in df.columns:
            df = df.rename(columns={"event_id": "exon_id"})
        else:
            first_col = df.columns[0]
            df = df.rename(columns={first_col: "exon_id"})
    if "branch_type" not in df.columns:
        df["branch_type"] = df["exon_id"].astype(str).str.split("|", n=1, regex=False).str[0]
    if "coord" not in df.columns:
        df["coord"] = df["exon_id"].astype(str).str.split("|", n=1, regex=False).str[-1]
    if "gene" not in df.columns:
        df["gene"] = df["gene_id"].astype(str) if "gene_id" in df.columns else ""
    if "TE_info" not in df.columns:
        df["TE_info"] = df["te_splice_site_repeat_TE"].astype(str) if "te_splice_site_repeat_TE" in df.columns else ""
    return df[["exon_id", "branch_type", "coord", "gene", "TE_info"]].drop_duplicates()


def _load_sig_events(table_path, padj, preselected=False):
    df = pd.read_csv(table_path, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=["exon_id", "branch_type", "coord", "gene", "TE_info"])
    if preselected:
        sig = df.copy()
    else:
        if "significant_usage" in df.columns:
            sig = df[pd.to_numeric(df["significant_usage"], errors="coerce").fillna(0).astype(int) == 1].copy()
        elif "padj" in df.columns:
            df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
            sig = df[df["padj"] < padj].copy()
        else:
            return pd.DataFrame(columns=["exon_id", "branch_type", "coord", "gene", "TE_info"])
    return _normalize_sig_schema(sig)


def _load_transcript_exons(gtf_path):
    transcript_data = {}
    exon_index = defaultdict(list)
    with open(gtf_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom = fields[0]
            start0 = int(fields[3]) - 1
            end0 = int(fields[4])
            strand = fields[6]
            attrs = _parse_attr(fields[8])
            tid = attrs.get("transcript_id")
            gid = attrs.get("gene_id", "")
            if not tid:
                continue
            if tid not in transcript_data:
                transcript_data[tid] = {"gene": gid, "chrom": chrom, "strand": strand, "exons": []}
            transcript_data[tid]["exons"].append((start0, end0))

    for tid, info in transcript_data.items():
        exons = sorted(info["exons"], key=lambda x: x[0])
        n = len(exons)
        if n == 0:
            continue
        for i, (start0, end0) in enumerate(exons):
            exon_type = "internal"
            if n == 1:
                exon_type = "single"
            elif info["strand"] == "+":
                if i == 0:
                    exon_type = "first"
                elif i == n - 1:
                    exon_type = "last"
            else:
                if i == 0:
                    exon_type = "last"
                elif i == n - 1:
                    exon_type = "first"
            exon_index[(info["chrom"], info["strand"])].append(
                {
                    "transcript_id": tid,
                    "gene": info["gene"],
                    "start": start0,
                    "end": end0,
                    "exon_type": exon_type,
                }
            )
    return transcript_data, exon_index


def _select_sig_transcripts(sig_events, exon_index):
    selected = defaultdict(lambda: {"gene": "", "events": []})
    for row in sig_events.itertuples(index=False):
        try:
            chrom, start, end, strand = _parse_coord(row.coord)
        except ValueError:
            continue
        for exon in exon_index.get((chrom, strand), []):
            overlap_len = min(end, exon["end"]) - max(start, exon["start"])
            if overlap_len <= 0:
                continue
            tid = exon["transcript_id"]
            selected[tid]["gene"] = exon["gene"]
            selected[tid]["events"].append(
                {
                    "event_id": row.exon_id,
                    "exon_id": row.exon_id,
                    "branch_type": row.branch_type,
                    "coord": row.coord,
                    "gene": row.gene if pd.notna(row.gene) else "",
                    "TE_info": row.TE_info if pd.notna(row.TE_info) else "",
                }
            )
    return selected


def _resolve_te_exon_usage_table(quant_dir, project):
    quant_folder = os.path.join(quant_dir, QUANT_DIR)
    expected = os.path.join(quant_folder, f"{project}.TE_overlap.exon_usage.tsv")
    if os.path.isfile(expected):
        return expected
    all_tables = sorted(glob.glob(os.path.join(quant_folder, "*.TE_overlap.exon_usage.tsv")))
    if len(all_tables) == 1:
        return all_tables[0]
    if len(all_tables) > 1:
        raise RuntimeError(
            f"Multiple TE-overlap exon usage tables found in {quant_folder}; please pass --project matching quant output."
        )
    raise FileNotFoundError(
        f"ncPred requires TE-overlap exon usage table from quant output: {expected}"
    )


def _load_event_supporting_transcripts(quant_dir, project, event_ids):
    usage_tx_path = _resolve_te_exon_usage_table(quant_dir, project)

    df = pd.read_csv(usage_tx_path, sep="\t", dtype=str).fillna("")
    required = {"exon_id", "transcript_id"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise RuntimeError(
            f"Invalid event-transcript table {usage_tx_path}; missing column(s): {', '.join(missing)}"
        )

    wanted = {str(x) for x in event_ids}
    event_to_txs = {event_id: set() for event_id in wanted}
    for _, row in df.iterrows():
        event_id = str(row["exon_id"])
        if event_id not in wanted:
            continue
        txs = [x.strip() for x in str(row["transcript_id"]).split(",") if x.strip()]
        event_to_txs[event_id].update(txs)

    return {event_id: sorted(txs) for event_id, txs in event_to_txs.items()}, usage_tx_path


def _select_sig_transcripts_from_support(sig_events, event_to_txs, transcript_data):
    selected = defaultdict(lambda: {"gene": "", "events": []})
    missing_events = []
    missing_transcripts = set()

    for row in sig_events.itertuples(index=False):
        event_id = str(row.exon_id)
        txs = event_to_txs.get(event_id, [])
        if not txs:
            missing_events.append(event_id)
            continue
        event_record = {
            "event_id": event_id,
            "exon_id": event_id,
            "branch_type": row.branch_type,
            "coord": row.coord,
            "gene": row.gene if pd.notna(row.gene) else "",
            "TE_info": row.TE_info if pd.notna(row.TE_info) else "",
        }
        for tid in txs:
            if tid not in transcript_data:
                missing_transcripts.add(tid)
                continue
            selected[tid]["gene"] = transcript_data[tid].get("gene", "")
            selected[tid]["events"].append(event_record)

    return selected, missing_events, missing_transcripts


def _merge_selected_transcript_maps(primary, fallback):
    merged = defaultdict(lambda: {"gene": "", "events": []})
    for source in (primary, fallback):
        for tid, payload in source.items():
            merged[tid]["gene"] = merged[tid]["gene"] or payload.get("gene", "")
            merged[tid]["events"].extend(payload.get("events", []))
    return merged


def _write_sig_gtf(gtf_input, gtf_output, selected_transcripts):
    with open(gtf_input, "r", encoding="utf-8") as src, open(gtf_output, "w", encoding="utf-8") as dst:
        for line in src:
            marker = 'transcript_id "'
            if marker not in line:
                continue
            tid = line.split(marker, 1)[1].split('"', 1)[0]
            if tid in selected_transcripts:
                dst.write(line)


def _build_plek_summary(transcripts_in_order, plek_result_file, selected_map, output_file):
    with open(plek_result_file, "r", encoding="utf-8") as fh:
        plek_lines = [line.strip() for line in fh if line.strip()]
    if len(plek_lines) != len(transcripts_in_order):
        raise ValueError("Mismatch between PLEK results and selected transcript list")

    records = []
    for tid, pred in zip(transcripts_in_order, plek_lines):
        evt_list = selected_map.get(tid, {}).get("events", [])
        genes = {selected_map.get(tid, {}).get("gene", "")}
        genes.update({e["gene"] for e in evt_list if e["gene"]})
        genes = {g for g in genes if g}

        diff_exon = sorted({str(e.get("branch_type", "other")) for e in evt_list})
        meta_exon = sorted({e["coord"] for e in evt_list})
        te_set = set()
        for e in evt_list:
            for te in str(e["TE_info"]).split(","):
                te = te.strip()
                if te:
                    te_set.add(te)

        records.append(
            {
                "Transcript": tid,
                "Prediction": pred,
                "Gene": ",".join(sorted(genes)),
                "DiffExon": ",".join(diff_exon),
                "Exon": ",".join(meta_exon),
                "MetaExon": ",".join(meta_exon),
                "TE": ",".join(sorted(te_set)),
            }
        )

    pd.DataFrame(records).to_csv(output_file, index=False)


def _count_fasta_records(fasta_path):
    n = 0
    with open(fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def _run_external_command(args, title, cmd, cwd=None):
    extend_log = getattr(args, "diff_extend_log", None)
    if extend_log:
        os.makedirs(os.path.dirname(extend_log), exist_ok=True)
        with open(extend_log, "a", encoding="utf-8") as handle:
            handle.write(f"========={title}=========\n")
            handle.write("command: " + " ".join(str(part) for part in cmd) + "\n")
            if cwd:
                handle.write(f"cwd: {cwd}\n")
            handle.write("\n")
            result = sp.run(cmd, cwd=cwd, stdout=handle, stderr=handle)
        if result.returncode != 0:
            raise RuntimeError(f"{title} failed with exit code {result.returncode}; see {extend_log}")
        return

    result = sp.run(cmd, cwd=cwd, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    if result.returncode != 0:
        message = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(f"{title} failed with exit code {result.returncode}: {message}")


def _build_plek_input_with_singleton_guard(fasta_path, out_dir):
    """
    PLEK2_model_v3 crashes on singleton FASTA due numpy shape edge case.
    For one-record input, duplicate the record as a dummy and trim predictions later.
    Returns (input_fasta_for_plek, real_record_count).
    """
    real_n = _count_fasta_records(fasta_path)
    if real_n <= 0:
        raise RuntimeError(f"No FASTA records found: {fasta_path}")
    if real_n > 1:
        return fasta_path, real_n

    guarded_path = os.path.join(out_dir, "selected_transcripts.plek_guarded.fasta")
    with open(fasta_path, "r", encoding="utf-8") as src:
        content = src.read().rstrip("\n")
    if not content.startswith(">"):
        raise RuntimeError(f"Invalid FASTA format: {fasta_path}")
    # Duplicate singleton record with a dummy id to avoid PLEK singleton bug.
    with open(guarded_path, "w", encoding="utf-8") as dst:
        dst.write(content + "\n")
        lines = content.splitlines()
        seq = "".join([ln.strip() for ln in lines[1:] if ln.strip()])
        dst.write(">__TEXTRA_PLEK_DUMMY__\n")
        dst.write(seq + "\n")
    return guarded_path, real_n


def _trim_plek_results(raw_result_path, out_path, keep_n):
    with open(raw_result_path, "r", encoding="utf-8") as fh:
        lines = [line.rstrip("\n") for line in fh if line.strip()]
    if len(lines) < keep_n:
        raise RuntimeError(
            f"PLEK output rows ({len(lines)}) smaller than expected real transcripts ({keep_n})."
        )
    with open(out_path, "w", encoding="utf-8") as fh:
        for line in lines[:keep_n]:
            fh.write(line + "\n")


def _write_selected_transcript_detail(selected_transcripts, selected_map, transcript_data, out_path):
    rows = []
    for tid in selected_transcripts:
        events = selected_map.get(tid, {}).get("events", [])
        genes = {selected_map.get(tid, {}).get("gene", "")}
        genes.update({str(event.get("gene", "")) for event in events if str(event.get("gene", "")).strip()})
        genes = sorted(gene for gene in genes if gene)
        rows.append(
            {
                "transcript_id": tid,
                "gene_id": ",".join(genes),
                "transcript_length": _transcript_length(transcript_data.get(tid, {})),
                "event_n": len(events),
                "event_ids": ",".join(sorted({str(event.get("event_id", "")) for event in events if str(event.get("event_id", "")).strip()})),
                "event_coords": ",".join(sorted({str(event.get("coord", "")) for event in events if str(event.get("coord", "")).strip()})),
                "TE_info": ",".join(sorted({str(event.get("TE_info", "")) for event in events if str(event.get("TE_info", "")).strip()})),
            }
        )
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)


def _transcript_length(transcript_info):
    return int(sum(max(0, int(end) - int(start)) for start, end in transcript_info.get("exons", [])))


def _resolve_diff_result_for_ncpred(de_out):
    candidates = [
        ("differential_significant_usage.tsv", True),
    ]
    for filename, preselected in candidates:
        path = os.path.join(de_out, filename)
        if os.path.isfile(path):
            return path, preselected
    raise FileNotFoundError(
        f"No supported differential result found in {de_out}. "
        "Expected differential_significant_usage.tsv."
    )


def ncPred_func(args):
    """Run coding-potential prediction on significant AFE/ALE-linked transcripts."""
    padj = args.padj
    prep_dir = args.prep
    genome_fasta = args.genome
    model = args.ncpred_model
    min_length = int(getattr(args, "min_length", 200))
    total_steps = int(getattr(args, "diff_total_steps", 2))

    ncPred_out = os.path.join(args.out_dir, "ncPred")
    os.makedirs(ncPred_out, exist_ok=True)
    log_message(
        "[INFO]",
        f"Step 2/{total_steps}: Coding potential prediction",
        bold=True,
        color="step",
    )
    log_message(
        "[INFO]",
        f"ncPred parameters: model={model}, min_length>={min_length}, genome={os.path.abspath(genome_fasta)}",
    )
    stale_outputs = [
        "selected_transcripts.detail.tsv",
        "selected_transcripts.fasta",
        "selected_transcripts.plek_guarded.fasta",
        "sig_transcripts.txt",
        "sig_transcripts.gtf",
        "sig_transcripts.fasta",
        "sig_transcripts.plek_guarded.fasta",
        "plek_final_result.csv",
        "results.trimmed",
    ]
    for name in stale_outputs:
        path = os.path.join(ncPred_out, name)
        if os.path.isfile(path):
            os.remove(path)

    DE_out = os.path.join(args.out_dir, "DE")
    diff_out, preselected = _resolve_diff_result_for_ncpred(DE_out)
    sig_events = _load_sig_events(diff_out, padj, preselected=preselected)
    log_message(
        "[INFO]",
        f"ncPred input: significant_events={len(sig_events):,}, source={os.path.abspath(diff_out)}",
    )
    if sig_events.empty:
        log_message("[WARNING]", "No significant differential TE-overlap exon events found for ncPred.", color="warning")
        return None

    gtf_input = resolve_consensus_gtf(prep_dir)
    if not os.path.isfile(gtf_input):
        raise FileNotFoundError(f"consensus transcript GTF not found: {gtf_input}")
    transcript_data, exon_index = _load_transcript_exons(gtf_input)
    event_to_txs, usage_tx_path = _load_event_supporting_transcripts(
        args.quant,
        args.project,
        sig_events["exon_id"].astype(str),
    )
    selected_map, missing_events, missing_transcripts = _select_sig_transcripts_from_support(
        sig_events=sig_events,
        event_to_txs=event_to_txs,
        transcript_data=transcript_data,
    )
    if missing_transcripts:
        log_message(
            "[WARNING]",
            f"{len(missing_transcripts)} supporting transcript(s) from {usage_tx_path} were not found in consensus GTF.",
            color="warning",
        )
    if missing_events:
        log_message(
            "[WARNING]",
            f"{len(missing_events)} significant event(s) lacked supporting transcripts in {usage_tx_path}; "
            "falling back to coordinate overlap for those events.",
            color="warning",
        )
        fallback_events = sig_events[sig_events["exon_id"].astype(str).isin(set(missing_events))].copy()
        selected_map = _merge_selected_transcript_maps(
            selected_map,
            _select_sig_transcripts(fallback_events, exon_index),
        )

    selected_transcripts = [
        tid
        for tid in sorted(selected_map.keys())
        if _transcript_length(transcript_data.get(tid, {})) >= min_length
    ]
    selected_map = {tid: selected_map[tid] for tid in selected_transcripts}
    if not selected_transcripts:
        log_message("[WARNING]", "No transcripts overlapped significant exon coordinates.", color="warning")
        return None
    log_message(
        "[INFO]",
        f"Selected transcripts for PLEK2: transcripts={len(selected_transcripts):,}",
    )

    if getattr(args, "detail", False) or getattr(args, "debug", False):
        detail_out = os.path.join(ncPred_out, "selected_transcripts.detail.tsv")
        _write_selected_transcript_detail(selected_transcripts, selected_map, transcript_data, detail_out)
        log_message("[SUCCESS]", f"Selected transcript detail: {os.path.abspath(detail_out)}", color="success")

    gtf_output = os.path.join(ncPred_out, "selected_transcripts.gtf")
    _write_sig_gtf(gtf_input, gtf_output, set(selected_transcripts))

    if shutil.which("gffread") is None:
        raise RuntimeError("gffread not found in PATH; required for --ncpred.")
    fasta_output = os.path.join(ncPred_out, "selected_transcripts.fasta")
    _run_external_command(args, "NCPRED-GFFREAD", ["gffread", gtf_output, "-g", genome_fasta, "-w", fasta_output])

    plek_path = getattr(args, "plek2_path", None)
    if not plek_path:
        plek_path, candidates = resolve_external_file("PLEK*", "PLEK2.py", anchor_file=__file__)
        if not plek_path:
            searched = ", ".join(candidates) if candidates else os.path.abspath(os.path.join("util", "external"))
            raise FileNotFoundError(
                "PLEK2.py not found for --ncpred. "
                "Expected PLEK2 under TEXTRA_EXTERNAL_DIR/PLEK*/PLEK2.py or util/external/PLEK*/PLEK2.py. "
                f"Searched external roots: {searched}"
            )
    plek_path = os.path.abspath(plek_path)
    plek_dir = os.path.dirname(plek_path)
    log_message("[INFO]", f"PLEK2 path: {plek_path}")

    fasta_output = os.path.abspath(fasta_output)
    plek_input_fasta, real_n = _build_plek_input_with_singleton_guard(fasta_output, ncPred_out)
    if os.path.abspath(plek_input_fasta) != fasta_output:
        log_message(
            "[INFO]",
            "Applied singleton PLEK guard by adding one dummy transcript to avoid PLEK2 shape bug.",
            color="step",
        )
    _run_external_command(args, "NCPRED-PLEK2", ["python", plek_path, "-i", plek_input_fasta, "-m", model], cwd=plek_dir)

    plek_result_raw = os.path.join(plek_dir, "results")
    plek_result = plek_result_raw
    if real_n == 1:
        plek_result = os.path.join(ncPred_out, "results.trimmed")
        _trim_plek_results(plek_result_raw, plek_result, keep_n=real_n)
    plek_final = os.path.join(ncPred_out, "plek_result.csv")
    _build_plek_summary(selected_transcripts, plek_result, selected_map, plek_final)
    log_message(
        "[INFO]",
        f"PLEK2 prediction completed: transcripts={len(selected_transcripts):,}",
    )
    if plek_result != plek_result_raw and os.path.isfile(plek_result):
        os.remove(plek_result)
    guarded_fasta = os.path.join(ncPred_out, "selected_transcripts.plek_guarded.fasta")
    if os.path.isfile(guarded_fasta):
        os.remove(guarded_fasta)
    if not getattr(args, "debug", False):
        if os.path.isfile(fasta_output):
            os.remove(fasta_output)
    if os.path.isfile(plek_result_raw):
        os.remove(plek_result_raw)
    stale_outputs = [
        "plek_final_result.csv",
        "sig_transcripts.gtf",
        "sig_transcripts.txt",
        "sig_transcripts.fasta",
    ]
    if not (getattr(args, "detail", False) or getattr(args, "debug", False)):
        stale_outputs.append("selected_transcripts.detail.tsv")
    if not getattr(args, "debug", False):
        stale_outputs.append("selected_transcripts.fasta")
    for name in stale_outputs:
        path = os.path.join(ncPred_out, name)
        if os.path.isfile(path):
            os.remove(path)

    log_message("[SUCCESS]", f"Selected transcripts: {os.path.abspath(gtf_output)}", color="success")
    log_message("[SUCCESS]", f"PLEK2 result: {os.path.abspath(plek_final)}", color="success")

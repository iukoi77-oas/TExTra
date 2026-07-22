"""Resolve quant command execution, reference checks, and read inputs."""

import os
import shutil
import subprocess as sp

import pandas as pd

from util.common.define_layout import ALIGNMENT_DIR, ASSEMBLY_DIR, TRANSCRIPT_GENE_ASSIGNMENT_TSV
from util.common.write_logs import log_message


def run_command(cmd, show_tool_output, log_path=None):
    """Run an external command, optionally capturing stdout/stderr to a log file."""
    if show_tool_output:
        sp.run(cmd, check=True)
        return
    try:
        result = sp.run(cmd, check=True, stdout=sp.PIPE, stderr=sp.PIPE, text=False)
    except sp.CalledProcessError as exc:
        if log_path:
            _write_command_log(log_path, cmd, exc.stdout, exc.stderr, returncode=exc.returncode)
        raise
    if log_path:
        _write_command_log(log_path, cmd, result.stdout, result.stderr)


def _write_command_log(log_path, cmd, stdout, stderr, returncode=None):
    stdout_text = stdout.decode("utf-8", errors="replace") if isinstance(stdout, (bytes, bytearray)) else str(stdout or "")
    stderr_text = stderr.decode("utf-8", errors="replace") if isinstance(stderr, (bytes, bytearray)) else str(stderr or "")
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "a", encoding="utf-8") as fh:
        fh.write("$ " + " ".join(str(x) for x in cmd) + "\n")
        if returncode is not None:
            fh.write(f"[returncode] {returncode}\n")
        if stdout_text:
            fh.write(stdout_text)
            if not stdout_text.endswith("\n"):
                fh.write("\n")
        if stderr_text:
            fh.write(stderr_text)
            if not stderr_text.endswith("\n"):
                fh.write("\n")


def is_valid_rsem_reference_prefix(ref_prefix):
    """Return True when an RSEM reference prefix has the expected index files."""
    if not ref_prefix:
        return False
    need_files = [
        f"{ref_prefix}.grp",
        f"{ref_prefix}.transcripts.fa",
        f"{ref_prefix}.seq",
    ]
    return all(os.path.isfile(p) for p in need_files)


def infer_rsem_strandedness(strand, rsem_strandedness="auto"):
    """Resolve RSEM strandedness from explicit option or shared strand mode."""
    v = str(rsem_strandedness).strip().lower()
    if v in {"none", "forward", "reverse"}:
        return v
    s = str(strand).strip().lower()
    if s in {"none", "0"}:
        return "none"
    if s in {"fr", "f", "1"}:
        return "forward"
    if s in {"rf", "r", "2"}:
        return "reverse"
    return "none"


def looks_like_fasta_path(path):
    """Return True when path text has a common FASTA extension."""
    p = str(path).strip().lower()
    return p.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"))


def looks_like_gzip_path(path):
    """Return True when path text has a gzip extension."""
    return str(path).strip().lower().endswith(".gz")


def looks_like_bam_path(path):
    """Return True when path text has a BAM extension."""
    return str(path).strip().lower().endswith(".bam")


def infer_star_path(rsem_star_path=None):
    """Resolve a STAR executable or directory path for RSEM STAR mode."""
    if rsem_star_path:
        path = os.path.abspath(str(rsem_star_path))
        if os.path.isfile(path):
            return os.path.dirname(path)
        if os.path.isdir(path):
            return path
        raise FileNotFoundError(f"Invalid RSEM STAR path: {rsem_star_path}")
    star_bin = shutil.which("STAR")
    if star_bin:
        return os.path.dirname(star_bin)
    return None


def resolve_rsem_aligner(aligner, has_star=False, has_bowtie2=False):
    """Resolve the RSEM aligner mode without changing existing fallback order."""
    mode = str(aligner).strip().lower()
    if mode in {"star", "bowtie2", "bowtie"}:
        return mode
    if bool(has_star):
        return "star"
    if bool(has_bowtie2):
        return "bowtie2"
    return "bowtie"


def pick_existing_col(df, candidates):
    """Return the first DataFrame column matching candidate names case-insensitively."""
    cols = {str(c).lower(): str(c) for c in df.columns}
    for c in candidates:
        k = str(c).lower()
        if k in cols:
            return cols[k]
    return None


def load_transcript_gene_map(args, transcript_gtf_path):
    """Load transcript-to-gene assignments produced by prep."""
    assignment_tsv = os.path.join(args.prep, ASSEMBLY_DIR, TRANSCRIPT_GENE_ASSIGNMENT_TSV)
    if not os.path.isfile(assignment_tsv):
        raise FileNotFoundError(f"required prep transcript-gene assignment table missing: {assignment_tsv}")
    df = pd.read_csv(assignment_tsv, sep="\t")
    required = {"transcript_id", "assigned_gene_id"}
    if not required.issubset(set(df.columns)):
        raise RuntimeError(
            f"invalid transcript-gene assignment table {assignment_tsv}; required columns: {sorted(required)}"
        )
    tx = df["transcript_id"].astype(str)
    gene = df["assigned_gene_id"].fillna("").astype(str)
    return dict(zip(tx, gene))


def load_prep_input_reads(input_path, sample_list, readtype):
    """Load quant reads directly from the original prep input TSV."""
    if not input_path:
        return {}
    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"input TSV not found: {input_path}")

    df = pd.read_csv(input_path, sep="\t", header=None, dtype=str).fillna("")
    if df.empty or df.shape[1] < 2:
        raise RuntimeError(f"invalid quant --input TSV: {input_path}")
    out = {}
    sample_set = set(str(sample).strip() for sample in sample_list if str(sample).strip())
    seen_samples = set()
    readtype = str(readtype).strip().lower()
    for _, row in df.iterrows():
        sample = str(row.iloc[0]).strip()
        if not sample or sample.lower() in {"na", "nan"}:
            continue
        if sample_set and sample not in sample_set:
            continue
        seen_samples.add(sample)
        rep_idx = 0
        for cell in row.iloc[1:].tolist():
            token = str(cell).strip()
            if not token or token.lower() in {"na", "nan", "-"}:
                continue
            rep_idx += 1
            quant_unit = f"{sample}_rep{rep_idx}"
            if "," in token:
                r1, r2 = [part.strip() for part in token.split(",", 1)]
                if looks_like_bam_path(r1) or looks_like_bam_path(r2):
                    bam_path = r1 if looks_like_bam_path(r1) else r2
                    out[quant_unit] = {"bam": bam_path}
                elif readtype == "paired":
                    if not r1 or not r2:
                        raise RuntimeError(f"paired quant --input entry requires read1,read2 for {quant_unit}")
                    out[quant_unit] = {"read1": r1, "read2": r2}
                else:
                    out[quant_unit] = {"read1": None, "read2": None, "single": [r1]}
            else:
                if looks_like_bam_path(token):
                    out[quant_unit] = {"bam": token}
                elif readtype == "paired":
                    raise RuntimeError(
                        f"paired quant --input entry requires read1,read2 for FASTQ replicate {quant_unit}"
                    )
                else:
                    out[quant_unit] = {"read1": None, "read2": None, "single": [token]}
    missing = sorted(sample_set - seen_samples)
    if missing:
        raise RuntimeError(f"quant --input missing sample(s): {', '.join(missing)}")
    if not out:
        raise RuntimeError(f"quant --input produced no read records for requested samples: {input_path}")
    return out


def convert_bam_to_reads_for_rsem(args, sample, bam_path, tag=None):
    """Convert an accepted-hit BAM to FASTQ inputs for reads-mode quantification."""
    samtools = shutil.which("samtools")
    if not samtools:
        raise RuntimeError("samtools not found in PATH; required for BAM->FASTQ conversion.")
    if not os.path.isfile(bam_path):
        raise FileNotFoundError(f"BAM file not found for sample '{sample}': {bam_path}")

    root = args.rsem_bam2fastq_dir or os.path.join(args.out_dir, "tmp", "rsem_bam2fastq")
    sample_dir = os.path.join(root, str(sample))
    os.makedirs(sample_dir, exist_ok=True)

    sort_threads = str(max(1, int(args.threads)))
    if tag is None:
        base = os.path.basename(str(bam_path))
        tag = os.path.splitext(base)[0]
    safe_tag = str(tag).replace("/", "_")
    sorted_bam = os.path.join(sample_dir, f"{safe_tag}.name_sorted.bam")
    conv_log = os.path.join(sample_dir, f"{sample}.bam2fastq.log")

    need_sort = bool(args.force) or (not os.path.isfile(sorted_bam))
    if need_sort:
        sort_cmd = [samtools, "sort", "-n", "-@", sort_threads, "-o", sorted_bam, os.path.abspath(bam_path)]
        run_command(sort_cmd, show_tool_output=bool(args.show_tool_output), log_path=conv_log)

    readtype = str(args.readtype).strip().lower()
    if readtype == "paired":
        r1 = os.path.join(sample_dir, f"{safe_tag}_R1.fastq.gz")
        r2 = os.path.join(sample_dir, f"{safe_tag}_R2.fastq.gz")
        need_fastq = bool(args.force) or (not os.path.isfile(r1)) or (not os.path.isfile(r2))
        if need_fastq:
            fastq_cmd = [
                samtools,
                "fastq",
                "-@",
                sort_threads,
                "-n",
                "-F",
                "0x900",
                "-1",
                r1,
                "-2",
                r2,
                "-0",
                "/dev/null",
                "-s",
                "/dev/null",
                sorted_bam,
            ]
            run_command(fastq_cmd, show_tool_output=bool(args.show_tool_output), log_path=conv_log)
        return {"read1": r1, "read2": r2, "single": []}

    single_fastq = os.path.join(sample_dir, f"{safe_tag}.single.fastq.gz")
    need_fastq = bool(args.force) or (not os.path.isfile(single_fastq))
    if need_fastq:
        fastq_cmd = [
            samtools,
            "fastq",
            "-@",
            sort_threads,
            "-n",
            "-F",
            "0x900",
            "-0",
            single_fastq,
            "-1",
            "/dev/null",
            "-2",
            "/dev/null",
            "-s",
            "/dev/null",
            sorted_bam,
        ]
        run_command(fastq_cmd, show_tool_output=bool(args.show_tool_output), log_path=conv_log)
    return {"read1": None, "read2": None, "single": [single_fastq]}


def auto_resolve_reads_from_bam(args, sample_list, bamfiles_dict):
    """Resolve quant reads by converting prep BAM outputs when fallback is enabled."""
    if not bool(args.rsem_auto_bam2fastq):
        return None
    if not bamfiles_dict:
        return None

    resolved = {}
    missing = []
    total_bams = sum(len(set(bamfiles_dict.get(sample, []))) for sample in sample_list)
    if total_bams:
        log_message(
            "[INFO]",
            f"Auto BAM->FASTQ fallback: samples={len(sample_list)}, BAMs={total_bams}",
            color="info",
        )
    for sample in sample_list:
        bams = sorted(set(bamfiles_dict.get(sample, [])))
        if not bams:
            missing.append(sample)
            continue
        log_message(
            "[INFO]",
            f"Auto BAM->FASTQ for sample '{sample}' from {len(bams)} BAM(s).",
            color="info",
            console=bool(getattr(args, "debug", False)),
        )
        seen_names = set(resolved.keys())
        for bam_path in bams:
            base_name = os.path.basename(str(bam_path))
            rep_name = os.path.splitext(base_name)[0]
            if rep_name.endswith("_accepted_hits"):
                rep_name = rep_name[: -len("_accepted_hits")]
            if not rep_name:
                rep_name = str(sample)
            uniq_name = str(rep_name)
            k = 2
            while uniq_name in seen_names:
                uniq_name = f"{rep_name}.bam{k}"
                k += 1
            seen_names.add(uniq_name)
            log_message(
                "[INFO]",
                f"  using BAM as sample '{uniq_name}': {bam_path}",
                color="info",
                console=bool(getattr(args, "debug", False)),
            )
            resolved[uniq_name] = convert_bam_to_reads_for_rsem(
                args, sample, bam_path, tag=uniq_name
            )
            log_message(
                "[INFO]",
                f"  BAM sample ready: {uniq_name}",
                color="info",
                console=bool(getattr(args, "debug", False)),
            )
    if missing:
        alignment_dir = os.path.join(args.prep, ALIGNMENT_DIR)
        raise RuntimeError(
            "BAM fallback enabled but no *_accepted_hits.bam found for sample(s): "
            f"{', '.join(missing)}. "
            f"Please check sample names in --samples and BAM layout under: {alignment_dir}"
        )
    return resolved


def resolve_sample_reads(args, sample_list, bamfiles_dict=None):
    """Resolve per-sample read inputs from prep input TSV, direct args, or BAM fallback."""
    input_map = load_prep_input_reads(getattr(args, "input", None), sample_list, args.readtype)
    if input_map:
        resolved = {}
        for quant_unit, entry in input_map.items():
            if entry.get("bam"):
                resolved[quant_unit] = convert_bam_to_reads_for_rsem(args, quant_unit, entry["bam"])
            elif str(args.readtype).lower() == "paired":
                r1 = str(entry.get("read1") or "").strip()
                r2 = str(entry.get("read2") or "").strip()
                if looks_like_bam_path(r1) or looks_like_bam_path(r2):
                    bam_path = r1 if looks_like_bam_path(r1) else r2
                    resolved[quant_unit] = convert_bam_to_reads_for_rsem(args, quant_unit, bam_path)
                else:
                    resolved[quant_unit] = entry
            else:
                singles = [str(x).strip() for x in list(entry.get("single") or []) if str(x).strip()]
                if len(singles) == 1 and looks_like_bam_path(singles[0]):
                    resolved[quant_unit] = convert_bam_to_reads_for_rsem(args, quant_unit, singles[0])
                else:
                    resolved[quant_unit] = entry
        return resolved

    auto_reads = auto_resolve_reads_from_bam(args, sample_list, bamfiles_dict or {})
    if auto_reads:
        return auto_reads

    if str(args.readtype).lower() == "paired":
        raise RuntimeError(
            "No read inputs provided. Supply --input, or ensure prep alignment BAMs exist for the default BAM fallback."
        )
    raise RuntimeError(
        "No read inputs provided. Supply --input, or ensure prep alignment BAMs exist for the default BAM fallback."
    )

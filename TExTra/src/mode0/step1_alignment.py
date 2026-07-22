"""Read alignment and BAM filtering helpers for TExTra prep."""

import os
import shutil

from joblib import Parallel, delayed
from util.common.write_logs import log_message
from util.common.define_layout import ALIGNMENT_DIR
from util.prep.run_commands import run_prep_command


STAR_SJDB_OVERHANG = 100
STAR_CHIM_SEGMENT_MIN = 10
STAR_BAM_SORTING_BINS = 200


def _is_gzipped_fastq(path):
    value = str(path or "").strip().lower()
    return value.endswith((".fq.gz", ".fastq.gz"))


def _is_bam_like(path):
    value = str(path or "").strip().lower()
    return value.endswith(".bam")


def align_func(input_info, args):
    """Run STAR alignment (or sort existing BAM inputs) and filter mapped reads."""
    log_message("[INFO]", "Step 1/3: Alignment", bold=True, color="step")

    alignment_dir = os.path.join(args.out_dir, ALIGNMENT_DIR)
    os.makedirs(alignment_dir, exist_ok=True)
    log_message("[INFO]", f"Created alignment directory at {alignment_dir}", color="info", console=getattr(args, "debug", False))

    input_info['mate1'] = input_info['mate1'].astype(str).str.strip()
    input_info['mate2'] = input_info['mate2'].astype(str).str.strip()

    log_message("[INFO]", "Resolving input replicates.", color="info", console=getattr(args, "debug", False))
    rows = input_info[['mate1', 'mate2', 'sample']].to_dict('records')
    bam_rows = []
    fastq_rows = []
    for row in rows:
        mate1 = str(row['mate1']).strip()
        mate2 = str(row['mate2']).strip()
        if _is_bam_like(mate1):
            bam_rows.append((mate1, row['sample']))
        else:
            fastq_rows.append(row)

    max_jobs = int(getattr(args, "njobs", 0) or args.threads)
    auto_index_dir = None
    if fastq_rows:
        if not args.index:
            args.index = os.path.join(alignment_dir, 'index')
            auto_index_dir = args.index
        os.makedirs(args.index, exist_ok=True)
        if not os.path.exists(os.path.join(args.index, 'SAindex')):
            log_message("[INFO]", "STAR index: generating", color="info")
            run_prep_command(
                [
                    "STAR",
                    "--runThreadN",
                    str(args.threads),
                    "--runMode",
                    "genomeGenerate",
                    "--genomeDir",
                    args.index,
                    "--genomeFastaFiles",
                    args.genome,
                    "--sjdbGTFfile",
                    args.gene,
                    "--sjdbOverhang",
                    str(getattr(args, "sjdbOverhang", STAR_SJDB_OVERHANG)),
                ],
                args,
                "STAR-index",
            )
            log_message("[INFO]", "STAR index: generated")
        else:
            log_message("[INFO]", "STAR index: existing", color="info")

    sorted_existing_bams = []
    if bam_rows:
        log_message("[INFO]", f"Sorting existing BAM inputs: replicates={len(bam_rows)}", color="info")
        sort_jobs = min(max_jobs, max(1, len(bam_rows)))
        sort_threads = max(1, args.threads // sort_jobs)
        sorted_existing_bams = Parallel(n_jobs=sort_jobs)(
            delayed(sort_bam)(bam, sample, alignment_dir, args, sort_threads) for bam, sample in bam_rows
        )

    aligned_bams = []
    if fastq_rows:
        align_jobs = min(max_jobs, max(1, len(fastq_rows)))
        star_threads = max(1, args.threads // align_jobs)
        log_message(
            "[INFO]",
            f"Aligning FASTQ replicates: replicates={len(fastq_rows)}, jobs={align_jobs}, threads/job={star_threads}",
            color="info",
        )
        aligned_bams = Parallel(n_jobs=align_jobs)(
            delayed(run_star)(row, alignment_dir, args, star_threads) for row in fastq_rows
        )

    bam_path_by_sample = {}
    for bam_path, sample in zip(sorted_existing_bams, [sample for _bam, sample in bam_rows]):
        bam_path_by_sample[str(sample)] = bam_path
    for bam_path, row in zip(aligned_bams, fastq_rows):
        bam_path_by_sample[str(row['sample'])] = bam_path

    sample_names = input_info['sample'].tolist()
    bam_files = [bam_path_by_sample[str(sample)] for sample in sample_names]
    log_message("[INFO]", f"Alignment completed: sorted BAMs={len(bam_files)}")
    output_bams = process_bam(bam_files, sample_names, args)

    if not args.debug:
        cleanup_counts = delete_temp_files(bam_files, alignment_dir, input_info['sample'], args)
        if auto_index_dir and os.path.isdir(auto_index_dir):
            shutil.rmtree(auto_index_dir, ignore_errors=True)
            cleanup_counts["directories"] += 1
            log_message(
                "[INFO]",
                f"Removed auto-generated STAR index: {os.path.relpath(auto_index_dir, args.out_dir)}",
                color="info",
                console=getattr(args, "detail", False),
            )
        log_message(
            "[INFO]",
            (
                "Cleaned alignment intermediates: "
                f"samples={len(sample_names)}, files={cleanup_counts['files']}, directories={cleanup_counts['directories']}"
            ),
            color="info",
        )

    return output_bams, input_info['sample']


def delete_temp_files(bam_files, alignment_dir, samples, args=None):
    """Delete intermediate files and keep only `{sample}_accepted_hits.bam`."""
    counts = {"files": 0, "directories": 0}
    for bam, sample in zip(bam_files, samples):
        sample_folder = os.path.dirname(bam)
        for file in os.listdir(sample_folder):
            file_path = os.path.join(sample_folder, file)
            if not file.endswith("accepted_hits.bam"):
                try:
                    if os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                        counts["directories"] += 1
                        log_message(
                            "[INFO]",
                            f"Deleted temporary folder: {file_path}",
                            color="info",
                            console=bool(getattr(args, "debug", False)),
                        )
                    else:
                        os.remove(file_path)
                        counts["files"] += 1
                        log_message(
                            "[INFO]",
                            f"Deleted temporary file: {file_path}",
                            color="info",
                            console=bool(getattr(args, "debug", False)),
                        )
                except Exception as e:
                    log_message("[ERROR]", f"Error deleting {file_path}: {str(e)}", color="error")
    return counts


def sort_bam(bam_file, sample, alignment_dir, args, sam_threads=None):
    """Sort BAM for downstream processing in a sample-specific folder."""
    sample_folder = os.path.join(alignment_dir, sample)
    os.makedirs(sample_folder, exist_ok=True)
    sorted_bam = os.path.join(sample_folder, "Aligned.sortedByCoord.out.bam")
    if sam_threads is None:
        sam_threads = max(1, int(getattr(args, 'threads', 1)))
    run_prep_command(
        ["samtools", "sort", "-@", str(sam_threads), "-o", sorted_bam, bam_file],
        args,
        f"samtools_sort-{sample}",
    )
    
    return sorted_bam


def run_star(row, alignment_dir, args, star_threads=None):
    """Run STAR for one sample and return sorted BAM path."""
    mate1 = row['mate1']
    mate2 = row['mate2']
    sample = row['sample']
    if star_threads is None:
        star_threads = max(1, int(getattr(args, 'threads', 1)))
    
    mate2 = '' if mate2 == '-' else mate2
    output_prefix = os.path.join(alignment_dir, f"{sample}/")
    os.makedirs(output_prefix, exist_ok=True)

    read_inputs = [mate1]
    if mate2:
        read_inputs.append(mate2)

    star_cmd = [
        "STAR",
        "--genomeDir",
        args.index,
        "--runThreadN",
        str(star_threads),
        "--twopassMode",
        "Basic",
        "--sjdbGTFfile",
        args.gene,
    ]
    if all(_is_gzipped_fastq(p) for p in read_inputs):
        star_cmd.extend(["--readFilesCommand", "zcat"])
    star_cmd.extend(
        [
            "--readFilesIn",
            mate1,
            *([mate2] if mate2 else []),
            "--outSAMtype",
            "BAM",
            "SortedByCoordinate",
            "--outFileNamePrefix",
            output_prefix,
            "--chimOutType",
            "Junctions",
            "--chimSegmentMin",
            str(STAR_CHIM_SEGMENT_MIN),
            "--outBAMsortingBinsN",
            str(STAR_BAM_SORTING_BINS),
            "--seedMultimapNmax",
            str(getattr(args, 'seedMultimapNmax', 50000)),
            "--winAnchorMultimapNmax",
            str(getattr(args, 'winAnchorMultimapNmax', 100)),
            "--outFilterMultimapNmax",
            str(getattr(args, 'outFilterMultimapNmax', 100)),
        ]
    )
    run_prep_command(star_cmd, args, f"STAR-{sample}")
    
    return os.path.join(output_prefix, 'Aligned.sortedByCoord.out.bam')


def process_bam(bam_files, samples, args):
    """Filter unmapped reads from BAM files."""
    log_message("[INFO]", "Filtering mapped reads.", color="info")
    max_jobs = int(getattr(args, "njobs", 0) or args.threads)
    filter_jobs = min(max_jobs, max(1, len(samples)))
    filter_threads = max(1, args.threads // filter_jobs)
    output_bams = Parallel(n_jobs=filter_jobs)(
        delayed(filter_bam)(bam, sample, args, filter_threads) for bam, sample in zip(bam_files, samples)
    )
    log_message("[INFO]", f"Accepted-hit BAMs: {len(output_bams)}")
    return output_bams


def filter_bam(bam_file, sample, args, sam_threads=None):
    """Generate `{sample}_accepted_hits.bam` from the aligned BAM."""
    bam_dir = os.path.dirname(bam_file)
    output_bam = os.path.join(bam_dir, f"{sample}_accepted_hits.bam")
    if sam_threads is None:
        sam_threads = max(1, int(getattr(args, 'threads', 1)))
    with open(output_bam, 'wb') as out_fh:
        run_prep_command(
            ["samtools", "view", "-@", str(sam_threads), "-b", "-F", "4", bam_file],
            args,
            f"samtools_filter-{sample}",
            stdout=out_fh,
        )
    return output_bam

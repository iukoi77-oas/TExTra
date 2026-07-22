"""Transcript assembly and consensus GTF construction for TExTra prep."""

import subprocess
import os
import glob
import shutil

from util.common.write_logs import log_message
from util.prep.map_transcripts import modify_gtf_with_mapping
from util.prep.run_commands import run_prep_command
from util.common.define_layout import ASSEMBLY_DIR, CONSENSUS_GTF, TRANSCRIPT_GENE_ASSIGNMENT_TSV
from util.common.external_tools import resolve_external_dir
from joblib import Parallel, delayed


STRINGTIE_MULTIMAP_FRACTION = 1.0
TRANSCRIPT_GENE_ASSIGNMENT_DETAIL_TSV = "transcript_gene_assignment.detail.tsv"


def _cleanup_assembly_dir(assembly_dir):
    """Keep only final merged transcript annotation in assembly directory."""
    keep_names = {CONSENSUS_GTF, TRANSCRIPT_GENE_ASSIGNMENT_TSV, TRANSCRIPT_GENE_ASSIGNMENT_DETAIL_TSV}
    for entry in os.listdir(assembly_dir):
        if entry in keep_names:
            continue
        path = os.path.join(assembly_dir, entry)
        if os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=True)
        else:
            try:
                os.remove(path)
            except FileNotFoundError:
                pass


def _validate_taco_executable(taco_path, source):
    if os.path.isfile(taco_path) and os.access(taco_path, os.X_OK):
        return taco_path
    message = (
        f"TACO executable taco_run is not executable or not found at {taco_path} ({source}). "
        "Please check the TACO location, or rerun prep with --merge-method stringtie."
    )
    log_message("[ERROR]", message, color="error")
    raise RuntimeError(message)


def _resolve_taco_path(args):
    """Resolve the external TACO executable without importing vendored Python modules."""
    if getattr(args, "taco_path", None):
        candidate = os.path.abspath(args.taco_path)
        taco_path = os.path.join(candidate, "taco_run") if os.path.isdir(candidate) else candidate
        return _validate_taco_executable(taco_path, "--taco-path")

    path_taco = shutil.which("taco_run")
    if path_taco:
        return _validate_taco_executable(path_taco, "PATH")

    external_dir = resolve_external_dir(__file__)
    taco_dirs = sorted(
        d for d in glob.glob(os.path.join(external_dir, 'taco*'))
        if os.path.isdir(d) and os.path.exists(os.path.join(d, 'taco_run'))
    )
    if not taco_dirs:
        message = (
            f"TACO executable taco_run was not found in --taco-path, PATH, or expected external tool directory: {external_dir}. "
            "Please check the TACO location, pass --taco-path, or rerun prep with --merge-method stringtie."
        )
        log_message("[ERROR]", message, color="error")
        raise RuntimeError(message)

    for taco_dir in taco_dirs:
        taco_path = os.path.join(taco_dir, 'taco_run')
        if os.path.exists(taco_path):
            return _validate_taco_executable(taco_path, taco_dir)

    message = (
        "TACO executable taco_run not found in any TACO directory: "
        f"{', '.join(taco_dirs)}. Please check the TACO location, pass --taco-path, "
        "or rerun prep with --merge-method stringtie."
    )
    log_message("[ERROR]", message, color="error")
    raise RuntimeError(message)


def _override_or_default(args, attr_name, default_value):
    value = getattr(args, attr_name, None)
    return default_value if value is None else value


def _has_user_merge_filter(args):
    return any(
        getattr(args, attr_name, None) is not None
        for attr_name in ["min_expr", "min_length", "min_frac"]
    )


def _decode_stderr(error):
    if getattr(error, "stderr", None):
        return error.stderr.decode("utf-8", errors="replace")
    return str(error)


def assemble_func(bamfiles, samples, args):
    """Assemble transcripts with StringTie and merge via StringTie or TACO."""
    log_message("[INFO]", "Step 2/3: Transcript assembly", bold=True, color="step")

    assembly_dir = os.path.join(args.out_dir, ASSEMBLY_DIR)
    if os.path.exists(assembly_dir):
        shutil.rmtree(assembly_dir)
    os.makedirs(assembly_dir)
    log_message("[INFO]", f"Created assembly directory at {assembly_dir}", color="info", console=getattr(args, "debug", False))
    tmp_dir = os.path.join(assembly_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    log_message("[INFO]", f"StringTie assembly: replicates={len(samples)}", color="info")

    if args.strand == 'rf':
        stringtie_cmd = ['stringtie', '--rf']
    elif args.strand == 'fr':
        stringtie_cmd = ['stringtie', '--fr']
    else:
        stringtie_cmd = ['stringtie']

    max_jobs = int(getattr(args, "njobs", 0) or args.threads)
    stringtie_jobs = min(max_jobs, max(1, len(samples)))
    stringtie_threads = max(1, args.threads // stringtie_jobs)
    output_gtf_results = Parallel(n_jobs=stringtie_jobs)(
        delayed(run_stringtie)(bam, sample, stringtie_cmd, tmp_dir, args, stringtie_threads) for bam, sample in zip(bamfiles, samples)
    )
    failed_samples = [sample for sample, gtf in output_gtf_results if gtf is None]
    if failed_samples:
        raise RuntimeError(
            "StringTie assembly failed for sample(s): " + ", ".join(sorted(map(str, failed_samples)))
        )
    output_gtfs = [gtf for _sample, gtf in output_gtf_results]
    if not output_gtfs:
        raise RuntimeError("No StringTie GTF outputs were produced.")

    output_gtfs_txt = os.path.join(assembly_dir, "output_gtfs.txt")
    with open(output_gtfs_txt, 'w') as f:
        for gtf in output_gtfs:
            f.write(f"{os.path.abspath(gtf)}\n")
    
    log_message(
        "[INFO]",
        f"StringTie assembly completed: assembled GTFs={len(output_gtfs)}",
    )
    log_message(
        "[INFO]",
        f"StringTie GTF list: {output_gtfs_txt}",
        color="info",
        console=getattr(args, "debug", False),
    )

    log_message("[INFO]", "Merging assembled transcripts.", color="info")

    use_stringtie_merge = str(args.merge_method).lower() == "stringtie"
    reference_guided = str(args.assembly_mode).lower() == "reference-guided"

    if use_stringtie_merge:
        log_message("[INFO]", "Merge method: StringTie", color="info")
        merge_command = [
            'stringtie', '--merge', '-i',
            '-p', str(args.threads),
            '-o', os.path.join(tmp_dir, "stringtie_merged.gtf")
        ] + output_gtfs

        if args.optimized or _has_user_merge_filter(args):
            if reference_guided:
                min_fpkm, min_len, min_frac = 3.48, 1000, 0.311
            else:
                min_fpkm, min_len, min_frac = 7.19, 1000, 0.345
            merge_args = []
            if args.optimized or getattr(args, 'min_expr', None) is not None:
                min_fpkm = _override_or_default(args, 'min_expr', min_fpkm)
                merge_args.extend(['-F', str(min_fpkm)])
            if args.optimized or getattr(args, 'min_length', None) is not None:
                min_len = _override_or_default(args, 'min_length', min_len)
                merge_args.extend(['-m', str(min_len)])
            if args.optimized or getattr(args, 'min_frac', None) is not None:
                min_frac = _override_or_default(args, 'min_frac', min_frac)
                merge_args.extend(['-f', str(min_frac)])
            merge_command.extend(merge_args)
            log_message(
                "[INFO]",
                f"StringTie merge filter parameters: {' '.join(merge_args)}",
                color="info",
            )

        output_merge_gtf = os.path.join(tmp_dir, "stringtie_merged.gtf")
    else:
        log_message("[INFO]", "Merge method: TACO", color="info")
        taco_path = _resolve_taco_path(args)
        
        log_message("[INFO]", f"TACO tool found at {taco_path}", color="info", console=getattr(args, "debug", False))
        
        merge_command = [
            f'{taco_path}', '-p', str(args.threads)
        ]

        if args.optimized or _has_user_merge_filter(args):
            if reference_guided:
                min_expr, min_len, iso_frac = 11.7, 500, 0.380
            else:
                min_expr, min_len, iso_frac = 24.2, 500, 0.345
            merge_args = []
            if args.optimized or getattr(args, 'min_expr', None) is not None:
                min_expr = _override_or_default(args, 'min_expr', min_expr)
                merge_args.extend(['--filter-min-expr', str(min_expr)])
            if args.optimized or getattr(args, 'min_length', None) is not None:
                min_len = _override_or_default(args, 'min_length', min_len)
                merge_args.extend(["--filter-min-length", str(min_len)])
            if args.optimized or getattr(args, 'min_frac', None) is not None:
                iso_frac = _override_or_default(args, 'min_frac', iso_frac)
                merge_args.extend(["--isoform-frac", str(iso_frac)])
            merge_command.extend(merge_args)
            log_message(
                "[INFO]",
                f"TACO merge filter parameters: {' '.join(merge_args)}",
                color="info",
            )

        output_gtfs_txt =os.path.abspath(output_gtfs_txt)
        merge_command.extend([output_gtfs_txt])
        output_merge_gtf = os.path.join(tmp_dir, 'output/assembly.gtf')

    try:
        merge_section = "stringtie_merge" if use_stringtie_merge else "taco_merge"
        run_prep_command(merge_command, args, merge_section, cwd=tmp_dir)
        log_message("[INFO]", "Transcript merge completed.")
    except subprocess.CalledProcessError as e:
        stderr_text = _decode_stderr(e)
        log_message("[ERROR]", f"Error during GTF merge. Error: {stderr_text}", color="error")
        raise RuntimeError(f"GTF merge failed: {stderr_text}") from e
    finally:
        try:
            os.remove(output_gtfs_txt)
        except FileNotFoundError:
            pass

    log_message("[INFO]", "Comparing merged transcripts to reference annotation.", color="info")
    output_gffcompare_prefix = os.path.join(tmp_dir, "gffcompare")
    run_prep_command([
        'gffcompare', '-r', str(args.gene),
        '-o', str(output_gffcompare_prefix),
        str(output_merge_gtf)
    ], args, "gffcompare")
    log_message("[INFO]", "Reference comparison completed.")

    log_message("[INFO]", "Resolving transcript-gene assignments.", color="info")
    output_final_gtf = os.path.join(assembly_dir, CONSENSUS_GTF)
    output_assignment_tsv = os.path.join(assembly_dir, TRANSCRIPT_GENE_ASSIGNMENT_TSV)
    output_assignment_detail_tsv = (
        os.path.join(assembly_dir, TRANSCRIPT_GENE_ASSIGNMENT_DETAIL_TSV)
        if getattr(args, "detail", False)
        else None
    )
    input_tmap = sorted(glob.glob(os.path.join(tmp_dir, '**', 'gffcompare*.tmap'), recursive=True))
    input_gtf = sorted(glob.glob(os.path.join(tmp_dir, 'gffcompare.annotated.gtf')))
    if not input_gtf:
        input_gtf = sorted(glob.glob(os.path.join(tmp_dir, '**', 'gffcompare*.gtf'), recursive=True))
    if not input_tmap:
        generated = sorted(glob.glob(os.path.join(tmp_dir, '**', 'gffcompare*'), recursive=True))
        preview = ", ".join(generated[:10]) if generated else "none"
        raise FileNotFoundError(
            f"gffcompare tmap not found in {tmp_dir} (pattern: **/gffcompare*.tmap). "
            f"Generated gffcompare files: {preview}"
        )
    if not input_gtf:
        generated = sorted(glob.glob(os.path.join(tmp_dir, '**', 'gffcompare*'), recursive=True))
        preview = ", ".join(generated[:10]) if generated else "none"
        raise FileNotFoundError(
            f"gffcompare annotated GTF not found in {tmp_dir} (pattern: gffcompare.annotated.gtf or **/gffcompare*.gtf). "
            f"Generated gffcompare files: {preview}"
        )
    modify_gtf_with_mapping(
        input_tmap[0],
        input_gtf[0],
        output_final_gtf,
        reference_gtf=args.gene,
        assignment_tsv=output_assignment_tsv,
        assignment_detail_tsv=output_assignment_detail_tsv,
        gene_assignment=getattr(args, "gene_assignment", "strict"),
        include_overlap=getattr(args, "include_overlap", False),
        detail=getattr(args, "detail", False),
        debug=getattr(args, "debug", False),
    )
    log_message(
        "[SUCCESS]",
        f"Consensus transcript GTF: {os.path.abspath(output_final_gtf)}",
        color="success",
    )
    log_message(
        "[SUCCESS]",
        f"Transcript-gene assignment: {os.path.abspath(output_assignment_tsv)}",
        color="success",
    )
    if output_assignment_detail_tsv:
        log_message(
            "[SUCCESS]",
            f"Transcript-gene assignment detail: {os.path.abspath(output_assignment_detail_tsv)}",
            color="success",
        )

    if not args.debug:
        _cleanup_assembly_dir(assembly_dir)
        
    return output_final_gtf


def run_stringtie(bam, sample, stringtie_cmd, tmp_dir, args, stringtie_threads=None):
    output_gtf = os.path.join(tmp_dir, f"{sample}.gtf")
    if stringtie_threads is None:
        stringtie_threads = max(1, int(getattr(args, 'threads', 1)))
    
    stringtie_command = [
        *stringtie_cmd, '-u',
        '-p', str(stringtie_threads),
        '-j', str(getattr(args, 'junction', 2.0)),
        '-M', str(STRINGTIE_MULTIMAP_FRACTION),
        '-o', output_gtf
    ]

    if str(args.assembly_mode).lower() == "reference-guided" and hasattr(args, 'gene'):
        stringtie_command.extend(['-G', args.gene])
    
    stringtie_command.append(bam)

    try:
        run_prep_command(stringtie_command, args, f"stringtie-{sample}")
        return sample, output_gtf
    except subprocess.CalledProcessError as e:
        log_message("[ERROR]", f"Error during StringTie assembly for {sample}. Error: {_decode_stderr(e)}", color="error")
        return sample, None

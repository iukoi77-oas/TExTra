"""Extract junction reads and junction BED inputs for qual HITindex."""

import os
import subprocess
import tempfile


def _write_filtered_junction_bam(input_bam, flag_args, outfile, sam_threads):
    header = subprocess.check_output(["samtools", "view", "-H", input_bam])
    with tempfile.NamedTemporaryFile(
        mode="wb",
        suffix=".sam",
        prefix=os.path.basename(outfile) + ".",
        dir=os.path.dirname(outfile),
        delete=False,
    ) as sam_tmp:
        sam_tmp.write(header)
        sam_tmp_path = sam_tmp.name

    try:
        with open(sam_tmp_path, "ab") as sam_out:
            view_cmd = ["samtools", "view", "-@", str(sam_threads), *flag_args, input_bam]
            with subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
                assert proc.stdout is not None
                for line in proc.stdout:
                    fields = line.split("\t", 6)
                    if len(fields) > 5 and "N" in fields[5]:
                        sam_out.write(line.encode("utf-8"))
                stderr = proc.stderr.read() if proc.stderr is not None else ""
                ret = proc.wait()
                if ret != 0:
                    raise subprocess.CalledProcessError(ret, view_cmd, stderr=stderr)

        unsorted_bam = outfile + ".unsorted.bam"
        subprocess.run(["samtools", "view", "-bS", "-o", unsorted_bam, sam_tmp_path], check=True)
        subprocess.run(
            ["samtools", "sort", "-@", str(sam_threads), "-T", outfile, "-o", outfile, unsorted_bam],
            check=True,
        )
    finally:
        if os.path.exists(sam_tmp_path):
            os.remove(sam_tmp_path)
        unsorted_bam = outfile + ".unsorted.bam"
        if os.path.exists(unsorted_bam):
            os.remove(unsorted_bam)


def extract_junction(input_bam, sample, output_dir, args):
    """Extract splice-junction reads from BAM into strand/read-specific BAM files."""
    sam_threads = max(1, int(getattr(args, 'threads', 1)))

    if args.readtype == 'paired' and args.strand != 'none':
        # For paired-end data, process read1
        outfile1 = os.path.join(output_dir, f'{sample}_read1.bam')
        _write_filtered_junction_bam(input_bam, ["-f", "64", "-F", "256"], outfile1, sam_threads)
        subprocess.run(["samtools", "index", "-@", str(sam_threads), outfile1], check=True)

        # For paired-end data, process read2
        outfile2 = os.path.join(output_dir, f'{sample}_read2.bam')
        _write_filtered_junction_bam(input_bam, ["-f", "128", "-F", "256"], outfile2, sam_threads)
        subprocess.run(["samtools", "index", "-@", str(sam_threads), outfile2], check=True)
    elif args.readtype == 'single' or args.strand == 'none':
        # For single-end data
        outfile = os.path.join(output_dir, f'{sample}_read.bam')
        _write_filtered_junction_bam(input_bam, ["-F", "256"], outfile, sam_threads)
        subprocess.run(["samtools", "index", "-@", str(sam_threads), outfile], check=True)

    return f'{output_dir}/{sample}'

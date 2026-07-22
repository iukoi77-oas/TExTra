"""Build metaexon BED inputs and HITindex support records for qual mode."""

import collections
import os

import pybedtools

from util.qual.define_hitindex_constants import HITINDEX_MIN_JUNCTION_READS, HITINDEX_SPLIT_OVERLAP_BP
from util.qual.annotate_te_events import _build_exon_event_key


def _build_te_overlap_unique_exon_rows(transcript_exon_rows, te_overlap_keys):
    if not transcript_exon_rows or not te_overlap_keys:
        return [], {}
    agg = {}
    for row in transcript_exon_rows:
        chrom = str(row.get("exon_chrom", ""))
        strand = str(row.get("exon_strand", ""))
        try:
            start = int(row.get("exon_start", -1))
            end = int(row.get("exon_end", -1))
        except (TypeError, ValueError):
            continue
        if not chrom or end <= start:
            continue
        exon_coord = f"{chrom}:{start}-{end}"
        event_key = _build_exon_event_key(exon_coord, strand)
        if event_key not in te_overlap_keys:
            continue
        rec = agg.setdefault(
            event_key,
            {
                "exon_chrom": chrom,
                "exon_start": int(start),
                "exon_end": int(end),
                "exon_strand": strand,
                "gene_ids": set(),
                "transcript_ids": set(),
            },
        )
        gene_id = str(row.get("gene_id", "")).strip()
        transcript_id = str(row.get("transcript_id", "")).strip()
        if gene_id:
            rec["gene_ids"].add(gene_id)
        if transcript_id:
            rec["transcript_ids"].add(transcript_id)

    unique_rows = []
    unique_meta = {}
    for event_key, rec in agg.items():
        gene_id = ",".join(sorted(rec["gene_ids"]))
        transcript_id = ",".join(sorted(rec["transcript_ids"]))
        transcript_n = int(len(rec["transcript_ids"]))
        row = {
            "event_key": event_key,
            "exon_chrom": rec["exon_chrom"],
            "exon_start": int(rec["exon_start"]),
            "exon_end": int(rec["exon_end"]),
            "exon_strand": rec["exon_strand"],
            "gene_id": gene_id,
            "transcript_id": transcript_id,
            "transcript_n": transcript_n,
        }
        unique_rows.append(row)
        unique_meta[event_key] = {
            "gene_id": gene_id,
            "transcript_id": transcript_id,
            "transcript_n": transcript_n,
        }

    unique_rows = sorted(
        unique_rows,
        key=lambda r: (str(r["exon_chrom"]), int(r["exon_start"]), int(r["exon_end"]), str(r["exon_strand"])),
    )
    return unique_rows, unique_meta


def _build_unique_exon_bed(unique_exon_rows, out_bed_path):
    if not unique_exon_rows:
        return None
    os.makedirs(os.path.dirname(out_bed_path), exist_ok=True)
    with open(out_bed_path, "w", encoding="utf-8") as out:
        for row in unique_exon_rows:
            chrom = str(row["exon_chrom"])
            start = int(row["exon_start"])
            end = int(row["exon_end"])
            strand = str(row["exon_strand"])
            gene_id = str(row.get("gene_id", ""))
            exon_coord = f"{chrom}:{start}-{end}"
            # Keep DEXSeq-style info slots so existing getExons/calculateMetric logic can be reused directly.
            tx_n = max(1, int(row.get("transcript_n", 1)))
            bed_name = f"{exon_coord};{gene_id};TXPT:{tx_n};FE:0;internal:1;LE:0;singleexon:0"
            out.write(f"{chrom}\t{start}\t{end}\t{bed_name}\t0\t{strand}\n")
    return out_bed_path


def _copy_bed_with_gene_prefix(input_bed_path, out_bed_path, gene_prefix):
    """Copy BED6 records while prefixing the semicolon-delimited gene slot in name."""
    os.makedirs(os.path.dirname(out_bed_path), exist_ok=True)
    with open(input_bed_path, "r", encoding="utf-8") as in_fh, open(out_bed_path, "w", encoding="utf-8") as out_fh:
        for line_no, line in enumerate(in_fh, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise RuntimeError(f"Invalid BED at {input_bed_path}:{line_no}: expected at least 6 columns.")
            name_parts = fields[3].split(";")
            if len(name_parts) < 2:
                raise RuntimeError(
                    f"Invalid BED name at {input_bed_path}:{line_no}: expected exon;gene;... format."
                )
            name_parts[1] = f"{gene_prefix}{name_parts[1]}"
            fields[3] = ";".join(name_parts)
            out_fh.write("\t".join(fields) + "\n")
    return out_bed_path


def _build_union_bed(primary_bed_path, secondary_bed_path, out_bed_path):
    """Concatenate primary and internally-prefixed secondary BED records for one intersectBed pass."""
    os.makedirs(os.path.dirname(out_bed_path), exist_ok=True)
    with open(out_bed_path, "w", encoding="utf-8") as out_fh:
        for bed_path in [primary_bed_path, secondary_bed_path]:
            with open(bed_path, "r", encoding="utf-8") as in_fh:
                for line in in_fh:
                    if line.strip():
                        out_fh.write(line if line.endswith("\n") else line + "\n")
    return out_bed_path


def _flatten_unique_exon_support(
    HITdict,
    sample,
    overlap_bp=HITINDEX_SPLIT_OVERLAP_BP,
    read_threshold=HITINDEX_MIN_JUNCTION_READS,
    te_event_anno_map=None,
    unique_event_meta=None,
):
    if te_event_anno_map is None:
        te_event_anno_map = {}
    if unique_event_meta is None:
        unique_event_meta = {}
    rows = []
    for gene, exons in HITdict.items():
        for exon_key, metrics in exons.items():
            exon_coord = str(exon_key)
            chrom = ""
            start = -1
            end = -1
            if ":" in exon_coord and "-" in exon_coord:
                try:
                    chrom, span = exon_coord.split(":", 1)
                    start_s, end_s = span.split("-", 1)
                    start = int(start_s)
                    end = int(end_s)
                except (TypeError, ValueError):
                    chrom = ""
                    start = -1
                    end = -1
            nleft = int(metrics.get("nleft", 0))
            nright = int(metrics.get("nright", 0))
            total = int(nleft + nright)
            event_key = _build_exon_event_key(exon_coord, metrics.get("strand", ""))
            te_ann = te_event_anno_map.get(event_key, {})
            event_meta = unique_event_meta.get(event_key, {})
            rows.append(
                {
                    "sample": str(sample),
                    "gene_id": str(event_meta.get("gene_id", str(gene))),
                    "transcript_id": str(event_meta.get("transcript_id", "")),
                    "exon": str(exon_coord),
                    "chrom": str(chrom),
                    "start": int(start),
                    "end": int(end),
                    "strand": str(metrics.get("strand", "")),
                    "nleft": int(nleft),
                    "nright": int(nright),
                    "total_junction_reads": int(total),
                    "junction_supported": int(total >= int(read_threshold)),
                    "junction_support_rule": f"split_overlap>={int(overlap_bp)}bp,total_junction_reads>={int(read_threshold)}",
                    "te_overlap_label": str(te_ann.get("label", "no_overlap")),
                    "te_overlap_n": int(te_ann.get("te_overlap_n", 0)),
                    "te_overlap_bp_max": int(te_ann.get("te_overlap_bp_max", 0)),
                    "te_overlap_frac_max": float(te_ann.get("te_overlap_frac_max", 0.0)),
                    "te_boundary_hit_any": int(te_ann.get("te_boundary_hit_any", 0)),
                    "te_overlap_pass_any": int(te_ann.get("te_overlap_pass_any", 0)),
                    "te_splice_site_repeat_TE": str(te_ann.get("te_splice_site_repeat_TE", "")),
                    "te_other_overlap_TE": str(te_ann.get("te_other_overlap_TE", "")),
                }
            )
    return rows


def parse_bed(bedfile):
    """Parse exon BED records into a nested gene->transcript structure."""
    genedict = collections.defaultdict(lambda: collections.defaultdict(list))

    with open(bedfile, 'r', encoding="utf-8") as f:
        for line_no, line in enumerate(f, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 6:
                raise ValueError(f"Invalid exon BED at {bedfile}:{line_no}: expected at least 6 columns.")
            chrom, start, end, name, score, strand = fields[:6]
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError as exc:
                raise ValueError(f"Invalid exon BED at {bedfile}:{line_no}: start/end must be integers.") from exc
            if end_i <= start_i:
                raise ValueError(f"Invalid exon BED at {bedfile}:{line_no}: end must be greater than start.")
            if strand not in {"+", "-"}:
                raise ValueError(f"Invalid exon BED at {bedfile}:{line_no}: strand must be '+' or '-'.")
            parts = name.split(":")
            if len(parts) != 3 or any(not part for part in parts):
                raise ValueError(
                    f"Invalid exon BED at {bedfile}:{line_no}: name must use gene:transcript:exon format."
                )
            gene, transcript, exon = parts[0], parts[1], parts[2]
            genedict[gene][transcript].append((chrom, start_i, end_i, strand, exon))

    if not genedict:
        raise ValueError(f"No valid exon records were found in BED file: {bedfile}")

    return genedict


def metaexon_bed(input_bed, annotation_dir, args, tmp_dir=None):
    """Build merged/buffered metaexon BED files from transcript exon BED."""
    if tmp_dir is None:
        tmp_dir = annotation_dir
    os.makedirs(annotation_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    outfh_path = os.path.join(annotation_dir, 'metaexon.bed')
    outfhbuffer_path = os.path.join(tmp_dir, f'buffer_ss3-{args.ss3buffer}_ss5-{args.ss5buffer}_metaexon.bed')
    outfhconst_path = os.path.join(tmp_dir, 'constituent_metaexon.bed')

    with open(outfh_path, 'w') as outfh, open(outfhbuffer_path, 'w') as outfhbuffer, open(outfhconst_path, 'w') as outfhconst:
        genedict = parse_bed(input_bed)
        exontypes = ['FE', 'internal', 'LE', 'singleexon']

        for gene, transcripts in genedict.items():
            exons = []
            exons_const = []
            ntxpt = 'TXPT:' + str(len(genedict[gene]))

            for transcript, transcript_exons in transcripts.items():
                exons_sorted = sorted(transcript_exons, key=lambda x: x[1])
                strand = exons_sorted[0][3]

                if len(exons_sorted) == 1:
                    exon_types = ['singleexon']
                else:
                    if strand == '+':
                        exon_types = ['FE'] + ['internal'] * (len(exons_sorted) - 2) + ['LE']
                    elif strand == '-':
                        exon_types = ['LE'] + ['internal'] * (len(exons_sorted) - 2) + ['FE']
                    else:
                        continue

                for i, exon in enumerate(exons_sorted):
                    chrom, start, end, strand, exon_id = exon
                    exon_type = exon_types[i]
                    exons.append(f"{chrom} {start} {end} {exon_type} . {strand}")
                    exons_const.append(f"{chrom} {start} {end} {start}-{end} . {strand}")

            if not exons:
                continue
            bedscratch = pybedtools.BedTool('\n'.join(exons), from_string=True)
            bedsortmerge = bedscratch.sort().merge(s=True, c="4,6", o="collapse,distinct")
            bedlist = str(bedsortmerge).split('\n')[:-1]

            bedscratch_const = pybedtools.BedTool('\n'.join(exons_const), from_string=True)
            bedsortmerge_const = bedscratch_const.sort().merge(s=True, c="4,6", o="collapse,distinct")
            bedlist_const = str(bedsortmerge_const).split('\n')[:-1]

            for ex in bedlist:
                exhere = ex.split('\t')
                exherename = f"{exhere[0]}:{exhere[1]}-{exhere[2]}"
                extypes = exhere[3].split(',')
                meta_strand = exhere[4].split(",")[0]
                ntypeset = [f"{x}:{extypes.count(x)}" for x in exontypes]

                # Apply splice-site buffers in transcript-strand orientation.
                if meta_strand == '+':
                    exstart = max(0, int(exhere[1]) - args.ss3buffer)
                    exend = int(exhere[2]) + args.ss5buffer
                elif meta_strand == '-':
                    exstart = max(0, int(exhere[1]) - args.ss5buffer)
                    exend = int(exhere[2]) + args.ss3buffer
                else:
                    continue

                outfh.write(f"{exhere[0]}\t{exhere[1]}\t{exhere[2]}\t{exherename};{gene};{ntxpt};{';'.join(ntypeset)}\t0\t{meta_strand}\n")
                outfhbuffer.write(f"{exhere[0]}\t{exstart}\t{exend}\t{exherename};{gene};{ntxpt};{';'.join(ntypeset)}\t0\t{meta_strand}\n")

            for ex in bedlist_const:
                exhere = ex.split('\t')
                exherename = f"{exhere[0]}:{exhere[1]}-{exhere[2]}"
                const_exon_type = exhere[3]
                const_strand = exhere[4].split(",")[0]
                if const_strand not in {"+", "-"}:
                    continue
                outfhconst.write(f"{exhere[0]}\t{exhere[1]}\t{exhere[2]}\t{exherename};{gene};{ntxpt};{const_exon_type}\t0\t{const_strand}\n")

    return outfhbuffer_path

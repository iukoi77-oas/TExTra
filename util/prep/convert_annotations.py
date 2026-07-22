"""Convert prep gene and TE annotations to BED outputs."""

import os
from collections import defaultdict

import pybedtools

from util.common.write_logs import log_message

def gtf_to_bed(gtf_path, bed_output_path, annotation_type):
    """Convert gene or TE annotation GTF/GFF records into normalized BED6."""
    if annotation_type == 'gene':
        _gene_gtf_to_bed(gtf_path, bed_output_path)
        return

    gtf = pybedtools.BedTool(gtf_path)
    file_extension = os.path.splitext(gtf_path)[1].lower()

    # Validate annotation_type
    if annotation_type not in {'TE', 'gene'}:
        raise ValueError("Unsupported annotation type. Use 'TE' or 'gene'.")

    filtered_simple = {'count': 0}
    dropped_missing_cf = {'count': 0}
    skip_classes = {"simple_repeat", "low_complexity", "satellite"}

    def _format_to_bed(feature):
        attributes = {}
        gene_id = None
        family_id = None
        class_id = None
        transcript_id = None
        exon_number = None

        try:
            attributes = feature.attrs
            gene_id = attributes.get("gene_id", "Unknown")
            family_id = attributes.get("family_id", attributes.get("family", attributes.get("Family", None)))
            class_id = attributes.get("class_id", attributes.get("class", attributes.get("Class", None)))
            transcript_id = attributes.get("transcript_id", None)
            exon_number = attributes.get("exon_number", None)
        except Exception:
            attributes = {}
            if not gene_id:
                gene_id = "Unknown"
        
        formatted_id = None

        # Format the ID based on annotation type
        if annotation_type == 'TE':  # TE annotation
            if file_extension == ".gtf":
                if not family_id or not class_id:
                    dropped_missing_cf['count'] += 1
                    return None
                if str(class_id).strip().lower() in skip_classes:
                    filtered_simple['count'] += 1
                    return None
                formatted = f"{gene_id}:{family_id}:{class_id}"
            elif file_extension == ".gff":
                ID_id = attributes.get("ID", None)
                classification_id = attributes.get("Classification", None)
                if classification_id:
                    parts = classification_id.split("/")
                    class_id = parts[0]
                    family_id = parts[1] if len(parts) > 1 else None
                if not family_id or not class_id:
                    dropped_missing_cf['count'] += 1
                    return None
                if str(class_id).strip().lower() in skip_classes:
                    filtered_simple['count'] += 1
                    return None
                formatted = f"{ID_id if ID_id else 'Unknown'}:{family_id}:{class_id}"
            formatted_id = f"{feature.chrom}|{feature.start}|{feature.end}|{formatted}|{feature.score}|{feature.strand}"
        if formatted_id is None:
            return None

        # Return BED format: seqname, start, end, formatted_id, score, strand
        return pybedtools.create_interval_from_list([
            feature.chrom,
            str(feature.start), 
            str(feature.end),
            formatted_id,
            feature.score if feature.score else ".",
            feature.strand
        ])

    # Apply the custom formatting function
    bed = gtf.each(_format_to_bed)
    sorted_bed = bed.sort()
    sorted_bed.saveas(bed_output_path)

    if annotation_type == 'TE':
        if dropped_missing_cf['count'] > 0:
            log_message("[WARNING]", f"Dropped {dropped_missing_cf['count']} TE records with missing Class/Family.", color="warning")
        if filtered_simple['count'] > 0:
            log_message("[WARNING]", f"Filtered {filtered_simple['count']} TE records from Simple_repeat/Low_complexity/Satellite.", color="warning")

    if os.path.exists(bed_output_path) and os.path.getsize(bed_output_path) == 0:
        log_message("[ERROR]", f"Please check that the reference gene/TE annotation files.", color="error")
        raise RuntimeError(f"Converted BED is empty: {bed_output_path}")


def _gene_gtf_to_bed(gtf_path, bed_output_path):
    """Convert transcript exon features to BED6 with strand-aware exon numbering fallback."""
    exons_by_tx = defaultdict(list)
    with open(gtf_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            attrs, _order = _parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id", "")
            transcript_id = attrs.get("transcript_id", "")
            if not gene_id or not transcript_id:
                continue
            try:
                start0 = int(fields[3]) - 1
                end0 = int(fields[4])
            except ValueError:
                continue
            if end0 <= start0:
                continue
            exons_by_tx[transcript_id].append(
                {
                    "chrom": fields[0],
                    "start0": start0,
                    "end0": end0,
                    "score": fields[5] if fields[5] else ".",
                    "strand": fields[6],
                    "gene_id": gene_id,
                    "transcript_id": transcript_id,
                    "exon_number": attrs.get("exon_number", ""),
                }
            )

    bed_rows = []
    for transcript_id, exons in exons_by_tx.items():
        strand = exons[0]["strand"] if exons else "."
        reverse = strand == "-"
        ordered = sorted(exons, key=lambda row: (row["start0"], row["end0"]), reverse=reverse)
        fallback_numbers = {id(row): idx for idx, row in enumerate(ordered, start=1)}
        for row in exons:
            exon_number = row["exon_number"] or str(fallback_numbers[id(row)])
            formatted_id = f"{row['gene_id']}:{transcript_id}:exon_{exon_number}"
            bed_rows.append(
                [
                    row["chrom"],
                    row["start0"],
                    row["end0"],
                    formatted_id,
                    row["score"],
                    row["strand"],
                ]
            )

    bed_rows.sort(key=lambda row: (row[0], int(row[1]), int(row[2]), row[5], row[3]))
    with open(bed_output_path, "w", encoding="utf-8") as out:
        for row in bed_rows:
            out.write("\t".join(map(str, row)) + "\n")

    if os.path.exists(bed_output_path) and os.path.getsize(bed_output_path) == 0:
        log_message("[ERROR]", f"Please check that the reference gene annotation file: {gtf_path}", color="error")
        raise RuntimeError(f"Converted gene BED is empty: {bed_output_path}")

def merge_and_sort_bed(bed_file1, bed_file2, output_bed_file):
    """Merge optional BED inputs, sort them, and write unique BED6 records."""
    bed1 = pybedtools.BedTool(bed_file1)
    if bed_file2 == None:
        merged_bed = bed1.sort()
    else:
        bed2 = pybedtools.BedTool(bed_file2)
        # Concatenate and sort the BED files
        merged_bed = bed1.cat(bed2, postmerge=False).sort()
    
    # Save the sorted BED file
    merged_bed = merged_bed.sort()
    with open(output_bed_file, "w", encoding="utf-8") as out:
        seen = set()
        for feature in merged_bed:
            row = "\t".join(feature.fields[:6])
            if row in seen:
                continue
            seen.add(row)
            out.write(row + "\n")

def _parse_gtf_attributes(attr_field):
    attrs = {}
    order = []
    for token in str(attr_field).split(";"):
        token = token.strip()
        if not token or " " not in token:
            continue
        key, value = token.split(" ", 1)
        key = key.strip()
        attrs[key] = value.strip().strip('"')
        if key not in order:
            order.append(key)
    return attrs, order


def _format_gtf_attributes(attrs, order):
    keys = list(order)
    for key in attrs:
        if key not in keys:
            keys.append(key)
    return " ".join(f'{key} "{attrs[key]}";' for key in keys if attrs.get(key) not in {None, ""})

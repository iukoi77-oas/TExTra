"""Map assembled transcripts to reference genes."""

import csv
from collections import Counter, defaultdict

from util.common.write_logs import log_message


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


def _load_reference_annotation(gtf_path):
    gene_to_strands = defaultdict(set)
    gene_to_name = {}
    transcript_to_gene = {}
    transcript_to_strand = {}
    transcript_to_name = {}
    with open(gtf_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            strand = fields[6].strip()
            if strand not in {"+", "-"}:
                continue
            attrs, _order = _parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id", "").strip()
            if not gene_id:
                continue
            gene_to_strands[gene_id].add(strand)
            gene_name = attrs.get("gene_name") or attrs.get("Name") or attrs.get("gene")
            if gene_name and gene_id not in gene_to_name:
                gene_to_name[gene_id] = gene_name
            transcript_id = attrs.get("transcript_id", "").strip()
            if transcript_id:
                transcript_to_gene[transcript_id] = gene_id
                transcript_to_strand[transcript_id] = strand
                if gene_name:
                    transcript_to_name[transcript_id] = gene_name
    return {
        "gene_to_strands": gene_to_strands,
        "gene_to_name": gene_to_name,
        "transcript_to_gene": transcript_to_gene,
        "transcript_to_strand": transcript_to_strand,
        "transcript_to_name": transcript_to_name,
    }


def _load_query_transcript_info(gtf_path):
    info = {}
    exon_numbers = defaultdict(set)
    exon_rows = Counter()
    order = []
    with open(gtf_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs, _attr_order = _parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id", "")
            if not transcript_id:
                continue
            if transcript_id not in info:
                order.append(transcript_id)
                info[transcript_id] = {
                    "transcript_id": transcript_id,
                    "original_gene_id": attrs.get("gene_id", ""),
                    "query_strand": fields[6],
                    "query_chrom": fields[0],
                    "query_start": fields[3],
                    "query_end": fields[4],
                    "class_code": attrs.get("class_code", ""),
                    "reference_transcript_id": attrs.get("cmp_ref", ""),
                    "reference_gene_name": attrs.get("cmp_ref_gene") or attrs.get("gene_name", ""),
                }
            if fields[2] == "transcript":
                info[transcript_id].update(
                    {
                        "original_gene_id": attrs.get("gene_id", info[transcript_id]["original_gene_id"]),
                        "query_strand": fields[6],
                        "query_chrom": fields[0],
                        "query_start": fields[3],
                        "query_end": fields[4],
                        "class_code": attrs.get("class_code", info[transcript_id]["class_code"]),
                        "reference_transcript_id": attrs.get("cmp_ref", info[transcript_id]["reference_transcript_id"]),
                        "reference_gene_name": attrs.get("cmp_ref_gene") or attrs.get("gene_name", info[transcript_id]["reference_gene_name"]),
                    }
                )
            elif fields[2] == "exon":
                exon_rows[transcript_id] += 1
                exon_number = attrs.get("exon_number")
                if exon_number:
                    exon_numbers[transcript_id].add(str(exon_number))
    for transcript_id, row in info.items():
        row["exon_count"] = len(exon_numbers[transcript_id]) if exon_numbers[transcript_id] else int(exon_rows[transcript_id])
    return info, order


def _read_tmap_candidates(input_tmap):
    candidates = defaultdict(list)
    with open(input_tmap, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            if fields[0] == "ref_gene_id" or fields[4] == "qry_id":
                continue
            ref_gene_id, ref_tx_id, class_code, qry_gene_id, qry_tx_id = fields[:5]
            candidates[qry_tx_id].append(
                {
                    "reference_gene_id": ref_gene_id,
                    "reference_transcript_id": ref_tx_id,
                    "class_code": class_code,
                    "query_gene_id": qry_gene_id,
                    "query_transcript_id": qry_tx_id,
                }
            )
    return candidates


def _pick_tmap_candidate(candidates, row=None, reference=None):
    priority = {
        "=": 0,
        "c": 1,
        "j": 2,
        "k": 3,
        "m": 4,
        "e": 5,
        "o": 6,
        "i": 7,
        "n": 8,
        "y": 9,
        "s": 10,
        "x": 11,
        "p": 12,
        "u": 13,
    }

    query_strand = (row or {}).get("query_strand", "")

    def _strand_rank(candidate):
        if not reference or query_strand not in {"+", "-"}:
            return 1
        ref_gene_id = candidate.get("reference_gene_id", "")
        if ref_gene_id in {"", "-"}:
            return 1
        reference_strand = _infer_reference_strand(candidate, reference)
        if reference_strand == query_strand:
            return 0
        if reference_strand in {"+", "-"}:
            return 2
        return 1

    return sorted(
        candidates,
        key=lambda c: (
            _strand_rank(c),
            priority.get(c.get("class_code", ""), 50),
            c.get("reference_gene_id", "-") == "-",
            c.get("reference_gene_id", ""),
            c.get("reference_transcript_id", ""),
        ),
    )[0]


def _infer_reference_strand(candidate, reference):
    ref_tx = candidate.get("reference_transcript_id", "")
    ref_gene = candidate.get("reference_gene_id", "")
    tx_strand = reference["transcript_to_strand"].get(ref_tx)
    if tx_strand:
        return tx_strand
    gene_strands = reference["gene_to_strands"].get(ref_gene, set())
    if len(gene_strands) == 1:
        return next(iter(gene_strands))
    if gene_strands:
        return ",".join(sorted(gene_strands))
    return ""


def _infer_unassigned_gene_type(assigned_gene_id, assignment_status, association_type, class_code, reference_gene_id):
    if assigned_gene_id:
        return ""
    if association_type == "antisense" or class_code in {"s", "x"} or assignment_status == "unassigned_strand_mismatch":
        return "antisense"
    if association_type == "intronic" or class_code == "i":
        return "intronic"
    if assignment_status == "intergenic_or_unmapped" or class_code == "u" or not reference_gene_id:
        return "intergenic"
    return "ambiguous"


def _build_assignment(row, candidate, reference, gene_assignment, include_overlap):
    confident_codes = {"=", "c", "j", "k"}
    flagged_codes = {"m"}
    flag_only_codes = {"e"}
    novel_locus_codes = {"i", "o", "n", "y"}
    excluded_codes = {"s", "x", "p"}

    original_gene_id = row.get("original_gene_id", "")
    final_gene_id = original_gene_id
    assigned_gene_id = ""
    associated_gene_id = ""
    antisense_gene_id = ""
    host_gene_analysis = "no"
    association_type = ""
    reason = ""

    ref_gene_id = candidate.get("reference_gene_id", "") if candidate else ""
    ref_tx_id = candidate.get("reference_transcript_id", "") if candidate else ""
    class_code = candidate.get("class_code", row.get("class_code", "")) if candidate else row.get("class_code", "")
    if not row.get("class_code") and class_code:
        row["class_code"] = class_code
    reference_strand = _infer_reference_strand(candidate, reference) if candidate and ref_gene_id != "-" else ""
    query_strand = row.get("query_strand", "")
    strand_consistent = "unknown"
    if query_strand in {"+", "-"} and reference_strand in {"+", "-"}:
        strand_consistent = "yes" if query_strand == reference_strand else "no"
    elif reference_strand:
        strand_consistent = "ambiguous"

    if not candidate or ref_gene_id in {"", "-"}:
        assignment_status = "intergenic_or_unmapped"
        reason = "no_reference_gene"
    elif class_code in confident_codes:
        if strand_consistent == "no":
            assignment_status = "unassigned_strand_mismatch"
            associated_gene_id = ref_gene_id
            reason = "confident_class_code_but_strand_mismatch"
        else:
            assignment_status = "assigned_confident"
            final_gene_id = ref_gene_id
            assigned_gene_id = ref_gene_id
            reason = "confident_class_code"
    elif class_code in flagged_codes:
        if strand_consistent == "no":
            assignment_status = "unassigned_strand_mismatch"
            associated_gene_id = ref_gene_id
            reason = "flagged_class_code_but_strand_mismatch"
        else:
            assignment_status = "assigned_flagged"
            final_gene_id = ref_gene_id
            assigned_gene_id = ref_gene_id
            reason = "flagged_class_code"
    elif class_code in flag_only_codes:
        assignment_status = "flagged_candidate"
        associated_gene_id = ref_gene_id
        association_type = "exonic_overlap"
        reason = "flag_only_class_code"
    elif class_code in novel_locus_codes:
        assignment_status = "novel_locus"
        associated_gene_id = ref_gene_id
        if class_code == "i":
            association_type = "intronic"
        elif class_code == "o":
            association_type = "exonic_overlap"
        elif class_code == "n":
            association_type = "retained_intron"
        elif class_code == "y":
            association_type = "contains_reference"
        if include_overlap and strand_consistent != "no":
            host_gene_analysis = "yes"
            reason = "novel_locus_allowed_for_host_gene_analysis"
        elif strand_consistent == "no":
            reason = "novel_locus_with_strand_mismatch"
        else:
            reason = "novel_locus_not_assigned_in_strict_mode"
    elif class_code in excluded_codes:
        assignment_status = "excluded"
        antisense_gene_id = ref_gene_id
        association_type = "antisense" if class_code in {"s", "x"} else "excluded_overlap"
        reason = f"excluded_class_code_{class_code}"
    elif class_code == "u":
        assignment_status = "intergenic_or_unmapped"
        reason = "intergenic_class_code"
    else:
        assignment_status = "unassigned_class_code"
        associated_gene_id = ref_gene_id
        reason = f"class_code_{class_code or 'missing'}_not_assigned"

    if gene_assignment != "strict":
        raise ValueError(f"Unsupported gene_assignment policy: {gene_assignment}")

    ref_gene_name = reference["gene_to_name"].get(ref_gene_id, "")
    if not ref_gene_name and ref_tx_id:
        ref_gene_name = reference["transcript_to_name"].get(ref_tx_id, "")
    if not ref_gene_name:
        ref_gene_name = row.get("reference_gene_name", "")
    reference_gene_id = "" if ref_gene_id == "-" else ref_gene_id
    unassigned_gene_type = _infer_unassigned_gene_type(
        assigned_gene_id=assigned_gene_id,
        assignment_status=assignment_status,
        association_type=association_type,
        class_code=class_code,
        reference_gene_id=reference_gene_id,
    )

    return {
        "transcript_id": row.get("transcript_id", ""),
        "original_gene_id": original_gene_id,
        "final_gene_id": final_gene_id,
        "assigned_gene_id": assigned_gene_id,
        "reference_gene_id": reference_gene_id,
        "reference_gene_name": ref_gene_name,
        "reference_transcript_id": "" if ref_tx_id == "-" else ref_tx_id,
        "associated_gene_id": associated_gene_id,
        "antisense_gene_id": antisense_gene_id,
        "class_code": class_code,
        "query_chrom": row.get("query_chrom", ""),
        "query_start": row.get("query_start", ""),
        "query_end": row.get("query_end", ""),
        "query_strand": query_strand,
        "reference_strand": reference_strand,
        "strand_consistent": strand_consistent,
        "exon_count": row.get("exon_count", 0),
        "assignment_status": assignment_status,
        "unassigned_gene_type": unassigned_gene_type,
        "association_type": association_type,
        "host_gene_analysis": host_gene_analysis,
        "reason": reason,
    }


def _assignment_public_row(row):
    explicit_reference_assignment = bool(row.get("assigned_gene_id"))
    return {
        "transcript_id": row.get("transcript_id", ""),
        "original_gene_id": row.get("original_gene_id", ""),
        "assigned_gene_id": row.get("final_gene_id", row.get("original_gene_id", "")),
        "assigned_gene_name": row.get("reference_gene_name", "") if explicit_reference_assignment else "",
        "unassigned_gene_type": row.get("unassigned_gene_type", ""),
    }


def _write_assignment_table(assignments, output_path):
    columns = [
        "transcript_id",
        "original_gene_id",
        "assigned_gene_id",
        "assigned_gene_name",
        "unassigned_gene_type",
    ]
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in assignments:
            writer.writerow(_assignment_public_row(row))


def _write_assignment_detail_table(assignments, output_path):
    columns = [
        "transcript_id",
        "original_gene_id",
        "assigned_gene_id",
        "assigned_gene_name",
        "reference_gene_id",
        "reference_gene_name",
        "reference_transcript_id",
        "class_code",
        "strand_consistent",
        "exon_count",
        "unassigned_gene_type",
    ]
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in assignments:
            out_row = dict(row)
            out_row.update(_assignment_public_row(row))
            writer.writerow(out_row)


def _annotate_attributes(attrs, order, assignment):
    original_gene_id = assignment.get("original_gene_id", "")
    final_gene_id = assignment.get("final_gene_id", original_gene_id)
    attrs["gene_id"] = final_gene_id
    attrs["original_gene_id"] = original_gene_id
    attrs["final_gene_id"] = final_gene_id
    attrs["gene_assignment_status"] = assignment.get("assignment_status", "")
    attrs["strand_consistent"] = assignment.get("strand_consistent", "")
    if assignment.get("reference_strand"):
        attrs["reference_strand"] = assignment["reference_strand"]
    if assignment.get("assigned_gene_id"):
        attrs["assigned_gene_id"] = assignment["assigned_gene_id"]
    if assignment.get("reference_gene_id"):
        attrs["reference_gene_id"] = assignment["reference_gene_id"]
    if assignment.get("associated_gene_id"):
        attrs["associated_gene_id"] = assignment["associated_gene_id"]
    if assignment.get("antisense_gene_id"):
        attrs["antisense_gene_id"] = assignment["antisense_gene_id"]
    if assignment.get("host_gene_analysis"):
        attrs["host_gene_analysis"] = assignment["host_gene_analysis"]
    if assignment.get("reason"):
        attrs["gene_assignment_reason"] = assignment["reason"]
    return _format_gtf_attributes(attrs, order)


def _format_count(value):
    return f"{int(value):,}"


def _format_percent(value, denominator):
    if not denominator:
        return "0.0%"
    return f"{(float(value) / float(denominator) * 100.0):.1f}%"


def _format_count_with_percent(value, denominator):
    return f"{_format_count(value)} ({_format_percent(value, denominator)})"


def modify_gtf_with_mapping(
    input_tmap,
    input_gtf,
    output_gtf,
    reference_gtf=None,
    assignment_tsv=None,
    assignment_detail_tsv=None,
    gene_assignment="strict",
    include_overlap=False,
    detail=False,
    debug=False,
):
    """Rewrite assembled transcript gene IDs using gffcompare mapping evidence."""
    if gene_assignment != "strict":
        raise ValueError(f"Unsupported gene_assignment policy: {gene_assignment}")

    query_info, transcript_order = _load_query_transcript_info(input_gtf)
    tmap_candidates = _read_tmap_candidates(input_tmap)
    reference = _load_reference_annotation(reference_gtf) if reference_gtf else {
        "gene_to_strands": defaultdict(set),
        "gene_to_name": {},
        "transcript_to_gene": {},
        "transcript_to_strand": {},
        "transcript_to_name": {},
    }

    assignment_by_tx = {}
    assignments = []
    for transcript_id in transcript_order:
        row = query_info[transcript_id]
        candidate = None
        if transcript_id in tmap_candidates:
            candidate = _pick_tmap_candidate(tmap_candidates[transcript_id], row=row, reference=reference)
        assignment = _build_assignment(row, candidate, reference, gene_assignment, include_overlap)
        assignment_by_tx[transcript_id] = assignment
        assignments.append(assignment)

    if assignment_tsv:
        _write_assignment_table(assignments, assignment_tsv)
    if assignment_detail_tsv:
        _write_assignment_detail_table(assignments, assignment_detail_tsv)

    input_feature_rows = 0
    output_feature_rows = 0
    replaced_gene_rows = 0
    with open(input_gtf, 'r') as in_file, open(output_gtf, 'w') as out_file:
        for line in in_file:
            if line.startswith("#"):
                out_file.write(line)
                continue
            
            if not line.strip():
                continue
            input_feature_rows += 1
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            attrs, order = _parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id", "")
            assignment = assignment_by_tx.get(transcript_id)
            if assignment:
                if assignment.get("assignment_status") == "excluded":
                    continue
                original_gene_id = attrs.get("gene_id", "")
                fields[8] = _annotate_attributes(attrs, order, assignment)
                if assignment.get("final_gene_id", original_gene_id) != original_gene_id:
                    replaced_gene_rows += 1
            
            # Write the modified line to the output file
            out_file.write('\t'.join(fields) + '\n')
            output_feature_rows += 1

    status_counts = Counter(row["assignment_status"] for row in assignments)
    class_counts = Counter(
        row.get("class_code", "missing") or "missing"
        for row in assignments
        if row.get("assignment_status") == "excluded"
    )
    excluded = int(status_counts.get("excluded", 0))
    kept = max(0, len(assignments) - excluded)
    reference_associated = (
        int(status_counts.get("assigned_confident", 0))
        + int(status_counts.get("assigned_flagged", 0))
        + int(status_counts.get("flagged_candidate", 0))
    )
    assigned_reference = int(status_counts.get("assigned_confident", 0)) + int(status_counts.get("assigned_flagged", 0))
    reference_overlap_not_reassigned = int(status_counts.get("flagged_candidate", 0))
    novel_locus = int(status_counts.get("novel_locus", 0))
    intergenic_or_unmapped = int(status_counts.get("intergenic_or_unmapped", 0))
    retained_other = max(0, kept - reference_associated - novel_locus - intergenic_or_unmapped)

    log_message(
        "[INFO]",
        "\n".join(
            [
                "Transcript assignment summary:",
                f"  Total transcripts: {_format_count(len(assignments))}",
                f"  Retained:          {_format_count_with_percent(kept, len(assignments))}",
                f"  Excluded:          {_format_count_with_percent(excluded, len(assignments))}",
            ]
        ),
    )
    if detail:
        retained_lines = [
            f"  Reference-associated: {_format_count_with_percent(reference_associated, kept)}",
            f"    Assigned to reference:            {_format_count(assigned_reference)}",
            f"    Reference overlap not reassigned: {_format_count(reference_overlap_not_reassigned)}",
            f"  Novel locus:          {_format_count_with_percent(novel_locus, kept)}",
            f"  Intergenic/unmapped:  {_format_count_with_percent(intergenic_or_unmapped, kept)}",
        ]
        if retained_other:
            retained_lines.append(f"  Other retained:       {_format_count_with_percent(retained_other, kept)}")

        excluded_lines = []
        for code in ["s", "x", "p"]:
            count = int(class_counts.get(code, 0))
            if count:
                excluded_lines.append(f"  class_code {code}:       {_format_count_with_percent(count, excluded)}")
        excluded_other = max(0, excluded - sum(int(class_counts.get(code, 0)) for code in ["s", "x", "p"]))
        if excluded_other:
            excluded_lines.append(f"  Other excluded:      {_format_count_with_percent(excluded_other, excluded)}")
        if not excluded_lines:
            excluded_lines.append("  none")

        detail_summary = "\n".join(
            [
                "Transcript classification details:",
                "",
                "Retained (percent of retained transcripts):",
                *retained_lines,
                "",
                "Excluded (percent of excluded transcripts):",
                *excluded_lines,
            ]
        )
        log_message("[INFO]", detail_summary)
    log_message(
        "[DEBUG]",
        (
            f"GTF assignment summary: transcripts={len(assignments)}, status_counts={dict(status_counts)}, "
            f"input_features={input_feature_rows}, output_features={output_feature_rows}, "
            f"replaced_gene_rows={replaced_gene_rows}, assignment_table={assignment_tsv or 'not_written'}"
        ),
        console=False,
    )
    if output_feature_rows == 0:
        message = (
            "Mapped consensus GTF is empty after gffcompare mapping. "
            f"input_tmap={input_tmap}, input_gtf={input_gtf}, output_gtf={output_gtf}"
        )
        log_message("[ERROR]", message, color="error")
        raise RuntimeError(message)
    elif input_feature_rows > 0 and output_feature_rows <= max(1, int(input_feature_rows * 0.2)):
        log_message(
            "[WARNING]",
            (
                "Mapped consensus GTF rows are much lower than input. "
                f"input_features={input_feature_rows}, output_features={output_feature_rows}"
            ),
            color="warning",
        )

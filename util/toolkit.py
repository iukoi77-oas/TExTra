import pybedtools
import re
import pandas as pd
import csv
import sys
import os
import subprocess
from collections import defaultdict, Counter
from tempfile import NamedTemporaryFile
from util.logger import *

# Filter the GTF file to keep only exons
def filter_exon_from_gtf(gtf_file, output_exon_gtf):
    gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)
    exon_df = gtf_df[gtf_df[2] == 'exon']
    exon_df.to_csv(output_exon_gtf, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)

# Convert GTF to BED
def gtf_to_bed(gtf_path, bed_output_path, annotation_type):
    gtf = pybedtools.BedTool(gtf_path)
    file_extension = os.path.splitext(gtf_path)[1].lower()

    # Validate annotation_type
    if annotation_type not in {'TE', 'gene'}:
        raise ValueError("Unsupported annotation type. Use 'TE' or 'gene'.")

    def format_to_bed(feature):
        attributes_column = feature[8]  # Extract the attributes column
        attributes = None
        gene_id = None
        family_id = None
        class_id = None
        transcript_id = None
        exon_number = None

        try:
            attributes = feature.attrs
            gene_id = attributes.get("gene_id", None)
            family_id = attributes.get("family_id", None)
            class_id = attributes.get("class_id", None)
            transcript_id = attributes.get("transcript_id", None)
            exon_number = attributes.get("exon_number", None)
        except Exception as e:
            pass
        
        formatted_id = None

        # Format the ID based on annotation type
        if annotation_type == 'TE':  # TE annotation
            if file_extension == ".gtf":
                if gene_id and family_id and class_id:
                    formatted = f"{gene_id}:{family_id}:{class_id}"
            elif file_extension == ".gff":
                ID_id = attributes.get("ID", None)
                classification_id = attributes.get("Classification", None)
                if classification_id:
                    parts = classification_id.split("/")
                    class_id = parts[0]
                    family_id = parts[1] if len(parts) > 1 else "Unknown"                
                if ID_id and family_id and class_id:
                    formatted = f"{ID_id}:{family_id}:{class_id}"
            formatted_id = f"{feature.chrom}|{feature.start}|{feature.end}|{formatted}|{feature.score}|{feature.strand}"
        elif annotation_type == 'gene':  # Gene annotation
            if gene_id and transcript_id and exon_number:
                formatted_id = f"{gene_id}:{transcript_id}:exon_{exon_number}"

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
    bed = gtf.each(format_to_bed)
    sorted_bed = bed.sort()
    sorted_bed.saveas(bed_output_path)

    if os.path.exists(bed_output_path) and os.path.getsize(bed_output_path) == 0:
        log_message("[ERROR]", f"Please check that the reference gene/TE annotation files.", color="error")
        sys.exit(1)

# Merge and sort BED files
def merge_and_sort_bed(bed_file1, bed_file2, output_bed_file):
    bed1 = pybedtools.BedTool(bed_file1)
    if bed_file2 == None:
        merged_bed = bed1.sort()
    else:
        bed2 = pybedtools.BedTool(bed_file2)
        # Concatenate and sort the BED files
        merged_bed = bed1.cat(bed2, postmerge=False).sort()
    
    # Save the sorted BED file
    merged_bed.saveas(output_bed_file)

# Process MSTRG/G.transcripts GTF
def modify_gtf_with_mapping(input_tmap, input_gtf, output_gtf):
    # Generate mapping
    mapping = defaultdict(list)
    with open(input_tmap, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            ref_gene_id = fields[0]
            qry_gene_id = fields[3]
            
            # Only consider MSTRG gene IDs
            if qry_gene_id.startswith("MSTRG") or qry_gene_id.startswith("G"):
                if ref_gene_id == "-":
                    mapping[qry_gene_id] = ['-']  # Map to itself
                elif ref_gene_id not in mapping[qry_gene_id]:
                    mapping[qry_gene_id].append(ref_gene_id)
    
    # Finalize mapping: Select the most frequent ref_gene_id or map to '-'
    for gene, ref_list in mapping.items():
        if len(ref_list) > 1:
            # Count occurrences and select the most frequent ref_gene_id
            most_common = Counter(ref_list).most_common(1)[0][0]
            mapping[gene] = [most_common]
        else:
            mapping[gene] = ref_list  # Keep the single element

    with open(input_gtf, 'r') as in_file, open(output_gtf, 'w') as out_file:
        for line in in_file:
            if line.startswith("#"):
                out_file.write(line)
                continue
            
            fields = line.strip().split('\t')
            gene_id_match = re.search(r'gene_id "([^"]+)"', fields[8])
            
            if gene_id_match:
                gene_id = gene_id_match.group(1)
                
                # If gene_id starts with "MSTRG" or "G", apply the mapping
                if gene_id.startswith("MSTRG") or gene_id.startswith("G"):
                    if gene_id in mapping:
                        new_gene_id = mapping[gene_id][0]  # Take the first mapped ref_gene_id
                        if new_gene_id == "-":
                            continue  # Skip this line if mapped to '-'
                        
                        # Replace the gene_id in the GTF line
                        fields[8] = re.sub(r'gene_id "([^"]+)"', f'gene_id "{new_gene_id}"', fields[8])
            
            # Write the modified line to the output file
            out_file.write('\t'.join(fields) + '\n')

# Extract junctions from BAM files
def extract_junction(input_bam, sample, output_dir, args):
    cmd_header = f'samtools view -H {input_bam} > {output_dir}/{sample}_header.txt'
    subprocess.call(cmd_header, shell=True)

    if args.readtype == 'paired' and args.strand != 'none':
        # For paired-end data, process read1
        outfile1 = os.path.join(output_dir, f'{sample}_read1.bam')
        cmd_junc1 = f"samtools view -@ {args.threads} -f 64 -F 256 {input_bam} | awk '{{if ($6 ~/N/) {{print $0}}}}' | cat {output_dir}/{sample}_header.txt - | samtools view -bS - | samtools sort - -T {outfile1} -o {outfile1}"
        cmd_index1 = f"samtools index -@ {args.threads} {outfile1}"
        for command in (cmd_junc1, cmd_index1):
            subprocess.call(command, shell=True)

        # For paired-end data, process read2
        outfile2 = os.path.join(output_dir, f'{sample}_read2.bam')
        cmd_junc2 = f"samtools view -@ {args.threads} -f 128 -F 256 {input_bam} | awk '{{if ($6 ~/N/) {{print $0}}}}' | cat {output_dir}/{sample}_header.txt - | samtools view -bS - | samtools sort - -T {outfile2} -o {outfile2}"
        cmd_index2 = f"samtools index -@ {args.threads} {outfile2}"
        for command in (cmd_junc2, cmd_index2):
            subprocess.call(command, shell=True)
    elif args.readtype == 'single' or args.strand == 'none':
        # For single-end data
        outfile = os.path.join(output_dir, f'{sample}_read.bam')
        cmd_junc = f"samtools view -@ {args.threads} -F 256 {input_bam} | awk '{{if ($6 ~/N/) {{print $0}}}}' | cat {output_dir}/{sample}_header.txt - | samtools view -bS - | samtools sort - -T {outfile} -o {outfile}"
        cmd_index = f"samtools index -@ {args.threads} {outfile}"
        for command in (cmd_header, cmd_junc, cmd_index):
            subprocess.call(command, shell=True)

    return f'{output_dir}/{sample}'

# Merge classified exons of biological replicates
def merged_exon(classified_exon_files, samples, merged_classified_exons_path, threshold):
    merged_df = pd.DataFrame()

    # Define mapping for splitting ID_position
    split_mapping = {
        "FirstInternal_medium": ["first", "internal"],
        "FirstInternal_high": ["first", "internal"],
        "InternalLast_medium": ["internal", "last"],
        "InternalLast_high": ["internal", "last"]
    }

    # Iterate through all files and process them
    for idx, file in enumerate(classified_exon_files):
        # Load the file
        df = pd.read_csv(file, usecols=['exon', 'gene', 'strand', 'ID_position'], sep='\t', header=0)
        sample_name = samples[idx]
        df[sample_name] = 1

        # Expand rows where ID_position matches split_mapping
        expanded_rows = []
        for _, row in df.iterrows():
            if row['ID_position'] in split_mapping:
                for new_position in split_mapping[row['ID_position']]:
                    new_row = row.copy()
                    new_row['ID_position'] = new_position
                    expanded_rows.append(new_row)
            else:
                expanded_rows.append(row)

        # Convert expanded rows back to DataFrame
        df = pd.DataFrame(expanded_rows)
        
        # Merge with the main DataFrame (on exon, gene, strand, and ID_position)
        if merged_df.empty:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df[['exon', 'gene', 'strand', 'ID_position', sample_name]], 
                                 on=['exon', 'gene', 'strand', 'ID_position'], how='outer')

    # Determine threshold count (minimum number of samples in which the exon should be present)
    threshold_count = threshold * len(samples)

    merged_df['exon_presence_count'] = merged_df[samples].sum(axis=1)
    # filtered_df = merged_df[merged_df['exon_presence_count'] >= threshold_count]

    # filtered_df = filtered_df.copy()
    filtered_df = merged_df[merged_df['exon_presence_count'] >= threshold_count].copy()

    filtered_df.loc[:, 'score'] = filtered_df['exon_presence_count'] / len(samples)
    filtered_df.loc[:, 'chr'] = filtered_df['exon'].str.split(':').str[0]
    filtered_df.loc[:, 'start'] = filtered_df['exon'].str.split(':').str[1].str.split('-').str[0].astype(int)
    filtered_df.loc[:, 'end'] = filtered_df['exon'].str.split(':').str[1].str.split('-').str[1].astype(int)
    filtered_df.loc[:, 'name'] = filtered_df['gene'] + ':' + filtered_df['ID_position'].astype(str)

    # filtered_df = merged_df[merged_df['exon_presence_count'] >= threshold_count].copy()

    # filtered_df['score'] = filtered_df['exon_presence_count'] / len(samples)
    # filtered_df['chr'] = filtered_df['exon'].str.split(':').str[0]
    # filtered_df[['start', 'end']] = filtered_df['exon'].str.split(':').str[1].str.split('-').apply(lambda x: pd.Series([int(x[0]), int(x[1])]))
    # filtered_df['name'] = filtered_df['gene'] + ':' + filtered_df['ID_position'].astype(str)


    # Save the result as a BED file
    bed_df = filtered_df[['chr', 'start', 'end', 'name', 'score', 'strand']]
    bed_df.to_csv(merged_classified_exons_path, sep='\t', header=False, index=False)

    return merged_classified_exons_path

# Process expressed transcripts
def process_expressed_transcripts(gtf_file, output_file, HIT_classified_exons):
    df = pybedtools.BedTool(gtf_file).to_dataframe()

    # Drop rows with invalid seqname
    df = df[df['seqname'].notna()]

    # Remove comment lines if any survived
    df = df[~df['seqname'].astype(str).str.startswith('#')]

    # Parse attributes
    df['gene_id'] = df['attributes'].str.extract(r'gene_id "([^"]+)"')
    df['transcript_id'] = df['attributes'].str.extract(r'transcript_id "([^"]+)"')

    df['cov'] = (
        df['attributes']
        .str.extract(r'cov "([^"]+)"')[0]
        .astype(float)
    )

    # Filter out exons with cov > 1
    df_exons = df[(df['feature'] == 'exon') & (df['cov'] > 0)]
    df_exons['exon_number'] = df['attributes'].str.extract(r'exon_number "([^"]+)"')

    # Function to classify exon types
    def classify_exons(group):
        # Sort exons by position (start)
        group = group.sort_values(by=['start'])

        # Create exon_type list based on the sorted exons
        if len(group) == 1:  # Single exon
            group['exon_type'] = 'single'
        else:
            group['exon_type'] = 'internal'  # Default: 'internal'
            group.iloc[0, group.columns.get_loc('exon_type')] = 'first'  
            group.iloc[-1, group.columns.get_loc('exon_type')] = 'last'

            # Handle strand-specific classification (reverse for negative strand)
            if group['strand'].iloc[0] == '-':
                group = group.iloc[::-1]  # Reverse the order for negative strand
                group.iloc[0, group.columns.get_loc('exon_type')] = 'first'
                group.iloc[-1, group.columns.get_loc('exon_type')] = 'last'

        return group

    # Apply exon classification for each transcript
    df_exons = df_exons.groupby('transcript_id').apply(classify_exons)
    df_exons['HIT_formatted_id'] = df_exons.apply(
        lambda row: f"{row['gene_id']}:{row['exon_type']}",
        axis=1
    )
    df_exons.head()
    df_exons['start'] = df_exons['start'].astype(int) - 1
    df_exons['end'] = df_exons['end'].astype(int)
    df_exons['cov'] = df_exons['cov'].astype(float)
    df_exons.head()

    # df_exons_filter = pd.merge(
    #     df_exons,  
    #     HIT_classified_exons[['seqname', 'start', 'end', 'HIT_formatted_id', 'strand']], 
    #     on=['seqname', 'start', 'end', 'HIT_formatted_id', 'strand'], 
    #     how='inner'
    # )
    # 1. Merge based on HIT_formatted_id first
    df_merged = pd.merge(
        df_exons,  
        HIT_classified_exons[['seqname', 'start', 'end', 'HIT_formatted_id', 'strand']], 
        on=['seqname', 'HIT_formatted_id', 'strand'], 
        how='inner', 
        suffixes=('_exon', '_hit')  # Add suffixes to distinguish columns from both dataframes
    )

    # 2. Apply filtering conditions: df_exons$start >= HIT_classified_exons$start 
    # and df_exons$end <= HIT_classified_exons$end
    df_exons_filter = df_merged[
        (df_merged['start_exon'] >= df_merged['start_hit']) & 
        (df_merged['end_exon'] <= df_merged['end_hit'])
    ].copy()

    # 3. Rename the columns to match the original naming convention
    df_exons_filter.rename(columns={
        'start_exon': 'start', 
        'end_exon': 'end'
    }, inplace=True)

    # Output the filtered result
    df_exons_filter.head()

    df_exons_filter['formatted_id'] = df_exons_filter.apply(
        lambda row: f"{row['gene_id']}:{row['transcript_id']}:exon_{row['exon_number']}:{row['exon_type']}",
        axis=1
    )

    df_exons_filter['start'] = df_exons_filter['start'].astype(int)
    df_exons_filter['end'] = df_exons_filter['end'].astype(int)
    df_exons_filter['cov'] = df_exons_filter['cov'].astype(float)
    df_exons_filter = df_exons_filter.sort_values(by=['transcript_id', 'start'])

    # Save to BED file
    df_exons_filter[['seqname', 'start', 'end', 'formatted_id', 'cov', 'strand']].to_csv(output_file, sep='\t', header=False, index=False)

# Scan for TEs undergoing exonization
def scan_TE_exon(bam_files, samples, te_bed_file, merge_novel_ref_gtf, merged_classified_exons, output, args):
    # Quantify transcripts
    with open(os.path.join(output, 'sample_list.txt'), 'w') as f:
        for bam, sample in zip(bam_files, samples):
            # sample_dir = os.path.join(output, sample)
            # os.makedirs(sample_dir, exist_ok=True)

            cmd = [
                'stringtie', '-e', '-M', str(1.0), '-u',
                '-p', str(args.threads),
                '-G', str(merge_novel_ref_gtf),
                '-o', f'{output}/{sample}_transcripts.gtf',
                '-A', f'{output}/{sample}_gene_abund.tab',
                # '-C', f'{output}/{sample}_known.cov_refs.gtf',
                # '-b', str(sample_dir),
                str(bam)
            ]

            if str(args.strand) == "rf":
                cmd.append('--rf')
            elif str(args.strand) == "fr":
                cmd.append('--fr')

            if args.is_long:
                cmd.extend(['-L', '-a', '3'])
            
            subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            f.write(f"{sample}\t{output}/{sample}_transcripts.gtf\n")
    
    subprocess.call(['prepDE.py', '-i', f'{output}/sample_list.txt',
                     '-g', f'{output}/gene_count_matrix.csv',
                     '-t', f'{output}/transcript_count_matrix.csv'], 
                     stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    
    # Scan exonization TEs
    # Step 1: Extract expressed exons for each sample
    all_expressed_exons = []
    # Load HIT_classified_exons using pybedtools
    HIT_classified_exons_bed = pybedtools.BedTool(merged_classified_exons)
    HIT_classified_exons_df = HIT_classified_exons_bed.to_dataframe(names=['seqname', 'start', 'end', 'HIT_formatted_id', 'score', 'strand'])
    
    with open(os.path.join(output, 'sample_list.txt'), 'r') as f:
        for line in f:
            sample, gtf_file = line.strip().split("\t")
            expressed_exon_file = os.path.join(output, f"{sample}_expressed_exon.bed")
            process_expressed_transcripts(gtf_file, expressed_exon_file, HIT_classified_exons_df)
            all_expressed_exons.append(expressed_exon_file)

    # Step 2: Take intersection of expressed exons from all samples
    all_exons = []
    for bed_file in all_expressed_exons:
        bed = pybedtools.BedTool(bed_file)
        for exon in bed:
            all_exons.append([exon.chrom, exon.start, exon.end, exon.name, exon.strand])

    df = pd.DataFrame(all_exons, columns=["chrom", "start", "end", "formatted_id", "strand"])
    exon_counts = df.groupby(["chrom", "start", "end", "formatted_id", "strand"]).size().reset_index(name="count")
    
    num_samples = len(all_expressed_exons)
    threshold_count = num_samples * args.threshold
    filtered_exons = exon_counts[exon_counts["count"] >= threshold_count]
    
    filtered_exons["score"] = filtered_exons["count"] / len(all_expressed_exons)
    filtered_exons = filtered_exons[["chrom", "start", "end", "formatted_id", "score", "strand"]]
    
    # Save as a BED file
    filtered_exons["exon_number"] = filtered_exons["formatted_id"].str.extract(r":exon_(\d+)")
    filtered_exons["transcript_id"] = filtered_exons["formatted_id"].str.extract(r"^([^:]+:[^:]+)") 
    filtered_exons["exon_number"] = pd.to_numeric(filtered_exons["exon_number"], errors='coerce').fillna(0).astype(int)
    sorted_exons = filtered_exons.sort_values(by=["transcript_id", "exon_number"], ascending=[True, True])  
    sorted_exons = sorted_exons[["chrom", "start", "end", "formatted_id", "score", "strand"]]

    exon_bed_file = os.path.join(output, "intersect_expressed_exon.bed")
    sorted_exons.to_csv(exon_bed_file, sep="\t", header=False, index=False)

    # Step 3: Intersect with TE.bed
    te_bed = pybedtools.BedTool(te_bed_file)
    filtered_bed = pybedtools.BedTool(exon_bed_file)
    overlaps = filtered_bed.intersect(te_bed, wa=True, wb=True, s=True)

    # Save total result to exon_TE.bed
    total_output = os.path.join(os.path.dirname(output), "exon_TE.bed")
    overlaps.saveas(total_output)
    
    # # Step 4: Classify overlapped TEs
    # output_paths = {
    #     "first": os.path.join(os.path.dirname(output), "first_TE.bed"),
    #     "last": os.path.join(os.path.dirname(output), "last_TE.bed"),
    #     "internal": os.path.join(os.path.dirname(output), "internal_TE.bed")
    # }

    # # Create dictionaries to store results and transcripts
    # results = {key: [] for key in output_paths}
    # transcripts = {key: [] for key in output_paths}

    # # Process the overlaps
    # for line in overlaps:
    #     fields = line.fields
    #     exon_name = fields[3]
    #     classified_result = exon_name.split(":")[3]
    #     classified_trans = exon_name.split(":")[1]
        
    #     # Append the result and transcript based on classification
    #     if classified_result in results:
    #         results[classified_result].append("\t".join(fields))
    #         transcripts[classified_result].append(classified_trans)
    
    # # Write the results to files (TE classification)
    # gtf_output_path = []
    # for classification, output_path in output_paths.items():
    #     # Skip empty results
    #     if not results[classification]:
    #         continue

    #     with open(output_path, "w") as file:
    #         file.write("\n".join(results[classification]) + "\n")
        
    #     # Save the transcript IDs for future GTF extraction
    #     transcript_list = transcripts[classification]
    #     transcript_file = os.path.join(output, f"{classification}_transcripts.txt")
    #     with open(transcript_file, "w") as tf:
    #         tf.write("\n".join(transcript_list))

    #     # Use fgrep to extract the corresponding lines from the original GTF file
    #     gtf_input = merge_novel_ref_gtf
    #     gtf_output = os.path.join(os.path.dirname(output), f"{classification}_classified.gtf")
    #     gtf_output_path.append(gtf_output)
    #     fgrep_command = f"fgrep -f {transcript_file} {gtf_input} > {gtf_output}"
    #     subprocess.run(fgrep_command, shell=True, check=True)
        
    #     # Extract sequences using gffread and the classified GTF file
    #     genome_fasta = args.genome
    #     fasta_output = os.path.join(os.path.dirname(output), f"{classification}_extracted_sequences.fasta")
    #     gffread_command = f"gffread {gtf_output} -g {genome_fasta} -w {fasta_output}"
    #     subprocess.run(gffread_command, shell=True, check=True)

    # return gtf_output_path

def generate_metaexon_gtf(sample_list, candidate_dir, out_dir, project):
    line_counter = Counter()

    for sample in sample_list:
        exon_bed = os.path.join(candidate_dir, sample, "exon_TE.bed")
        if os.path.isfile(exon_bed):
            with open(exon_bed) as f:
                unique_lines = set(line.strip() for line in f)
                line_counter.update(unique_lines)
        else:
            print(f"⚠️ Skipping missing file: {exon_bed}")

    num_samples = len(sample_list)
    shared_lines = [line for line, count in line_counter.items() if count == num_samples]

    with NamedTemporaryFile(mode="w+", delete=False) as temp_bed:
        for line in shared_lines:
            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue
            selected_fields = fields[:6] + [fields[9]]
            temp_bed.write("\t".join(selected_fields) + "\n")
        temp_bed_path = temp_bed.name

    # Sort using bedtools
    # sorted_bed_path = temp_bed_path + ".sorted"
    sorted_bed_path = os.path.join(out_dir, f"{project}_TEexon.bed")
    with open(sorted_bed_path, "w") as sorted_bed:
        subprocess.run(["bedtools", "sort", "-i", temp_bed_path], stdout=sorted_bed)

    # Merge using bedtools
    metaexon_bed = os.path.join(out_dir, f"{project}_metaexon.bed")
    merge_cmd = [
        "bedtools", "merge", "-s",
        "-c", "4,6,7", "-o", "distinct",
        "-i", sorted_bed_path
    ]
    with open(metaexon_bed, "w") as out_bed:
        subprocess.run(merge_cmd, stdout=out_bed)

    gtf_out = os.path.join(out_dir, f"{project}_metaexons.gtf")
    with open(metaexon_bed) as f_in, open(gtf_out, "w") as f_out:
        for line in f_in:
            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue
            chrom, start, end, names, strand, te = fields[0], str(int(fields[1])+1), fields[2], fields[3], fields[4], fields[5]
            exon_id = names.split(",")[0]
            gene_id = exon_id.split(":")[0]
            attr = f'gene_id "{gene_id}"; exon_id "{names}"; te_id "{te}"'
            f_out.write(f"{chrom}\t{project}\texon\t{start}\t{end}\t.\t{strand}\t.\t{attr}\n")

    os.remove(temp_bed_path)

    return gtf_out

# Generate corrected TE-exon expression matrix
def generate_TEexon(metaexon_gtf, featurecount_exon_txt, tecount_file_list):

    # 1. Load GTF and extract metaexon ↔ TE ID mapping
    gtf_cols = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    gtf_df = pd.read_csv(metaexon_gtf, sep='\t', names=gtf_cols, comment='#')

    gtf_df['metaexon'] = gtf_df['chr'].astype(str) + ':' + gtf_df['start'].astype(str) + '-' + gtf_df['end'].astype(str) + ':' + gtf_df['strand']
    gtf_df['TEid'] = gtf_df['attribute'].str.extract(r'te_id "([^"]+)"')

    # Handle multiple TE IDs per exon
    gtf_df = gtf_df.copy()
    gtf_df = gtf_df.assign(TEid=gtf_df['TEid'].str.split(',')).explode('TEid').reset_index(drop=True)
    gtf_df['TEid'] = gtf_df['TEid'].str.replace(r'^((?:[^|]*\|){4})[^|]*(\|[^|]*)$', r'\1.\2', regex=True)

    gene_te_map = gtf_df[['metaexon', 'TEid']].dropna()

    # 2. Load exon count table and construct metaexon column
    exon_df_raw = pd.read_csv(featurecount_exon_txt, sep='\t', comment='#')
    sample_names = [col.split('/')[-1].replace('_accepted_hits.bam', '') for col in exon_df_raw.columns if col.endswith('.bam')]
    exon_df_raw.columns = exon_df_raw.columns[:6].tolist() + sample_names

    exon_df_raw['metaexon'] = exon_df_raw['Chr'].astype(str) + ':' + exon_df_raw['Start'].astype(str) + '-' + exon_df_raw['End'].astype(str) + ':' + exon_df_raw['Strand']
    exon_df = exon_df_raw[['metaexon'] + sample_names].copy()
    exon_df.set_index('metaexon', inplace=True)

    # 3. Load TE count files for each sample
    tecount_dict = {}
    for f in tecount_file_list:
        sample = f.split('/')[-1].replace('_TEcounts.txt', '')
        te_df = pd.read_csv(f, sep='\t')
        te_df = te_df[['TE_ID', 'uniq_counts', 'tot_counts']].copy()
        te_df.set_index('TE_ID', inplace=True)
        tecount_dict[sample] = te_df

    # 4. Compute corrected TEexon counts
    teexon_df = pd.DataFrame(index=exon_df.index)

    for sample in sample_names:
        if sample not in tecount_dict:
            print(f"Warning: {sample} TEcount not found, skipping.")
            continue

        exon_counts = exon_df[sample]
        te_counts = tecount_dict[sample]

        # Merge TE counts with gene mapping and aggregate
        merged = gene_te_map.merge(te_counts, left_on='TEid', right_index=True, how='left')
        print(merged.head())
        merged_grouped = merged.groupby('metaexon')[['uniq_counts', 'tot_counts']].sum()

        # Correct exon counts: TEexon = exon - sum(uniq_counts) + sum(tot_counts)
        correction = -merged_grouped['uniq_counts'].fillna(0) + merged_grouped['tot_counts'].fillna(0)
        correction = correction.reindex(exon_counts.index).fillna(0)

        teexon_df[sample] = exon_counts + correction

    return teexon_df

# Refine DESeq2 output by integrating metaexon BED information.
def refine_result(deseq_out, metaexon_bed):
    deseq_df = pd.read_csv(deseq_out, sep="\t")

    deseq_df["coord"] = deseq_df.index

    bed_cols = ["chr", "start", "end", "gene_info", "strand", "TE_info"]
    bed_df = pd.read_csv(metaexon_bed, sep="\t", header=None, names=bed_cols)

    bed_df["coord"] = bed_df["chr"] + ":" + (bed_df["start"]+1).astype(str) + "-" + bed_df["end"].astype(str) + ":" + bed_df["strand"]

    def extract_genes(gene_str):
        genes = [g.split(":")[0] for g in gene_str.split(",")]
        return ",".join(sorted(set(genes)))

    bed_df["genes"] = bed_df["gene_info"].apply(extract_genes)

    bed_clean = bed_df[["coord", "genes", "TE_info"]]

    merged = pd.merge(deseq_df, bed_clean, on="coord", how="left")
    merged.set_index("coord", inplace=True)
    merged.to_csv(deseq_out, sep="\t", index=True)

# Integrate PLEK predictions, transcript info, metaexon info, and TE annotation into a single CSV.
def generate_plek_annotation_result(plek_result_file, transcript_file, meta_sig, te_bed_path, output_file):

    # Load PLEK predictions and transcript IDs
    with open(plek_result_file, 'r') as f:
        plek_results = [line.strip() for line in f]

    with open(transcript_file, 'r') as f:
        all_transcripts = [line.strip() for line in f]

    if len(plek_results) != len(all_transcripts):
        raise ValueError("Mismatch between PLEK results and transcript list")

    # Parse meta_sig exon_info into transcript-level mapping
    exon_map = {}
    for exon_info, metaexon in zip(meta_sig["exon_info"], meta_sig["coord"]):
        for item in exon_info.split(","):
            try:
                gene, transcript, exon, loc = item.split(":")
                exon_id = f"{exon}:{loc}"
                full_id = f"{gene}:{transcript}:{exon}:{loc}"
                if transcript not in exon_map:
                    exon_map[transcript] = {
                        "Gene": set(),
                        "DiffExon": set(),
                        "MetaExon": set(),
                        "FullExonIDs": set()
                    }
                exon_map[transcript]["Gene"].add(gene)
                exon_map[transcript]["DiffExon"].add(exon_id)
                exon_map[transcript]["MetaExon"].add(metaexon)
                exon_map[transcript]["FullExonIDs"].add(full_id)
            except ValueError:
                continue

    # Parse TE BED annotation
    te_map = {}
    with open(te_bed_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 7:
                full_exon_id = parts[3]
                te_name = parts[6]
                te_map.setdefault(full_exon_id, set()).add(te_name)

    # Compile integrated results
    records = []
    for tx, pred in zip(all_transcripts, plek_results):
        gene = ",".join(sorted(exon_map.get(tx, {}).get("Gene", [])))
        diff_exon = ",".join(sorted(exon_map.get(tx, {}).get("DiffExon", [])))
        meta_exon = ",".join(sorted(exon_map.get(tx, {}).get("MetaExon", [])))

        # Aggregate TE names
        full_ids = exon_map.get(tx, {}).get("FullExonIDs", [])
        te_list = set()
        for full_id in full_ids:
            te_list.update(te_map.get(full_id, []))
        te_str = ",".join(sorted(te_list))

        records.append({
            "Transcript": tx,
            "Prediction": pred,
            "Gene": gene,
            "DiffExon": diff_exon,
            "MetaExon": meta_exon,
            "TE": te_str
        })

    result_df = pd.DataFrame(records)
    result_df.to_csv(output_file, index=False)
    print(f"Integrated annotation saved to: {output_file}")

    return result_df

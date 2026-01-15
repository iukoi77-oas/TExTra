import os
import pybedtools
import collections
from util.toolkit import *
from util.HITindex_toolkit import *

# Parse raw exon bed
def parse_bed(bedfile):
    # dictionary to store genes
    # structure: genedict[gene_name][transcript_name][exon_#][exon_coords]
    genedict = collections.defaultdict(lambda: collections.defaultdict(list))

    with open(bedfile, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom, start, end, name, score, strand = fields[:6]
            parts = name.split(":")
            
            if len(parts) == 3:
                # process gene:transcript:exon
                gene, transcript, exon = parts[0], parts[1], parts[2]
                genedict[gene][transcript].append((chrom, int(start), int(end), strand, exon))
            else:
                continue 
    
    return genedict

# Merge exons into metaexon
# def metaexon_bed(input_bed, outdir, args):
#     # Output file paths
#     outfh_path = os.path.join(outdir, 'metaexon.bed')
#     outfhbuffer_path = os.path.join(outdir, f'buffer_ss3-{args.ss3buffer}_ss5-{args.ss5buffer}_metaexon.bed')
#     outfhconst_path = os.path.join(outdir, 'constituent_metaexon.bed')

#     # Open output files for writing
#     with open(outfh_path, 'w') as outfh, open(outfhbuffer_path, 'w') as outfhbuffer, open(outfhconst_path, 'w') as outfhconst:
#         # Parse input BED file into a gene dictionary
#         genedict = parse_bed(input_bed)

#         # List of possible exon types
#         exontypes = ['FE', 'internal', 'LE', 'singleexon']
        
#         # Iterate over each gene in the dictionary
#         for gene, transcripts in genedict.items():
#             exons = []  # List to store exons
#             exons_const = []  # List to store constituent exons
#             ntxpt = 'TXPT:' + str(len(genedict[gene]))

#             # Iterate over each transcript for the current gene
#             for transcript, transcript_exons in transcripts.items():
#                 # Sort exons by their start position
#                 exons_sorted = sorted(transcript_exons, key=lambda x: x[1])
#                 strand = exons_sorted[0][3]  # Get strand direction (assume all exons in the same transcript have the same strand)

#                 # Determine exon types based on strand and exon count
#                 if len(exons_sorted) == 1:
#                     exon_types = ['singleexon']
#                 else:
#                     if strand == '+':
#                         exon_types = ['FE'] + ['internal'] * (len(exons_sorted) - 2) + ['LE']
#                     elif strand == '-':
#                         exon_types = ['LE'] + ['internal'] * (len(exons_sorted) - 2) + ['FE']
#                     else:
#                         continue  # Skip if strand is undefined

#                 # Process each exon for the current transcript
#                 for i, exon in enumerate(exons_sorted):
#                     chrom, start, end, strand, exon_id = exon  # Extract exon information
#                     exon_type = exon_types[i]  # Get the corresponding exon type
#                     exons.append(f"{chrom} {start} {end} {exon_type} . {strand}")  # Add to exons list
#                     exons_const.append(f"{chrom} {start} {end} {start}-{end} . {strand}")  # Add to constituent exons list
            
#             # Merge overlapping exons into single exons, while collapsing exon type information
#             bedscratch = pybedtools.BedTool('\n'.join(exons), from_string=True)
#             bedsort = bedscratch.sort()  # Sort exons
#             bedsortmerge = bedsort.merge(s=True, c=4, o='collapse')  # Merge overlapping exons, collapsing exon type counts
#             bedlist = str(bedsortmerge).split('\n')[:-1]  # Get the merged exon list

#             # Same process for the constituent exons
#             bedscratch_const = pybedtools.BedTool('\n'.join(exons_const), from_string=True)
#             bedsort_const = bedscratch_const.sort()
#             bedsortmerge_const = bedsort_const.merge(s=True, c=4, o='collapse')
#             bedlist_const = str(bedsortmerge_const).split('\n')[:-1]

#             # Write out the merged exons
#             for ex in bedlist:
#                 exhere = ex.split('\t')  # Split each line into components
#                 exherename = f"{exhere[0]}:{exhere[1]}-{exhere[2]}"  # Create the exon name (e.g., chrom:start-end)
#                 extypes = exhere[3].split(',')  # Get the types of exons for counting
#                 ntypeset = []  # List to store counts of exon types
                
#                 # Count the occurrence of each exon type
#                 for x in exontypes:
#                     ntypeset.append(f"{x}:{extypes.count(x)}")
                
#                 # Calculate buffered regions based on strand direction
#                 if strand == '+':
#                     exstart = int(exhere[1]) - args.ss3buffer  # 3' end buffer
#                     exend = int(exhere[2]) + args.ss5buffer  # 5' end buffer
#                 elif strand == '-':
#                     exstart = int(exhere[1]) - args.ss5buffer  # 5' end buffer for negative strand
#                     exend = int(exhere[2]) + args.ss3buffer  # 3' end buffer for negative strand
                
#                 # Write the merged exons to the output file
#                 outfh.write(f"{exhere[0]}\t{exhere[1]}\t{exhere[2]}\t{exherename};{gene};{ntxpt};{';'.join(ntypeset)}\t0\t{strand}\n")
#                 outfhbuffer.write(f"{exhere[0]}\t{exstart}\t{exend}\t{exherename};{gene};{ntxpt};{';'.join(ntypeset)}\t0\t{strand}\n")

#             # Write out the constituent exons (non-merged)
#             for ex in bedlist_const:
#                 exhere = ex.split('\t')
#                 exherename = f"{exhere[0]}:{exhere[1]}-{exhere[2]}"
#                 exactual = exhere[3]  # Exon type (e.g., 'FE', 'LE', etc.)
#                 # Write the constituent exons to the corresponding file
#                 outfhconst.write(f"{exhere[0]}\t{exhere[1]}\t{exhere[2]}\t{exherename};{gene};{ntxpt};{exactual}\t0\t{strand}\n")
    
#     return outfhbuffer_path

def metaexon_bed(input_bed, outdir, args):
    # Define output file paths
    outfh_path = os.path.join(outdir, 'metaexon.bed')  # File for merged exons
    outfhbuffer_path = os.path.join(outdir, f'buffer_ss3-{args.ss3buffer}_ss5-{args.ss5buffer}_metaexon.bed')  # File for buffered exons
    outfhconst_path = os.path.join(outdir, 'constituent_metaexon.bed')  # File for non-merged exons

    # Open output files for writing
    with open(outfh_path, 'w') as outfh, open(outfhbuffer_path, 'w') as outfhbuffer, open(outfhconst_path, 'w') as outfhconst:
        # Parse the input BED file into a dictionary where keys are genes and values are transcript-to-exon mappings
        genedict = parse_bed(input_bed)

        # Define the possible exon types
        exontypes = ['FE', 'internal', 'LE', 'singleexon']

        # Iterate over each gene in the dictionary
        for gene, transcripts in genedict.items():
            exons = []  # List to store exon information for merging
            exons_const = []  # List to store raw exon information
            ntxpt = 'TXPT:' + str(len(genedict[gene]))  # Number of transcripts for this gene

            # Iterate over transcripts of the current gene
            for transcript, transcript_exons in transcripts.items():
                # Sort exons based on start position
                exons_sorted = sorted(transcript_exons, key=lambda x: x[1])
                strand = exons_sorted[0][3]  # Get strand information (assumed to be the same for all exons in a transcript)

                # Determine exon types based on strand and exon count
                if len(exons_sorted) == 1:
                    exon_types = ['singleexon']
                else:
                    if strand == '+':
                        exon_types = ['FE'] + ['internal'] * (len(exons_sorted) - 2) + ['LE']
                    elif strand == '-':
                        exon_types = ['LE'] + ['internal'] * (len(exons_sorted) - 2) + ['FE']
                    else:
                        continue  # Skip undefined strand cases

                # Process each exon
                for i, exon in enumerate(exons_sorted):
                    chrom, start, end, strand, exon_id = exon  # Extract exon details
                    exon_type = exon_types[i]  # Assign exon type
                    exons.append(f"{chrom} {start} {end} {exon_type} . {strand}")  # Store exon in list
                    exons_const.append(f"{chrom} {start} {end} {start}-{end} . {strand}")  # Store raw exon details

            # Merge overlapping exons while preserving exon type information
            bedscratch = pybedtools.BedTool('\n'.join(exons), from_string=True)
            bedsortmerge = bedscratch.sort().merge(s=True, c=4, o='collapse')  # Sort and merge overlapping exons
            bedlist = str(bedsortmerge).split('\n')[:-1]  # Convert merged exons into a list

            # Process constituent exons without merging
            bedscratch_const = pybedtools.BedTool('\n'.join(exons_const), from_string=True)
            bedsortmerge_const = bedscratch_const.sort().merge(s=True, c=4, o='collapse')
            bedlist_const = str(bedsortmerge_const).split('\n')[:-1]

            # Write merged exons to output file
            for ex in bedlist:
                exhere = ex.split('\t')  # Split into components
                exherename = f"{exhere[0]}:{exhere[1]}-{exhere[2]}"  # Format exon name
                extypes = exhere[3].split(',')  # Get list of exon types
                ntypeset = [f"{x}:{extypes.count(x)}" for x in exontypes]  # Count occurrences of each exon type

                # Adjust exon start and end positions based on strand direction and buffer values
                if strand == '+':
                    exstart = max(0, int(exhere[1]) - args.ss3buffer)  # Buffer 3' end
                    exend = int(exhere[2]) + args.ss5buffer  # Buffer 5' end
                elif strand == '-':
                    exstart = max(0, int(exhere[1]) - args.ss5buffer)  # Buffer 5' end for negative strand
                    exend = int(exhere[2]) + args.ss3buffer  # Buffer 3' end for negative strand
                
                # Write merged exon to the main output file
                outfh.write(f"{exhere[0]}\t{exhere[1]}\t{exhere[2]}\t{exherename};{gene};{ntxpt};{';'.join(ntypeset)}\t0\t{strand}\n")
                outfhbuffer.write(f"{exhere[0]}\t{exstart}\t{exend}\t{exherename};{gene};{ntxpt};{';'.join(ntypeset)}\t0\t{strand}\n")

            # Write out the constituent exons (non-merged)
            for ex in bedlist_const:
                exhere = ex.split('\t')
                exherename = f"{exhere[0]}:{exhere[1]}-{exhere[2]}"
                exactual = exhere[3]  # Exon type (e.g., 'FE', 'LE', etc.)
                outfhconst.write(f"{exhere[0]}\t{exhere[1]}\t{exhere[2]}\t{exherename};{gene};{ntxpt};{exactual}\t0\t{strand}\n")

    return outfhbuffer_path  # Return the buffered exon file path


# Run the HITindex pipeline
def HITindex_pipeline(input_buffer_bed, bamfiles, samples, output_hitindex_dir, args):
    tmp_dir = os.path.join(output_hitindex_dir, 'tmp')
    os.makedirs(tmp_dir, exist_ok=True)

    script_dir = os.path.dirname(os.path.realpath(__file__))
    paramfile_path = os.path.join(script_dir, 'HIT_identity_parameters.txt')
    exon_outnames = []
    
    # Loop through BAM files and corresponding sample names
    for bam, sample in zip(bamfiles, samples):
        # Calculating HITindex metrics
        juncbam = extract_junction(bam, sample, tmp_dir, args)
        print(f'============ Processing for {sample} ============')
        beddict = getExons(input_buffer_bed)
        getJuncBed(input_buffer_bed, juncbam, args.readtype, args.strand)
        startdict, enddict = getExonReads(juncbam, 10)
        HITdict = calculateMetric(beddict, startdict, enddict, 2)
        HITdict, param = edge_flagging(HITdict)
        exon_outname = f'{output_hitindex_dir}/{sample}.exon'
        writeMetrics(HITdict, param, exon_outname)
        running_genmodel(exon_outname)

        # Exon Classification
        HITcombo = readMetrics(exon_outname)
        paramdict = readParameters(paramfile_path)
        HITcombo = significance(HITcombo, 1000, float(paramdict['HIThybrid']), exon_outname)
        call_terminal(HITcombo, paramdict, exon_outname)

        exon_outnames.append(exon_outname)

    return exon_outnames
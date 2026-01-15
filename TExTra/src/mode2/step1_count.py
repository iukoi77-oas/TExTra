import os
import re
import glob
import tempfile
import subprocess as sp
from datetime import datetime

from util.logger import *
from util.SQuIRE_toolkit import *
from util.test import match_reads_python

strand_map = {
    "none": 0,
    "rf": 1,  # First-strand
    "fr": 2,  # Second-strand
    "r": 1,   # First-strand
    "f": 2    # Second-strand
}

def count_func(rmsk_bed, args, sample, replicate):
    tempfolder=args.tempfolder
    strandedness=strand_map.get(args.strand.lower(), 0)  # default: 0（unstranded）
    EM=args.EM
    verbosity=True
    debug=False
        
    count_dir = os.path.join(args.out_dir, "quantification")
    os.makedirs(count_dir, exist_ok=True)

    align_dir = os.path.join(args.prep, f"alignment/{replicate}")
    logfile = find_file(align_dir,".final.out", replicate, 1, False)
    bamfile = find_file(align_dir,"_accepted_hits.bam", replicate, 1, True)
    if not bamfile:
        if replicate:
            raise Exception("Cannot find bamfile matching " + replicate )
        else:
            raise Exception("Cannot find bamfile in map_folder" )
    if not replicate:
        filename = os.path.basename(bamfile)
        replicate = os.path.splitext(filename)[0]

    # rmsk_bed = '/home/yanj/TE_exon/script/test1/clean/GRCh38_GENCODE_rmsk_TE.bed'
    # clean_out = os.path.join(args.out_dir, "clean")
    convert_out = os.path.join(args.prep, "convert")
    # rmsk_bed = find_file(convert_out, "TE_anno.bed", False, 1, True)
    # copies = find_file(clean_out, "_copies.txt", False, 1, True)
    if not rmsk_bed:
        raise Exception("Cannot find bedfile in clean_folder" )
    # if not copies:
    #     raise Exception("Cannot find copies.txt file in clean_folder")

    if not tempfolder:
        tempfolder = os.path.join(count_dir, "tmp")
        os.makedirs(tempfolder, exist_ok=True)

    paired_end = True if args.readtype == "paired" else False

    outgtf_ref = os.path.join(count_dir, f"{replicate}.gtf")
    abund_ref = outgtf_ref.replace(".gtf","_abund.txt")
    stringtie_dir = os.path.join(args.qual, f"candidate/{sample}/tmp")
    outgtf_ref_temp = find_file(stringtie_dir, "_transcripts.gtf", replicate, 1, True)
    abund_ref_temp = find_file(stringtie_dir, "_gene_abund.tab", replicate, 1, True)

    outgtf_copy = shutil.copy(outgtf_ref_temp, tempfolder)
    abund_ref_copy = shutil.copy(abund_ref_temp, tempfolder)

    sort_coord(outgtf_copy,outgtf_ref,1,4,debug)
    sort_coord_header(abund_ref_copy,abund_ref,3,5,debug)	       
    genename_dict={}
    filter_abund(abund_ref,genename_dict,False)
    genecounts=count_dir + "/" + replicate + "_refGenecounts.txt"

    cmd = f"samtools view {bamfile} | awk '{{print length($10)}}' | head -100"
    process = sp.run(cmd, shell=True, capture_output=True, text=True)

    read_lengths = [int(line) for line in process.stdout.strip().split("\n") if line.isdigit()]
    if read_lengths:
        readlength = max(set(read_lengths), key=read_lengths.count)
        print(f"Most common read length: {readlength}")
    else:
        readlength = None
        print("No reads found.")

    filter_tx(outgtf_ref, genename_dict,readlength,genecounts)

    #### OPEN OUTPUTS & WRITE HEADER INFORMATION#############
    if verbosity:
        print("Creating temporary files "+ str(datetime.now()) ,file = sys.stderr)
    counts_temp = tempfile.NamedTemporaryFile(delete=False, dir = tempfolder, prefix="count" +  ".tmp")
    countsfilepath = count_dir + "/" + replicate + "_TEcounts.txt"
    counts_file_header = open(countsfilepath +".header",'w')

    counts_file_header.writelines("tx_chr" + "\t"  + "tx_start"  + "\t" + "tx_stop"  + "\t" + "TE_ID" + "\t" + "fpkm"  + "\t" + "tx_strand" + "\t" + "Sample" + "\t" + "alignedsize" + "\t" + "TE_chr" + "\t" + "TE_start" + "\t" + "TE_stop" + "\t" + "TE_name" + "\t" + "milliDiv" + "\t" + "TE_strand" + "\t" + "uniq_counts" + "\t" + "tot_counts"  + "\t" + "tot_reads" +"\t" + "score" + "\n" )

    counts_file_header.close()


#####CREATE TEMPFILES #######
    if verbosity:
        print("Creating unique and multiple alignment bedfiles "+ str(datetime.now()) ,file = sys.stderr)

    if not paired_end:
        single_bam = bamfile
        if verbosity:
            print("Intersecting bam file with TE bedfile "+ str(datetime.now()) ,file = sys.stderr)
        #intersect bam files with TE bed files
        single_bed_tempfile1 = make_tempfile(replicate,"single_bed_1", tempfolder)
        intersect_flank(single_bam, rmsk_bed, single_bed_tempfile1,debug)
        if verbosity:
            print("Combining adjacent TEs with same read alignment "+ str(datetime.now()) ,file = sys.stderr)
        #reduce reads   #Find reads aligned to same position but different TE_IDs (overlapping flanks) and merge
        single_reduced_tempfile1 = make_tempfile(replicate,"single_reduced_1", tempfolder)
        single_reduced_tempfile1_sorted =single_reduced_tempfile1  + "_sorted"
        sort_coord(single_bed_tempfile1,single_reduced_tempfile1_sorted,1,2,debug)
        reduce_reads(single_reduced_tempfile1_sorted, single_reduced_tempfile1,debug)
        if verbosity:
            print("Getting genomic coordinates of read"+ str(datetime.now()) ,file = sys.stderr)
        #get genomic coordinates and RNA strand for all alignments
        single_coords_tempfile1= make_tempfile(replicate,"single_coords_1", tempfolder)
        get_coords(single_reduced_tempfile1,1,strandedness,single_coords_tempfile1,debug)

        # os.unlink(single_bed_tempfile1)

        if verbosity:
            print("Identifying and labeling unique and multi reads"+ str(datetime.now()) ,file = sys.stderr)
        single_labeled_tempfile1 = make_tempfile(replicate,"single_labeled_1", tempfolder)
        single_labeled_tempfile2 = make_tempfile(replicate,"single_labeled_2", tempfolder)
        label_files(single_coords_tempfile1,single_labeled_tempfile1,"single",debug)
        label_files(single_labeled_tempfile1,single_labeled_tempfile2,"R1",debug)

        #find unique single alignments
        first_tempfile1 = make_tempfile(replicate,"first_1", tempfolder)
        unique_tempfile1 = make_tempfile(replicate,"unique_1", tempfolder)
        multi_tempfile1 = make_tempfile(replicate,"multi_1", tempfolder)

        find_uniq(single_labeled_tempfile2,first_tempfile1,unique_tempfile1, multi_tempfile1,debug)


        #label uniq, multi, or single
        multi_bed = make_tempfile(replicate,"multi_bed", tempfolder)
        unique_bed = make_tempfile(replicate,"unique_bed", tempfolder)

        label_files(unique_tempfile1, unique_bed, "uniq",debug)
        label_files(multi_tempfile1, multi_bed, "multi",debug)

        aligned_libsize = getlibsize(logfile, bamfile,multi_bed,unique_bed,paired_end,debug)

    if paired_end:
        #intersect bam files with TE bed files
        if verbosity:
            print("Identifying properly paired reads "+ str(datetime.now()) ,file = sys.stderr)
        paired_bam = bamfile
        proper_bam = make_tempfile(replicate,"proper_bam", tempfolder)
        nonproper_bam = make_tempfile(replicate,"nonproper_bam", tempfolder)
        find_properpair(paired_bam, proper_bam,nonproper_bam)

        if verbosity:
            print("Intersecting bam files with TE bedfile "+ str(datetime.now()) ,file = sys.stderr)

        proper_bed = make_tempfile(replicate,"proper_bed", tempfolder)
        nonproper_bed = make_tempfile(replicate,"nonproper_bed", tempfolder)
        intersect_flank(proper_bam, rmsk_bed, proper_bed,debug)
        intersect_flank(nonproper_bam, rmsk_bed, nonproper_bed,debug)

        proper_labeled_tempfile = make_tempfile(replicate,"proper_labeled", tempfolder)
        nonproper_labeled_tempfile = make_tempfile(replicate,"nonproper_labeled", tempfolder)
        label_files(proper_bed,proper_labeled_tempfile,"proper",debug)
        label_files(nonproper_bed,nonproper_labeled_tempfile,"nonproper",debug)

        proper_nonproper_labeled_tempfile = make_tempfile(replicate,"proper_nonproper_labeled", tempfolder)
        combine_files(proper_labeled_tempfile,nonproper_labeled_tempfile,proper_nonproper_labeled_tempfile,debug)

        if verbosity:
            print("Splitting into read1 and read 2 "+ str(datetime.now()) ,file = sys.stderr)
        paired_bed_tempfile1 = make_tempfile(replicate,"paired_1.bed",tempfolder)
        paired_bed_tempfile2 = make_tempfile(replicate,"paired_2.bed",tempfolder)
        split_paired(proper_nonproper_labeled_tempfile,paired_bed_tempfile1,paired_bed_tempfile2,debug)
        paired_bed_tempfile1_sorted = paired_bed_tempfile1 + "_sorted"
        paired_bed_tempfile2_sorted = paired_bed_tempfile2 + "_sorted"
        if not debug:
            os.unlink(proper_bam)
            os.unlink(nonproper_bam)
        #reduce reads   #Find reads aligned to same position but different TE_IDs (overlapping flanks) and merge
        if verbosity:
            print("Combining adjacent TEs with same read alignment "+ str(datetime.now()) ,file = sys.stderr)
        paired_reduced_tempfile1 = make_tempfile(replicate,"paired_reduced_1", tempfolder)
        paired_reduced_tempfile2 = make_tempfile(replicate,"paired_reduced_2", tempfolder)
        sort_coord(paired_bed_tempfile1,paired_bed_tempfile1_sorted,1,2,debug)
        sort_coord(paired_bed_tempfile2,paired_bed_tempfile2_sorted,1,2,debug)

        reduce_reads(paired_bed_tempfile1_sorted, paired_reduced_tempfile1,debug)
        reduce_reads(paired_bed_tempfile2_sorted, paired_reduced_tempfile2,debug)

        #get genomic coordinates and RNA strand for all alignments
        if verbosity:
            print("Getting genomic coordinates of read "+ str(datetime.now()) ,file = sys.stderr)
        paired_coords_tempfile1= make_tempfile(replicate,"paired_coords_1", tempfolder)
        paired_coords_tempfile2= make_tempfile(replicate,"paired_coords_2", tempfolder)
        get_coords(paired_reduced_tempfile1,1,strandedness,paired_coords_tempfile1,debug)
        get_coords(paired_reduced_tempfile2,2,strandedness,paired_coords_tempfile2,debug)

        paired_labeled_tempfile1 = make_tempfile(replicate,"paired_labeled_1", tempfolder)
        paired_labeled_tempfile2 = make_tempfile(replicate,"paired_labeled_2", tempfolder)
        label_files(paired_coords_tempfile1,paired_labeled_tempfile1,"R1",debug)
        label_files(paired_coords_tempfile2,paired_labeled_tempfile2,"R2",debug)

        #remove /1 and /2 from read ID column
        paired_fixed_tempfile1 = make_tempfile(replicate,"paired_fixed_1", tempfolder)
        paired_fixed_tempfile2 = make_tempfile(replicate,"paired_fixed_2", tempfolder)
        fix_paired(paired_labeled_tempfile1,paired_labeled_tempfile2, paired_fixed_tempfile1,paired_fixed_tempfile2,debug)

        #find unique single alignments
        if verbosity:
            print("Identifying and labeling unique and multi reads "+ str(datetime.now()) ,file = sys.stderr)
        first_tempfile1 = make_tempfile(replicate,"first_1", tempfolder)
        unique_tempfile1 = make_tempfile(replicate,"unique_1", tempfolder)
        multi_tempfile1 = make_tempfile(replicate,"multi_1", tempfolder)


        first_tempfile2 = make_tempfile(replicate,"first_2", tempfolder)
        unique_tempfile2 = make_tempfile(replicate,"unique_2", tempfolder)
        multi_tempfile2 = make_tempfile(replicate,"multi_2", tempfolder)

        find_uniq(paired_fixed_tempfile1,first_tempfile1,unique_tempfile1, multi_tempfile1,debug)
        find_uniq(paired_fixed_tempfile2,first_tempfile2, unique_tempfile2, multi_tempfile2,debug)

        #label uniq, multi, or paired

        unique_tempfile1_labeled = make_tempfile(replicate,"unique_labeled_1", tempfolder)
        multi_tempfile1_labeled = make_tempfile(replicate,"multi_labeled_1", tempfolder)

        unique_tempfile2_labeled = make_tempfile(replicate,"unique_labeled_2", tempfolder)
        multi_tempfile2_labeled = make_tempfile(replicate,"multi_labeled_2", tempfolder)

        label_files(unique_tempfile1, unique_tempfile1_labeled, "uniq",debug)
        label_files(unique_tempfile2, unique_tempfile2_labeled, "uniq",debug)
        label_files(multi_tempfile1, multi_tempfile1_labeled, "multi",debug)
        label_files(multi_tempfile2, multi_tempfile2_labeled, "multi",debug)

        paired_tempfile1_ulabeled = make_tempfile(replicate,"paired_ulabeled_1", tempfolder)
        paired_tempfile2_ulabeled = make_tempfile(replicate,"paired_ulabeled_2", tempfolder)
        combine_files(unique_tempfile1_labeled,multi_tempfile1_labeled, paired_tempfile1_ulabeled,debug)
        combine_files(unique_tempfile2_labeled,multi_tempfile2_labeled, paired_tempfile2_ulabeled,debug)

        #combine pairs
        if verbosity:
            print("Matching paired-end mates and merging coordinates "+ str(datetime.now()) ,file = sys.stderr)
        paired_unmatched1= make_tempfile(replicate,"paired_unmatched_1", tempfolder)
        paired_unmatched2 = make_tempfile(replicate,"paired_unmatched_2", tempfolder)
        paired_matched_tempfile = make_tempfile(replicate,"paired_matched", tempfolder)
        match_reads_python(paired_tempfile1_ulabeled,paired_tempfile2_ulabeled,strandedness,paired_matched_tempfile,paired_unmatched1, paired_unmatched2)
        # match_reads(paired_tempfile1_ulabeled,paired_tempfile2_ulabeled,strandedness,paired_matched_tempfile,paired_unmatched1, paired_unmatched2,debug) #match pairs between paired files

        #sort matched
        matched_tempfile_sorted = make_tempfile(replicate,"paired_matched_sorted", tempfolder)
        sort_temp(paired_matched_tempfile,4,matched_tempfile_sorted,debug)

        #combine start and stop of paired reads
        matched_bed = make_tempfile(replicate,"matched_bed", tempfolder)
        merge_coords(matched_tempfile_sorted,matched_bed,debug)

        # os.unlink(matched_tempfile_sorted)

        combined_unmatched = make_tempfile(replicate,"combined_unmatched", tempfolder)
        combine_files(paired_unmatched1, paired_unmatched2, combined_unmatched,debug)

        #Find single reads that are matched outside of TE but are still proper pair
        if verbosity:
            print("Adding properly paired reads that have mates outside of TE into matched file"+ str(datetime.now()) ,file = sys.stderr)
        proper_single = make_tempfile(replicate,"proper_single",tempfolder)
        combined_unmatched2 = make_tempfile(replicate,"combined_unmatched2", tempfolder)

        find_proper(combined_unmatched,combined_unmatched2,proper_single,debug) #outputs only nonproper pairs in combined_unmatched2, and single reads that are part of proper pairs in proper_single

        combined_matched = make_tempfile(replicate,"combined_matched", tempfolder)
        combine_files(matched_bed,proper_single,combined_matched,debug)
        # os.unlink(proper_single)
        ###Remove single alignments of reads that have paired matches using other valid alignments
        if verbosity:
            print("Removing single-end reads that have matching paired-end mates at other alignment locations"+ str(datetime.now()) ,file = sys.stderr)
        only_unmatched = make_tempfile(replicate,"only_unmatched", tempfolder)
        remove_repeat_reads(combined_matched,combined_unmatched2,only_unmatched,debug)

        #combine matched and unmatched alignments
        combined_bed = make_tempfile(replicate,"combined_bed", tempfolder)
        combine_files(combined_matched,only_unmatched,combined_bed,debug)

        if verbosity:
            print("Identifying and labeling unique and multi fragments"+ str(datetime.now()) ,file = sys.stderr)
        first_tempfile = make_tempfile(replicate,"first", tempfolder)
        paired_uniq_tempfile = make_tempfile(replicate,"paired_uniq_tempfile",tempfolder)
        unique_bed = make_tempfile(replicate,"unique_bed", tempfolder)
        multi_bed_pre = make_tempfile(replicate,"multi_bed_pre", tempfolder)
        find_uniq(combined_bed,first_tempfile,unique_bed,multi_bed_pre,debug)

        #find unique pairs in multi_bed
        multi_bed = make_tempfile(replicate,"multi_bed", tempfolder)
        if verbosity:
            print("Identifying multi read pairs with one end unique"+ str(datetime.now()) ,file = sys.stderr)
        find_paired_uniq(multi_bed_pre,paired_uniq_tempfile,multi_bed,unique_bed,debug)

        aligned_libsize = getlibsize(logfile, bamfile,multi_bed,unique_bed,paired_end,debug)

    ######## COUNT READ(S) #########################
    read_multidict={}  #dictionary to store TE_IDs for each read alignment
    read_locdict = {}  #dictionary to store genomic location of each alignment

    if verbosity:
        print("counting unique alignments "+ str(datetime.now()) ,file = sys.stderr)

    unique_bedfile = open(unique_bed,'r')
    avg_fraglength=uniquecount(unique_bedfile,RepCalc_dict,read_locdict)

    if verbosity:
        print("counting multi alignments "+ str(datetime.now()),file = sys.stderr)

    multi_bedfile = open(multi_bed, 'r')
    multicount(multi_bedfile,RepCalc_dict,read_multidict,read_locdict)

    unique_bedfile.close()
    multi_bedfile.close()

    if verbosity:
        print("Adding Tag information to aligned TEs "+ str(datetime.now()) ,file = sys.stderr)

    for TE_ID,RepClass in RepCalc_dict.items():
        RepClass.calcuniqRep()

    if verbosity:
                print("Calculating multialignment assignments "+ str(datetime.now()) ,file = sys.stderr)

    comparedict(read_multidict,RepCalc_dict)
    iteration=0
    if EM == "auto":
            notconverged=True
            prev_read_change=1
            prev_count_change = 0
            max_count_change = 0
            while notconverged:
                iteration +=1
                changed_count = 0
                total_TE =0
                total_TE_0 = 0
                total_TE_1 = 0
                total_TE_10 = 0
                avg_changed_count_pct =0
                max_count_change=0
                total_TE_10_1pct =0
                if verbosity:
                    print("Running expectation-maximization calculation for iteration:" + str(iteration) + " " + str(datetime.now()) ,file = sys.stderr)

                for TE_ID,RepClass in RepCalc_dict.items():
                    TE_changecount = RepClass.calcmultiRep(iteration)
                    max_count_change = max(TE_changecount,max_count_change)
                    changed_count +=TE_changecount
                    total_TE +=1
                    if TE_changecount > 0:
                        total_TE_0 +=1
                    if TE_changecount >= 1:
                        total_TE_1 += 1
                    if TE_changecount >= 1 and RepClass.counts_tot >= 10:
                        total_TE_10 += 1
                    if TE_changecount >= 1 and RepClass.counts_tot >= 10 and (TE_changecount/RepClass.counts_tot) > 0.01:
                        total_TE_10_1pct += 1
                    avg_changed_count_pct = changed_count/total_TE

                if verbosity:
                    print("Average change in TE count:" + str(avg_changed_count_pct) + " " + str(datetime.now()) ,file = sys.stderr)
                    print("Max change in TE count:" + str(max_count_change) + " " + str(datetime.now()) ,file = sys.stderr)
                    print("Number changed TE:" + str(total_TE_0) + " " + str(datetime.now()) ,file = sys.stderr)
                    print("Number TEs changed by at least 1 count:" + str(total_TE_1) + " " + str(datetime.now()) ,file = sys.stderr)
                    print("Number TEs changed by at least 1 count with at least 10 counts:" + str(total_TE_10) + " " + str(datetime.now()) ,file = sys.stderr)
                    print("Number TEs changed by at least 1 count with at least 10 counts and > 1pct total count:" + str(total_TE_10_1pct) + " " + str(datetime.now()) ,file = sys.stderr)
                new_read_change = estdict(read_multidict,RepCalc_dict)
                if total_TE_10_1pct == 0 and iteration > 1:
                    notconverged = False
                else:
                    prev_read_change = new_read_change
            if verbosity:
                print("Finished running expectation-maximization calculation after iteration:" + str(iteration) + " " + str(datetime.now()) ,file = sys.stderr)

    elif int(EM) > 0:
        notconverged=True
        prev_read_change=1
        prev_count_change = 0
        max_count_change = 0
        while iteration < int(EM):
            iteration +=1
            changed_count = 0
            total_TE =0
            total_TE_0 = 0
            total_TE_1 = 0
            total_TE_10 = 0
            avg_changed_count_pct =0
            max_count_change=0
            total_TE_10_1pct =0
            if verbosity:
                print("Running expectation-maximization calculation for iteration:" + str(iteration) + " " + str(datetime.now()) ,file = sys.stderr)

            for TE_ID,RepClass in RepCalc_dict.items():
                TE_changecount = RepClass.calcmultiRep(iteration)
                max_count_change = max(TE_changecount,max_count_change)
                changed_count +=TE_changecount
                total_TE +=1
                if TE_changecount > 0:
                    total_TE_0 +=1
                if TE_changecount >= 1:
                    total_TE_1 += 1
                if TE_changecount >= 1 and RepClass.counts_tot >= 10:
                    total_TE_10 += 1
                avg_changed_count_pct = changed_count/total_TE
                if TE_changecount >= 1 and RepClass.counts_tot >= 10 and (TE_changecount/RepClass.counts_tot) > 0.01:
                    total_TE_10_1pct += 1
            if verbosity:
                print("Average change in TE count:" + str(avg_changed_count_pct) + " " + str(datetime.now()) ,file = sys.stderr)
                print("Max change in TE count:" + str(max_count_change) + " " + str(datetime.now()) ,file = sys.stderr)
                print("Number changed TE:" + str(total_TE_0) + " " + str(datetime.now()) ,file = sys.stderr)
                print("Number TEs changed by at least 1 count:" + str(total_TE_1) + " " + str(datetime.now()) ,file = sys.stderr)
                print("Number TEs changed by at least 1 count with at least 10 counts:" + str(total_TE_10) + " " + str(datetime.now()) ,file = sys.stderr)
                print("Number TEs changed by at least 1 count with at least 10 counts and > 1pct total count:" + str(total_TE_10_1pct) + " " + str(datetime.now()) ,file = sys.stderr)
            new_read_change = estdict(read_multidict,RepCalc_dict)

    if verbosity:
                print("Writing counts "+ str(datetime.now()) ,file = sys.stderr)

    read_multidict.clear()


#     if copies:
#         temp_subF = tempfile.NamedTemporaryFile(delete=False, dir = tempfolder, prefix="count" +  ".SFtmp")
#         subF_filepath = count_dir + "/" + replicate + "_subFcounts.txt"
#         subF_file_header = open(subF_filepath + ".header",'w')
# #		subF_file_header.writelines("Sample" + "\t" + "aligned_libsize" + "\t" + "Subfamily:Family:Class" + "\t" + "copies" + "\t" + "exp_copies" + "\t" + "uniq_counts" + "\t" + "tot_counts" + "\t" + "avg_conf"  + "\t" + "tot_sense" + "\t" + "tot_antisense" + "\n")
#         subF_file_header.writelines("Sample" + "\t" + "aligned_libsize" + "\t" + "Subfamily:Family:Class" + "\t" + "copies"  + "\t" + "fpkm" + "\t" + "uniq_counts" + "\t" + "tot_counts" + "\t" + "tot_reads" + "\t" + "score"  + "\n")

#         subF_file_header.close()
#         subF_dict = {}

    counts_temp = open(counts_temp.name, 'w')
    for TE_ID,RepClass in RepCalc_dict.items(): #for each TE_ID
        RepClass.writeRep(aligned_libsize, counts_temp, replicate,strandedness,iteration)
    ##########Sort by highest total counts before writing and add to subF dictionary

        # if copies:
        #     subF = get_subF(TE_ID)
        #     subF_list = split_subF(subF)
        #     if subF not in subF_dict:
        #         subF_dict[subF]=subfamily(subF,subF_reads[subF])
        #         subF_dict[subF].add_TE_count(RepClass,strandedness)
        #     else:
        #         subF_dict[subF].add_TE_count(RepClass,strandedness)

    #Close dictionaries from memory
    counts_temp.close()

    sort_counts(counts_temp.name,counts_file_header.name,countsfilepath,5,debug) #sort on 5th column (fpkm)
    read_locdict.clear()
    RepCalc_dict.clear()

    # if copies:
    #     if verbosity:
    #             print("Writing subfamily counts "+ str(datetime.now()) ,file = sys.stderr)
    #     temp_subF = open(temp_subF.name, 'w')
    #     with open(copies,'r') as copiesfile: #copiesfile is sorted
    #         copiesfile.readline() #skip header
    #         for line in copiesfile:
    #             line = line.rstrip()
    #             line_tabs = line.split("\t")
    #             line_subF = line_tabs[0]
    #             if line_subF in subF_dict:
    #                 subF_dict[line_subF].add_copy_info(line_tabs)
    #             #### Write lines ###
    #             # temp_subF.writelines(replicate + "\t" + str(aligned_libsize) + "\t" + line_subF + "\t" + line_copies + "\t" + str(uniq) + "\t" + "{0:.2f}".format(multi) + "\t" + str(subF_conf) + "\t" + "{0:.2f}".format(sense_reads) + "\t"  +  "{0:.2f}".format(antisense_reads) + "\n")
    #                 subF_dict[line_subF].write_subfamily(temp_subF,replicate,aligned_libsize,iteration)
    #             else:
    #                 subF_dict[line_subF]=subfamily(line_subF,0)
    #                 subF_dict[line_subF].add_copy_info(line_tabs)
    #                 subF_dict[line_subF].write_subfamily(temp_subF,replicate,aligned_libsize,iteration)

    #     copiesfile.close()
    #     temp_subF.close()

    #     sort_counts(temp_subF.name, subF_file_header.name, subF_filepath, 6, debug) #Sort by 7th field (multi)

    if not debug:
        os.unlink(unique_bed)
        os.unlink(multi_bed)

    ####### STOP TIMING SCRIPT #######################
    if verbosity:
        print("finished writing outputs at "+ str(datetime.now()),file = sys.stderr)

    # candidate_dir = os.path.join(args.qual, 'candidate', sample)
    # exonTEfile = find_file(candidate_dir, "exon_TE.bed", False, 1, True)
    # exonTE_output = countsfilepath.replace("_TEcounts.txt","_exonTE.txt")
    # generate_TEexon(countsfilepath, exonTEfile, exonTE_output)
    # print("finished exonTE writing outputs at "+ str(datetime.now()),file = sys.stderr)
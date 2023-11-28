#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to subset the SAM file to include only reads (and header) that have pairs
within 100 bp of each other.

Created on Wed Oct 25 16:06:25 2023

@author: kevinmyers
"""

#open SAM file to process
with open("/Users/kevinmyers/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/CtOE-EC-33_combined_out_mapped_filtered_lengthFiltered.sam", 'r') as f:
    out_list = []
    list_final = []
    for line in f:
        if line.startswith("@"):
            list_final.append(line)
        else:
            read_name = line.split("\t")[0]
            flag = line.split("\t")[1]
            rname = line.split("\t")[2]
            first_read_start_pos = line.split("\t")[3]
            mapq = line.split("\t")[4]
            cigar = line.split("\t")[5]
            rnext = line.split("\t")[6]
            second_read_start_pos = line.split("\t")[7]
            int_dist = line.split("\t")[8]
            seq = line.split('\t')[9]
            qual = line.rstrip('\n').split('\t')[10]
            out_list.append(f"{read_name}\t{flag}\t{rname}\t{first_read_start_pos}\t{mapq}\t{cigar}\t{rnext}\t{second_read_start_pos}\t{int_dist}\t{seq}\t{qual}\n")
    with open("/Users/kevinmyers/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/CtOE-EC-33_combined_less_than_100.txt",'r') as g:
        for line in g:
            for each in out_list:
                if line.split('\t')[0] in each:
                    list_final.append(each)
                    
    with open("/Users/kevinmyers/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/CtOE-EC-33_combined_100bp_filter.sam", 'w') as r:
        for each in list_final:
            r.write(each)
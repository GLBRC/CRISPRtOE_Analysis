#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_table_for_histogram.py

Description
-----------

This is a simple script to convert BED files to tables ready to go for histogram analysis
with Histogram_CRISPRtOE.R

This script is hard coded but can be updated for future use

Created on Thu Nov  2 09:29:12 2023

@author: kevinmyers
"""

import re
import sys
import argparse

"""
These steps are only needed if using genome wide Tn-seq data as a control

#Filter the distances, for this between 0 and 500 - ONLY NEEDED IF USING GENOME WIDE Tn-seq DATA AS CONTROL
out_table = []
with open("whole_genome_tnseq_sort_dist.bed", 'r') as f:
    for line in f:
        value = line.rstrip().split("|")[2]
        if int(value) > 0 and int(value) < 500:
            out_table.append(f"{value}\n")
with open("whole_genome_tnseq_sort_dist.txt", 'w') as r:
    for line in out_table:
        r.write(line)

#Determine the distance between transposon insertion site and gene start if upstream of gene - ONLY NEEDED IF USING GENOME WIDE Tn-seq DATA AS CONTROL
whole_genome_table = []
with open("/Users/kevinmyers/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/distance_measurements/whole_genome_tnseq_sort_in_gene_dist.bed.txt", 'r') as f:
    for line in f:
        strand = line.split("\t")[8]
        value = line.rstrip().split("\t")[9]
        if strand == "+" and int(value) > 0:
            dist_value = 100 - int(value)
            whole_genome_table.append(f"{dist_value}\n")
        if strand == "-" and int(value) < 0:
            dist_value = 101 - abs(int(value))
            whole_genome_table.append(f"{dist_value}\n")
        else:
            continue
with open("/Users/kevinmyers/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/distance_measurements/whole_genome_tnseq_dist.txt", "w") as r:
    for line in whole_genome_table:
        r.write(line)
"""
def bed_process(bedfile)
#Determine the distance between transposon insertion site and gene start if inside genes
    inside_table = []
    output_file = bedfile.replace(".bed","_forHistogram.txt")
    with open(bedfile, 'r') as f:
        for line in f:
            strand = line.split("\t")[8]
            value = line.rstrip().split("\t")[9]
            if strand == "+" and int(value) > 0:
                dist_value = 100 - int(value)
                inside_table.append(f"{dist_value}\n")
            if strand == "-" and int(value) < 0:
                dist_value = 101 - abs(int(value))
                inside_table.append(f"{dist_value}\n")
            else:
                continue
    with open(output_file, "w") as r:
        for line in inside_table:
            r.write(line)

def main():
   
    cmdparser = argparse.ArgumentParser(description="Converting a BED file into a table ready for a histogram in R.",
                                        usage='%(prog)s -f <bed file to process> [optional arguments: -d]', prog='make_table_for_histogram.py'  )                                  
    cmdparser.add_argument('-f', '--file',  action='store', dest='FILE' , help='bedfile to process', metavar='')
    cmdparser.add_argument('-d', '--detail',  action='store_true', dest='DETAIL' , help='Print a more detailed description of the program.')
    cmdResults = vars(cmdparser.parse_args())
        
    BAM_files = []
    
    if cmdResults['DETAIL']:
        print("\nScript Name: make_table_for_histogram.py")
        print("\tThis script was written my Kevin Myers (kmyers2@wisc.edu).")
        print("\nPurpose : Convert a BED file to a text file that is ready for processing in R to make a histogram")
        print("\nInput : BED file to process.")
        sys.exit(1)

    if cmdResults['FILE'] is not None:
        bedfile = cmdResults['FILE']
    else:
        print("Please provide a BED file to process.\n")
        cmdparser.print_help()
        sys.exit(1)
            
    bed_process(bedfile)    
        
if __name__ == "__main__":
    main()
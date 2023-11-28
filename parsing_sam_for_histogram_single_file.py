#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parsing_sam_for_histogram_single_file.py

Description
-----------

Simple script to parse a SAM file for processing with histograms

Analyzing a SAM file to identify the read length and the distance between protospacer location and transposon insertion site

Requires Python3 and argparse, os, re, sys

Created on Fri Oct 20 13:25:33 2023

@author: kevinmyers
"""

import argparse
import os
import re
import sys

def process(each):
    """
    Open and process the SAM file. Will write mulitple files
    """
    out_text = each.replace(".sam", "_position_list.txt")
    with open(each, 'r') as f:
        spacer_length_dict = {}
        sam_out_dict = {}
        out_sam = each.replace('.sam', '_length_filtered.sam')
        for _ in range(5):
            next(f)
        for line in f:
            read_name = line.split("\t")[0]
            flag = line.split("\t")[1]
            rname = line.split("\t")[2]
            first_read_start_pos = line.split("\t")[3]
            mapq = line.split("\t")[4]
            cigar = line.split("\t")[5]
            rnext = line.split("\t")[6]
            second_read_start_pos = line.split("\t")[7]
            int_dist = line.split("\t")[8]
            read_length = len(line.split("\t")[9])
            seq = line.split('\t')[9]
            qual = line.rstrip('\n').split('\t')[10]
            int_dist2 = int(int_dist)
            read_length2 = int(read_length)
            #need to separate if value is positive or negative to determine how to calculate the distance
            if int(int_dist) > 0:
                distance_to_use = int(int_dist) - 32
            else:
                distance_to_use = abs(int(int_dist)) - 32
            
            #add the values to a dictionary
            if read_name in spacer_length_dict:
                spacer_length_dict[read_name].append(f"{distance_to_use}\t{read_length2}\t{int_dist2}")
            else:
                spacer_length_dict[read_name]=[]
                spacer_length_dict[read_name].append(f"{distance_to_use}\t{read_length2}\t{int_dist2}")
        
        #add dictionary values to a list for further processing and ease of export
        spacer_list = []
        for key,val in spacer_length_dict.items():
            spacer_list.append(f"{key}\t{val}\n")
            
        spacer_list2 = []
        for line in spacer_list:
            line2 = re.sub("\]", "", line)
            line2 = re.sub("\[", "", line2)
            line2 = re.sub("['']", "", line2)
            line2 = re.sub(", ", "\t", line2)
            line2 = line2.replace('\\t', '\t')
            line2 = line2.replace('\t\n', '\n')
            spacer_list2.append(line2)    
            
        distance_values = []
        for line in spacer_list2:
            read_name = line.split("\t")[0]
            distance1 = line.split("\t")[1]
            #distance2 = line.split("\t")[4]
            if int(distance1) > 0:
                distance_values.append(f"{read_name}\t{distance1}\n")
            else:
                continue
    
        with open(out_text, "w") as f:
            f.write("read_name\tdistance_between_reads\n")
            for each in distance_values:
                f.write(each)
        
def main():
   
    cmdparser = argparse.ArgumentParser(description="Processing SAM files for CRISPRtOE.",
                                        usage='%(prog)s -f <SAM file to process> [optional arguments: -d]', prog='parsing_sam_for_histogram.py'  )                                  
    cmdparser.add_argument('-f', '--file',  action='store', dest='FILE' , help='SAM file to process.', metavar='')
    cmdparser.add_argument('-d', '--detail',  action='store_true', dest='DETAIL' , help='Print a more detailed description of the program.')
    cmdResults = vars(cmdparser.parse_args())
    
    cwd   = os.getcwd()
    
    sam_files = []

    if cmdResults['DETAIL']:
        print("\nScript Name: parsing_sam_for_histogram.py")
        print("\tThis script was written my Kevin Myers (kmyers2@wisc.edu).")
        print("\nPurpose : Parse SAM file and determine distance from protospacer location and transposon insertion site.")
        print("\nInput : SAM file for processing.")
        sys.exit(1)
    
    if cmdResults['FILE'] is not None:
        each = cmdResults['FILE']
    else:
        print("Please provide a file with SAM files, one per line.\n")
        cmdparser.print_help()
        sys.exit(1)
    
    #with open(SAM_files, 'r') as f:
    #    for line in f:
    #        line2 = line.rstrip('\n')
    #        sam_files.append(line2)
            
    print(each)
    process( each )
        
if __name__ == "__main__":
    main()
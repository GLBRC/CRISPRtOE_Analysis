#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to parse BAM file and bedfile counts to get site of transposon insertion
and count from bedfile counts. Output is two files:
    one is a table with strand information
    the other is a WIG file for visualization
    
input is a text file of BAM files and BED files, separated by a tab, one per line

will parse BAM and BED files, find set and parse to make
table out put with strand information and counts for each position and a wig file.

requires Python3 and argparse, os, pysam, re, sys

Created on Mon Oct 30 11:28:32 2023

@author: kevinmyers
"""
import argparse
import os
import pysam
import re
import sys

def parse( bam, bedFile ):
    """
    Read in BAM file use bedFile to determine counts
    """
    data = pysam.AlignmentFile(bam, "r")
    
    data_list = []
    for line in data:
        seq = line.seq
        pos = line.pos
        tlen = line.tlen
        if tlen < 0:
            strand = "forward"
            position = pos + len(seq)
            data_list.append(f"{position}\t{strand}\n")
        else:
            strand = "reverse"
            data_list.append(f"{pos}\t{strand}\n")
                
    out_set = set(data_list)
    out_list = list(out_set)
    
    final_count_list = []
    wig_out = []
    final_dict1 = {}
    final_dict2 = {}
    #open BedFile and compare location and count data
    #write WIG and table outputs
    with open(bedFile, "r") as t:
        for line in t:
            pos = line.split("\t")[1]
            count = line.rstrip().split("\t")[3]
            final_dict1[int(pos)] = count
        for each in out_list:
            position = each.split("\t")[0]
            strand = each.rstrip("\n").split("\t")[1]
            final_dict2[int(position)] = strand

        #construct sets from each dictionary keys and find intersection    
        pileup_key_set = set(final_dict1.keys())
        samfile_key_set = set(final_dict2.keys())
        same_pos = pileup_key_set.intersection(samfile_key_set)
        
        for each in sorted(same_pos):
            if each in final_dict1 and each in final_dict2:
                final_count_list.append(f"{each}\t{final_dict1[each]}\t{final_dict2[each]}\n")
                wig_out.append(f"{each}\t{final_dict1[each]}\n")
    
    out_table = bam.replace(".bam", "_insertion_site_table.txt")
    with open(out_table, 'w') as o:
        for each in final_count_list:
            o.write(each)
    
    out_wig = bam.replace(".bam", "_insertion_site.wig")            
    with open(out_wig, 'w') as o:
        o.write("track type=wiggle_0\nvariableStep chrom=NC_000913.3 span=1\n")
        for each in wig_out:
            o.write(each)

def main():
   
    cmdparser = argparse.ArgumentParser(description="Parsing BAM and mpileup for transposon insertion locations.",
                                        usage='%(prog)s -f <list of BAM and BED files to process> [optional arguments: -d]', prog='site_of_insertion.py'  )                                  
    cmdparser.add_argument('-f', '--file',  action='store', dest='FILE' , help='list of BAM files to process, one per line.', metavar='')
    cmdparser.add_argument('-d', '--detail',  action='store_true', dest='DETAIL' , help='Print a more detailed description of the program.')
    cmdResults = vars(cmdparser.parse_args())
    
    cwd   = os.getcwd()
    
    BAM_files = []

    if cmdResults['DETAIL']:
        print("\nScript Name: site_of_insertion.py")
        print("\tThis script was written my Kevin Myers (kmyers2@wisc.edu).")
        print("\nPurpose : Script to parse BAM file and bedfile counts to get site of transposon insertion and count from bedfile counts.")
        print("\nInput : Text file with list of BAM files and associated BED file, separated by a tab, one per line.")
        sys.exit(1)
    
    if cmdResults['FILE'] is not None:
        BAM_files = cmdResults['FILE']
    else:
        print("Please provide a file with BAM files, one per line.\n")
        cmdparser.print_help()
        sys.exit(1)
    
    with open(BAM_files, 'r') as f:
        for line in f:
            bam = line.rstrip().split("\t")[0]
            bedFile = line.rstrip().split("\t")[1]
            #bedFile = bam.replace("lessThan25.bam", "_BinSize1_CPM_norm.bed")
            
            print(f"Analyzing {bam} and {bedFile} files now\n")
            
            parse(bam, bedFile)
            
        
if __name__ == "__main__":
    main()
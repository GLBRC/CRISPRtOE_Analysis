#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
counting_spacers.py

Description
-----------

This is a script to count instances of protospacers in filtered and cleaned
BAM files for CRISPRtOE experiments.

The user must supply the protospacer sequences and IDs along with a text file
listing the BAM files to analyze, one per line.

The script will count the spacers for each BAM file, write an individual results
file (*_spacer_count.txt), and then combine the results for all BAM files examined.

Notes
-----

    Input:
        
        The protospacer file must be set up in two columns as follows:

            Target_ID        Target_to_Search
            Spacer1          ATCGATCGATCGATCGATCGATCG
            Spacer2          ATCGATCGATCGATCGATCGATCG
            Spacer3          CGATCGATCGATCGATCGATCGAT
            
        The text file listing the BAM files should be formated with one BAM file per
        line, no header.
    
    Output:
        
        Individual spacer_count.txt files for each BAM tested.
        Combined spacer counts for all BAMs for downstream analysis.
        
    Dependencies:
        
        - Python3
        - Python modules argparse, defaultdict, glob, pysam, re, sys        

Created on Fri Oct 27 13:19:40 2023

@author: kevinmyers
"""
import argparse
from collections import defaultdict
import glob
import pysam
import re
import sys
        
gene_list = []
header = []
merge = defaultdict(dict)

def count_spacers(bam, spacer_table, out_file):
    """

    Parameters
    ----------
    bam : BAM file
        cleaned BAM file for analysis.
    spacer_table : TEXT file
        Two column file with protospacer ID in first column and sequence in second.
    out_file : TEXT file
        spacer counts for the BAM file examined.

    """
    def matchGene(seq, search_terms):
        """
        Generate the search_terms dictionary

        """
        if seq in search_terms:
            return search_terms[seq]
        
    search_terms = {}
    geneList = []
    
    # Genearate a list of all spacaerIDs
    with open(spacer_table, 'r') as f:
        for _ in range(1):
            next(f)
        for line in f:
            (k, v) = line.rstrip().split('\t')
            search_terms[v]=k
            if k not in geneList:
                geneList.append(k)
                
    geneListSet = set(geneList)
                
    results = {}
    
    # Open BAM and perform search and count everytime a match is found
    ec33 = pysam.AlignmentFile(bam, "r")
    for line in ec33:
        seq = line.seq
        geneName = matchGene(seq, search_terms)
        if isinstance(geneName, list):
            geneName = ','.join(geneName)
        if geneName:
            if geneName in results:
                results[geneName] += 1
            else:
                results[geneName] = 1
                
    resultsSet = set(results.keys())
    
    not_in_results = geneListSet.difference(resultsSet)
    
    # Report 0 if protospacer sequence is not found
    for each in not_in_results:
        results[each] = 0    
    
    results_sorted = {}
    
    results_sorted = dict(sorted(results.items()))
    
    # Write output file for each BAM
    with open(out_file, 'w') as f:
        for key,value in results_sorted.items():
            f.write(f"{key}\t{value}\n")
            
def getGeneList(spacer_table):
    """
    Get a list of spacer IDs
    
    """
    with open(spacer_table) as table:
        for _ in range(1):
            next(table)
        for line in table:
            spacerID = line.split("\t")[0]
            gene_list.append(spacerID)
                
def loadFiles(files):
    """

    Get file names for headers

    """
    for f in files:
        fileName = re.sub('_spacer_count.txt','',f)
        header.append(fileName)
        with open(f, 'r') as data:
            for l in data:
                l = l.rstrip()
                ct = l.split("\t")
                merge[f][ct[0]] = ct[1]
        

def main():
   
    cmdparser = argparse.ArgumentParser(description="Counting CRISPRtOE spacer instances in BAM files.",
                                        usage='%(prog)s -f <list of BAM files to process> -t <table of spacers and IDs> [optional arguments: -d]', prog='counting_spacers.py'  )                                  
    cmdparser.add_argument('-f', '--file',  action='store', dest='FILE' , help='list of BAM files to process, one per line.', metavar='')
    cmdparser.add_argument('-t', '--table',  action='store', dest='TABLE' , help='Table of IDs and spacer sequences to search, tab delimited.', metavar='')    
    cmdparser.add_argument('-d', '--detail',  action='store_true', dest='DETAIL' , help='Print a more detailed description of the program.')
    cmdResults = vars(cmdparser.parse_args())
        
    BAM_files = []
    
    if cmdResults['DETAIL']:
        print("\nScript Name: counting_spacers.py")
        print("\tThis script was written my Kevin Myers (kmyers2@wisc.edu).")
        print("\nPurpose : Count matches to protospacer sequences in a list of BAM files. Then combine the results into a single results file.")
        print("\nInput : File of protospacer IDs and sequences (tab deliminted with header) and text file list of BAM files to process (no header).")
        sys.exit(1)

    if cmdResults['FILE'] is not None:
        BAM_files = cmdResults['FILE']
    else:
        print("Please provide a file with BAM files, one per line.\n")
        cmdparser.print_help()
        sys.exit(1)
        
    if cmdResults['TABLE'] is not None:
        spacer_table = cmdResults['TABLE']
    else:
        print("Please provide a file IDs and spacer sequences, tab deliminted, one per line.\n")
        cmdparser.print_help()
        sys.exit(1) 
    
    with open(BAM_files, 'r') as f:
        for line in f:
            bam = line.rstrip('\n')
            out_file = bam.replace("_combined_out_mapped_filtered_lengthFiltered_100bpInt_filtered_sort.bam", "_spacer_count.txt")
            
            print(f"Analyzing {bam} file now...")
            
            count_spacers(bam, spacer_table, out_file)
            
    count_files = glob.glob("*_spacer_count.txt")
    
    if count_files:       
        getGeneList(spacer_table)
        loadFiles(count_files)
        
        with open('spacer_count_output_file.txt', 'w') as outFile:
            outFile.write("Spacer_ID\t")
            outHeader = "\t".join(header) + "\n"
            outFile.write(outHeader)
            
            for gene in gene_list:
                out = []
                out.append(gene)
                for file in count_files:
                    out.append(merge[file][gene])
                outLine = "\t".join(out) + "\n"
                outFile.write(outLine)
        
    
        
if __name__ == "__main__":
    main()
        

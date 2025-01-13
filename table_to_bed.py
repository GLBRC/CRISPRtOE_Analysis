#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quick script to convert the Table files from "site_of_insertion.py" script to
BED file format for comparison with genome.

Input is list of table files in an input_file.txt

Output is a BED file.

Created on Thu Nov  2 09:02:06 2023

@author: kevinmyers
"""

import argparse
import os
import sys

def process( input_file, output_file ):  
    """
    Convert file to bed file
    """ 
    pre_bed_list = []
    with open(input_file, 'r') as f:
        for line in f:
            refname = line.split("\t")[0]
            start_pos = line.split("\t")[1]
            end_pos = int(start_pos)+1
            strand = line.rstrip().split("\t")[3]
            if strand == "forward":
                pre_bed_list.append(f"{refname}\t{start_pos}\t{end_pos}\t.\t.\t+\n")
            else:
                pre_bed_list.append(f"{refname}\t{start_pos}\t{end_pos}\t.\t.\t-\n")
                
    with open(output_file, 'w') as f:
        for each in pre_bed_list:
            f.write(each)
        
def main():
   
    cmdparser = argparse.ArgumentParser(description="Converting table to BED file format.",
                                        usage='%(prog)s -f <list of table files to process> [optional arguments: -d]', prog='table_to_bed.py'  )                                  
    cmdparser.add_argument('-f', '--file',  action='store', dest='FILE' , help='list of table files to process, one per line.', metavar='')
    cmdparser.add_argument('-d', '--detail',  action='store_true', dest='DETAIL' , help='Print a more detailed description of the program.')
    cmdResults = vars(cmdparser.parse_args())
    
    cwd   = os.getcwd()
    
    table_files = []

    if cmdResults['DETAIL']:
        print("\nScript Name: table_to_bed.py")
        print("\tThis script was written my Kevin Myers (kmyers2@wisc.edu).")
        print("\nPurpose : Quick script to convert the Table files from site_of_insertion.py script to BED file format for comparison with genome.")
        print("\nInput : Text file with list of table files, one per line.")
        sys.exit(1)
    
    if cmdResults['FILE'] is not None:
        table_files = cmdResults['FILE']
    else:
        print("Please provide a file with table files, one per line.\n")
        cmdparser.print_help()
        sys.exit(1)
    
    with open(table_files, 'r') as f:
        for line in f:
            input_file = line.rstrip('\n')
            output_file = input_file.replace("_table.txt", ".bed")
            
            print(f"Analyzing {input_file} file now...")
            
            process(input_file, output_file)
            
        
if __name__ == "__main__":
    main()
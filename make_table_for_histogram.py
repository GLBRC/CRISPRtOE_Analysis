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

#Filter the distances, for this between 0 and 500
out_table = []
with open("whole_genome_tnseq_sort_dist.bed", 'r') as f:
    for line in f:
        value = line.rstrip().split("|")[2]
        if int(value) > 0 and int(value) < 500:
            out_table.append(f"{value}\n")
with open("whole_genome_tnseq_sort_dist.txt", 'w') as r:
    for line in out_table:
        r.write(line)
        
#Determine the distance between transposon insertion site and gene start if inside genes
inside_table = []
with open("/Users/kevinmyers/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/distance_measurements/combined_inside_gene.bed.txt", 'r') as f:
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
with open("/Users/kevinmyers/Documents/Peters_Lab/CRISPRtOE_Project/combined_analysis/distance_measurements/combined_inside_gene_dist.txt", "w") as r:
    for line in inside_table:
        r.write(line)

#Determine the distance between transposon insertion site and gene start if upstream of gene
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
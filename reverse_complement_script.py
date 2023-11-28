#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quick script to reverse complement a series of sequences and write the output file

Requires BioPython and PYthon3

Is hardcoded, but can be edited for other files

Created on Wed Oct 25 15:33:25 2023

@author: kevinmyers
"""

from Bio.Seq import Seq

with open("/Users/kevinmyers/Documents/Peters_Lab/CRISPRtOE_Project/rev_targets_for_RC.txt", 'r') as f:
    rev_complement_list = []
    for line in f:
        target_id = line.split("\t")[0]
        org_seq = line.rstrip("\n").split("\t")[1]
        org_seq2 = Seq(org_seq)
        revcomp_seq = org_seq2.reverse_complement()
        rev_complement_list.append(f"{target_id}\t{revcomp_seq}\n")
        
    with open("/Users/kevinmyers/Documents/Peters_Lab/CRISPRtOE_Project/rev_targets_reverse_complement.txt", 'w') as r:
        for each in rev_complement_list:
            r.write(each)
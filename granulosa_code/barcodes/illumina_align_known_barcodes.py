#!/usr/bin/env python3

import sys

import pandas as pd
import gzip
from Bio import SeqIO, pairwise2

###########################
#Setup for finding barcodes. Will use alignment
#NOTE: this ignores the read quality data that's present in the fastq
#Ignoring the read quality data speeds things up, but alignments may not be optimal.
#However, for this purpose (finding the best match) it's probably fine.

align_params = (2, -1, -1, -.5)
#match score, mismatch score, gap open penalty, gap extension penalty
#for now, these parameters are chosen arbitrarily. Need to optimize them.

align_thresh = 30 #threshold for if the alignment is rejected. Depends on align_params
#With the current settings a perfect match of a 20-mer has a score of 40.
perfect_match = 40

###########################

known_barcodes = pd.read_csv(sys.argv[2])
known_barcodes["Counts"] = 0

with gzip.open(sys.argv[1], mode = "rt") as sanger_file:
    for read in SeqIO.parse(sanger_file, "fastq"):
    
        sequence = str(read.seq)
        
        max_score = float("-Inf")
        best_index = -1
        #Check to see if each barcode matches.
        for index, barcode in enumerate(known_barcodes["Barcode"]):
            if pd.isna(barcode): continue
            score = pairwise2.align.globalms(sequence, barcode, *align_params, penalize_end_gaps = False, score_only = True)
            if score > max_score:
                max_score = score
                best_index = index
            if score >= perfect_match: break #to speed things up, exit the loop if there's a perfect match
        
        if max_score > align_thresh:
            known_barcodes.at[best_index, "Counts"] += 1
        else:
            known_barcodes.loc[known_barcodes["Plasmid"] == "Unmatched", ["Counts"]] += 1
            sys.stderr.write(sequence + '\n')
                
countfile_name = sys.argv[1].split("_merged")[0] + "_counts_aligned.csv"

print("Saving counts to "+countfile_name)

known_barcodes.to_csv(countfile_name)
sys.exit(0)
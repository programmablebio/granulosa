#!/usr/bin/env python3

import sys

import pandas as pd
import gzip
from Bio import SeqIO

###########################
#Setup for finding barcodes. Will only find EXACT matches.

known_barcodes = pd.read_csv(sys.argv[2])
known_barcodes["Counts"] = 0

with gzip.open(sys.argv[1], mode = "rt") as sanger_file:
    for read in SeqIO.parse(sanger_file, "fastq"):
    
        sequence = str(read.seq)
        
        for index, barcode in enumerate(known_barcodes["Barcode"]):
            #print(barcode)
            if barcode in sequence:
                known_barcodes.at[index, "Counts"] += 1
                continue
                
countfile_name = sys.argv[1].split("_merged")[0] + "_counts_known.csv"

print("Saving counts to "+countfile_name)

known_barcodes.to_csv(countfile_name)

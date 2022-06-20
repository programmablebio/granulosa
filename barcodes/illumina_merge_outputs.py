#!/usr/bin/env python3

import sys
import os

import pandas as pd


initialized = False #has the script read the first file?

target_folder = sys.argv[1]

merged_df = None
temp_df = None

for subfolder in os.walk(target_folder):
    for file in subfolder[2]:
        if file.endswith("counts_aligned.csv"): #find all the output files
            csv_name = os.path.join(subfolder[0],file)
            sample_name = file.split("_counts_aligned")[0]
            
            print("Reading file: " + csv_name)
            if not initialized:
                merged_df = pd.read_csv(csv_name)
                merged_df.rename(columns={"Counts": sample_name}, inplace=True)
                initialized = True
            else:
                temp_df = pd.read_csv(csv_name, usecols = ["Plasmid", "Counts"])
                temp_df.rename(columns={"Counts": sample_name}, inplace=True)
                merged_df = pd.concat([merged_df, temp_df[sample_name]], axis = 1)#, ignore_index = True)
                #merged_df = merged_df.join(temp_df.set_index('Plasmid'), on='Plasmid')
        
merged_df.to_csv(sys.argv[2], index = False)

sys.exit(0)
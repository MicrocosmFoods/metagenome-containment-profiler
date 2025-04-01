#!/usr/bin/env python3

import pandas as pd
import glob
import sys
import os

def process_tsv_files(input_dir, output_file):
    # List to store all dataframes
    dfs = []
    
    # Process each TSV file
    for tsv_file in glob.glob(os.path.join(input_dir, "*_profile.tsv")):
        # Read the TSV
        df = pd.read_csv(tsv_file, sep='\t')
        
        # Extract sample name from filename
        sample_name = os.path.basename(tsv_file).replace('_profile.tsv', '')
        
        # Add sample_name column
        df['sample_name'] = sample_name
        
        # Append to list
        dfs.append(df)
    
    # Combine all dataframes
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Filter for adjusted_ANI >= 98
    filtered_df = combined_df[combined_df['Adjusted_ANI'] >= 98]
    
    # Save to file
    filtered_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    process_tsv_files(input_dir, output_file)
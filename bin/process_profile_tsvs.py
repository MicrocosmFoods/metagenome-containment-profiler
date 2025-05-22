#!/usr/bin/env python3

import pandas as pd
import glob
import os
import argparse

def process_tsv_files(input_dir, output_file, ani_threshold=95):
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
    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)
        
        # Filter for adjusted_ANI >= threshold
        filtered_df = combined_df[combined_df['Adjusted_ANI'] >= ani_threshold]
        
        # Save to file
        filtered_df.to_csv(output_file, sep='\t', index=False)
        print(f"Processed {len(dfs)} profile files. Output saved to {output_file}")
    else:
        print("No profile files found to process.")

def parse_args():
    parser = argparse.ArgumentParser(description='Process Sylph profile TSV files and combine them.')
    parser.add_argument('input_dir', help='Directory containing profile TSV files')
    parser.add_argument('output_file', help='Path to output combined TSV file')
    parser.add_argument('--ani_threshold', type=float, default=95,
                        help='Adjusted ANI threshold for filtering (default: 95)')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    process_tsv_files(args.input_dir, args.output_file, args.ani_threshold)
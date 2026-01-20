#!/usr/bin/env python3
"""
Combine insert size features and motif features for classification
"""
import sys
import pandas as pd

def main(insert_file, motif_file, metadata_file, output_file):
    # Load insert size features
    insert_df = pd.read_csv(insert_file, sep=' ')
    
    # Load motif features
    motif_df = pd.read_csv(motif_file)
    
    # Merge on sample_id
    combined_df = pd.merge(insert_df, motif_df, on='sample_id')
    
    # Add disease labels
    metadata_df = pd.read_csv(metadata_file)
    disease_map = dict(zip(metadata_df['Run'], metadata_df['disease_status']))
    
    combined_df['disease_status'] = combined_df['sample_id'].map(disease_map)
    
    # Save
    combined_df.to_csv(output_file, index=False)
    
    print(f"Combined features created:")
    print(f"  Samples: {len(combined_df)}")
    print(f"  Total features: {len(combined_df.columns) - 2}")  # Exclude sample_id and disease_status
    print(f"  Insert size features: 8 (mean, median, stddev, min, max, q25, q75, n_fragments)")
    print(f"  Motif features: 40")
    print(f"  Output saved to: {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: combine_features.py <insert_file> <motif_file> <metadata_file> <output_file>")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

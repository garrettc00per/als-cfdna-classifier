#!/usr/bin/env python3
"""
Combine insert size, methylation, and motif features with disease labels
(Excluding positional features as they degrade performance)
"""
import sys
import pandas as pd

def main(insert_file, meth_file, motif_file, metadata_file, output_file):
    # Load features (NOTE: deliberately excluding position_file)
    insert_df = pd.read_csv(insert_file, sep=' ')
    meth_df = pd.read_csv(meth_file, sep=' ')
    motif_df = pd.read_csv(motif_file)
    
    # Merge features (no positions)
    combined_df = insert_df.merge(meth_df, on='sample_id')
    combined_df = combined_df.merge(motif_df, on='sample_id')
    
    # Add disease labels
    metadata_df = pd.read_csv(metadata_file)
    disease_map = dict(zip(metadata_df['Run'], metadata_df['disease_status']))
    combined_df['disease_status'] = combined_df['sample_id'].map(disease_map)
    
    # Save
    combined_df.to_csv(output_file, index=False)
    
    print(f"Combined features (excluding positions):")
    print(f"  Samples: {len(combined_df)}")
    print(f"  Insert size features: 8")
    print(f"  Methylation features: 8")
    print(f"  Motif features: 40")
    print(f"  Total features: {len(combined_df.columns) - 2}")
    print(f"  Output: {output_file}")

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

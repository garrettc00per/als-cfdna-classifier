#!/usr/bin/env python3
"""
Convert end motif counts to features for classification
"""
import sys
import os
import pandas as pd
import re
from collections import defaultdict

def parse_motif_file(filepath):
    """Parse a single motif file and return motif counts"""
    motifs = {'START': {}, 'END': {}}
    current_section = None
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if '=== START MOTIFS ===' in line:
                current_section = 'START'
            elif '=== END MOTIFS ===' in line:
                current_section = 'END'
            elif line and current_section:
                # Parse count and motif: "   3336 AAAA"
                parts = line.split()
                if len(parts) == 2:
                    count = int(parts[0])
                    motif = parts[1]
                    motifs[current_section][motif] = count
    
    return motifs

def main(motif_dir, output_file):
    # Get all motif files
    motif_files = [f for f in os.listdir(motif_dir) if f.endswith('_end_motifs.txt')]
    
    # Collect all unique motifs across all samples
    all_start_motifs = set()
    all_end_motifs = set()
    sample_data = {}
    
    for filename in motif_files:
        sample_id = filename.replace('_end_motifs.txt', '')
        filepath = os.path.join(motif_dir, filename)
        
        motifs = parse_motif_file(filepath)
        sample_data[sample_id] = motifs
        
        all_start_motifs.update(motifs['START'].keys())
        all_end_motifs.update(motifs['END'].keys())
    
    print(f"Found {len(sample_data)} samples")
    print(f"Found {len(all_start_motifs)} unique start motifs")
    print(f"Found {len(all_end_motifs)} unique end motifs")
    
    # Create feature matrix
    # Use top N most common motifs to avoid too many features
    # For now, let's use top 20 start and top 20 end motifs
    
    # Get top motifs by summing across all samples
    start_totals = defaultdict(int)
    end_totals = defaultdict(int)
    
    for sample_id, motifs in sample_data.items():
        for motif, count in motifs['START'].items():
            start_totals[motif] += count
        for motif, count in motifs['END'].items():
            end_totals[motif] += count
    
    # Get top 20 motifs
    top_start = sorted(start_totals.items(), key=lambda x: x[1], reverse=True)[:20]
    top_end = sorted(end_totals.items(), key=lambda x: x[1], reverse=True)[:20]
    
    top_start_motifs = [m for m, c in top_start]
    top_end_motifs = [m for m, c in top_end]
    
    print(f"\nTop 10 start motifs: {top_start_motifs[:10]}")
    print(f"Top 10 end motifs: {top_end_motifs[:10]}")
    
    # Build feature table
    rows = []
    for sample_id in sorted(sample_data.keys()):
        motifs = sample_data[sample_id]
        
        # Calculate total counts for normalization
        total_start = sum(motifs['START'].values())
        total_end = sum(motifs['END'].values())
        
        row = {'sample_id': sample_id}
        
        # Add start motif frequencies (normalized)
        for motif in top_start_motifs:
            count = motifs['START'].get(motif, 0)
            freq = count / total_start if total_start > 0 else 0
            row[f'start_{motif}'] = freq
        
        # Add end motif frequencies (normalized)
        for motif in top_end_motifs:
            count = motifs['END'].get(motif, 0)
            freq = count / total_end if total_end > 0 else 0
            row[f'end_{motif}'] = freq
        
        rows.append(row)
    
    # Create dataframe and save
    df = pd.DataFrame(rows)
    df.to_csv(output_file, index=False)
    print(f"\nMotif features saved to: {output_file}")
    print(f"Feature dimensions: {len(df)} samples x {len(df.columns)-1} features")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: create_motif_features.py <motif_dir> <output_file>")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2])

#!/usr/bin/env python3
"""
Visualize fragment start and end position distributions
"""
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os

def plot_position_distributions(position_dir, metadata_file, output_file):
    """Plot positional distributions for ALS vs Control"""
    
    # Load metadata
    metadata_df = pd.read_csv(metadata_file)
    disease_map = dict(zip(metadata_df['Run'], metadata_df['disease_status']))
    
    # Read all position bin files
    bin_files = glob.glob(os.path.join(position_dir, '*_position_bins.txt'))
    
    als_bins = {}
    ctrl_bins = {}
    
    for bin_file in bin_files:
        sample_id = os.path.basename(bin_file).replace('_position_bins.txt', '')
        disease = disease_map.get(sample_id, 'unknown')
        
        if disease == 'unknown':
            continue
        
        # Read bins (format: bin_number count)
        with open(bin_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split()
                if len(parts) == 2:
                    bin_num = int(parts[0])
                    count = int(parts[1])
                    
                    if disease == 'als':
                        als_bins[bin_num] = als_bins.get(bin_num, 0) + count
                    else:
                        ctrl_bins[bin_num] = ctrl_bins.get(bin_num, 0) + count
    
    # Create figure
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    # Plot 1: Coverage across chr21 (1Mb bins)
    ax = axes[0]
    
    # Get bin ranges
    max_bin = max(max(als_bins.keys()), max(ctrl_bins.keys()))
    bins = range(0, max_bin + 1)
    
    als_counts = [als_bins.get(b, 0) for b in bins]
    ctrl_counts = [ctrl_bins.get(b, 0) for b in bins]
    
    # Normalize by total fragments
    als_total = sum(als_counts)
    ctrl_total = sum(ctrl_counts)
    als_norm = [c / als_total * 100 for c in als_counts]
    ctrl_norm = [c / ctrl_total * 100 for c in ctrl_counts]
    
    # Plot as lines
    positions_mb = [b for b in bins]  # Already in Mb
    ax.plot(positions_mb, als_norm, color='#e74c3c', linewidth=2, label='ALS', alpha=0.8)
    ax.plot(positions_mb, ctrl_norm, color='#3498db', linewidth=2, label='Control', alpha=0.8)
    
    ax.set_xlabel('Genomic Position on chr21 (Mb)', fontsize=12)
    ax.set_ylabel('Fragment Density (%)', fontsize=12)
    ax.set_title('Fragment Start Position Distribution Across chr21', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(alpha=0.3)
    
    # Plot 2: Regional distribution (3 regions)
    ax = axes[1]
    
    # Read position summary for regional counts
    position_summary = pd.read_csv(os.path.join(position_dir, '../position_summary.txt'), sep=' ')
    position_summary['disease_status'] = position_summary['sample_id'].map(disease_map)
    
    # Calculate proportions
    position_summary['region1_prop'] = position_summary['region1_count'] / position_summary['total_fragments']
    position_summary['region2_prop'] = position_summary['region2_count'] / position_summary['total_fragments']
    position_summary['region3_prop'] = position_summary['region3_count'] / position_summary['total_fragments']
    
    als_data = position_summary[position_summary['disease_status'] == 'als']
    ctrl_data = position_summary[position_summary['disease_status'] == 'ctrl']
    
    # Bar plot
    regions = ['Region 1\n(0-16 Mb)', 'Region 2\n(16-32 Mb)', 'Region 3\n(32-48 Mb)']
    als_means = [als_data['region1_prop'].mean(), als_data['region2_prop'].mean(), als_data['region3_prop'].mean()]
    ctrl_means = [ctrl_data['region1_prop'].mean(), ctrl_data['region2_prop'].mean(), ctrl_data['region3_prop'].mean()]
    
    x = range(len(regions))
    width = 0.35
    
    ax.bar([i - width/2 for i in x], ctrl_means, width, label='Control', color='#3498db', alpha=0.8)
    ax.bar([i + width/2 for i in x], als_means, width, label='ALS', color='#e74c3c', alpha=0.8)
    
    ax.set_ylabel('Proportion of Fragments', fontsize=12)
    ax.set_xlabel('Genomic Region', fontsize=12)
    ax.set_title('Regional Distribution of Fragment Starts', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(regions)
    ax.legend(fontsize=11)
    ax.grid(alpha=0.3, axis='y')
    
    # Add p-values
    from scipy import stats
    for i, region in enumerate(['region1_prop', 'region2_prop', 'region3_prop']):
        t_stat, p_val = stats.ttest_ind(ctrl_data[region], als_data[region])
        y_pos = max(ctrl_means[i], als_means[i]) * 1.05
        ax.text(i, y_pos, f'p={p_val:.3f}', ha='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Position distribution plot saved to: {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: plot_position_distributions.py <position_dir> <metadata_file> <output_file>")
        sys.exit(1)
    
    plot_position_distributions(sys.argv[1], sys.argv[2], sys.argv[3])

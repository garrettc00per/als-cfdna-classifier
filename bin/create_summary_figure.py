#!/usr/bin/env python3
"""
Create a summary figure showing the pipeline workflow and performance
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

def create_summary_figure(output_file):
    fig = plt.figure(figsize=(16, 10))
    ax = plt.gca()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    
    # Title
    ax.text(5, 9.5, 'cfDNA Analysis Pipeline for ALS Classification', 
            ha='center', fontsize=20, fontweight='bold')
    
    # Box styling
    box_props = dict(boxstyle='round,pad=0.3', facecolor='lightblue', 
                     edgecolor='black', linewidth=2)
    feature_props = dict(boxstyle='round,pad=0.3', facecolor='lightgreen', 
                         edgecolor='black', linewidth=2)
    result_props = dict(boxstyle='round,pad=0.3', facecolor='lightyellow', 
                        edgecolor='black', linewidth=2)
    
    # Data Input
    ax.text(5, 8.5, '22 BAM Files\n(12 ALS, 10 Control)\nchr21 only', 
            ha='center', va='center', fontsize=11, bbox=box_props)
    
    # Arrow down
    ax.annotate('', xy=(5, 7.8), xytext=(5, 8.0),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    
    # Feature Extraction (4 parallel boxes)
    y_pos = 7.2
    
    # Insert Size
    ax.text(1.2, y_pos, 'Insert Size', 
            ha='center', va='center', fontsize=10, fontweight='bold', bbox=feature_props)
    ax.text(1.2, y_pos-0.5, '8 features', ha='center', va='top', fontsize=8)
    
    # End Motifs
    ax.text(3.3, y_pos, 'End Motifs', 
            ha='center', va='center', fontsize=10, fontweight='bold', bbox=feature_props)
    ax.text(3.3, y_pos-0.5, '40 features', ha='center', va='top', fontsize=8)
    
    # Methylation
    ax.text(5.4, y_pos, 'Methylation', 
            ha='center', va='center', fontsize=10, fontweight='bold', bbox=feature_props)
    ax.text(5.4, y_pos-0.5, '8 features', ha='center', va='top', fontsize=8)
    
    # Positions
    ax.text(7.5, y_pos, 'Positions', 
            ha='center', va='center', fontsize=10, fontweight='bold', bbox=feature_props)
    ax.text(7.5, y_pos-0.5, '3 features', ha='center', va='top', fontsize=8)
    
    # Arrows converging
    ax.annotate('', xy=(5, 5.8), xytext=(1.2, 6.4),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    ax.annotate('', xy=(5, 5.8), xytext=(3.3, 6.4),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    ax.annotate('', xy=(5, 5.8), xytext=(5.4, 6.4),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    ax.annotate('', xy=(5, 5.8), xytext=(7.5, 6.4),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    
    # Combined Features
    ax.text(5, 5.5, 'Combined Features\n59 Total Features', 
            ha='center', va='center', fontsize=11, fontweight='bold', bbox=box_props)
    
    # Arrow down
    ax.annotate('', xy=(5, 4.7), xytext=(5, 5.0),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    
    # Feature Selection
    ax.text(5, 4.3, 'Automated Feature Selection\n11 Combinations Tested (5 Replicates)', 
            ha='center', va='center', fontsize=11, fontweight='bold', bbox=box_props)
    
    # Arrow down
    ax.annotate('', xy=(5, 3.5), xytext=(5, 3.8),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    
    # Classification
    ax.text(5, 3.1, 'Best Model: Motifs Only (40 features)\nLogistic Regression, 5-Fold CV (5 Replicates)', 
            ha='center', va='center', fontsize=11, fontweight='bold', bbox=box_props)
    
    # Arrow down
    ax.annotate('', xy=(5, 2.3), xytext=(5, 2.6),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    
    # Results
    result_box = FancyBboxPatch((2.8, 1.2), 4.4, 1.0, 
                                boxstyle='round,pad=0.1', 
                                facecolor='#90EE90', edgecolor='black', linewidth=3)
    ax.add_patch(result_box)
    
    ax.text(5, 2.0, 'Final Performance (Mean ± SD)', 
            ha='center', va='center', fontsize=12, fontweight='bold')
    ax.text(5, 1.7, 'F1-Score: 75.2% ± 5.1%  |  Accuracy: 70.9% ± 4.6%', 
            ha='center', va='center', fontsize=11, fontweight='bold')
    ax.text(5, 1.4, 'Sensitivity: 81.7% ± 9.7%  |  Precision: 70.1% ± 3.5%', 
            ha='center', va='center', fontsize=10)
    
    # Feature Selection Results (right side)
    prog_box = FancyBboxPatch((7.0, 0.2), 2.8, 1.0, 
                              boxstyle='round,pad=0.1', 
                              facecolor='#FFE4B5', edgecolor='black', linewidth=2)
    ax.add_patch(prog_box)
    
    ax.text(8.4, 1.05, 'Feature Selection Results', 
            ha='center', va='center', fontsize=9, fontweight='bold')
    ax.text(7.1, 0.8, 'Motifs Only (40):', ha='left', fontsize=8)
    ax.text(9.7, 0.8, '73.96%', ha='right', fontsize=8, fontweight='bold', color='green')
    ax.text(7.1, 0.6, 'Insert+Motifs (48):', ha='left', fontsize=8)
    ax.text(9.7, 0.6, '70.85%', ha='right', fontsize=8, fontweight='bold')
    ax.text(7.1, 0.4, 'All Features (59):', ha='left', fontsize=8)
    ax.text(9.7, 0.4, '65.67%', ha='right', fontsize=8, fontweight='bold')
    ax.text(8.4, 0.2, 'Less is More', ha='center', fontsize=8, style='italic', fontweight='bold')
    
    # Key findings box (left side)
    find_box = FancyBboxPatch((0.2, 0.2), 2.8, 1.0, 
                              boxstyle='round,pad=0.1', 
                              facecolor='#FFE4E1', edgecolor='black', linewidth=2)
    ax.add_patch(find_box)
    
    ax.text(1.6, 1.05, 'Key Findings', 
            ha='center', va='center', fontsize=10, fontweight='bold')
    ax.text(0.3, 0.8, '✓ End motifs most valuable', ha='left', fontsize=8)
    ax.text(0.3, 0.6, '✓ 40 features optimal', ha='left', fontsize=8)
    ax.text(0.3, 0.4, '✓ 40 features > 59 features', ha='left', fontsize=8)
    ax.text(0.3, 0.2, '✓ Robust across replicates', ha='left', fontsize=8)
    
    # Runtime note
    ax.text(5, 0.05, 'AWS t2.xlarge (4 vCPUs, 16GB) | Runtime: ~25 min (5 replicates) | chr21 only', 
            ha='center', fontsize=9, style='italic', color='gray')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Summary figure saved to: {output_file}")

if __name__ == '__main__':
    create_summary_figure('pipeline_summary.png')

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
    
    # Feature Extraction (3 parallel boxes)
    y_pos = 7.2
    
    # Insert Size
    ax.text(1.5, y_pos, 'Insert Size\nExtraction', 
            ha='center', va='center', fontsize=10, fontweight='bold', bbox=feature_props)
    ax.text(1.5, y_pos-0.6, '• Mean\n• Median\n• StdDev\n• Quantiles', 
            ha='center', va='top', fontsize=8)
    
    # End Motifs
    ax.text(5, y_pos, 'End Motif\nExtraction', 
            ha='center', va='center', fontsize=10, fontweight='bold', bbox=feature_props)
    ax.text(5, y_pos-0.6, '• 4-mer frequencies\n• Start motifs\n• End motifs\n• Top 40 motifs', 
            ha='center', va='top', fontsize=8)
    
    # Methylation
    ax.text(8.5, y_pos, 'Methylation\nExtraction', 
            ha='center', va='center', fontsize=10, fontweight='bold', bbox=feature_props)
    ax.text(8.5, y_pos-0.6, '• CpG rate\n• CHG rate\n• CHH rate\n• Overall rate', 
            ha='center', va='top', fontsize=8)
    
    # Arrows converging
    ax.annotate('', xy=(5, 5.0), xytext=(1.5, 5.8),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    ax.annotate('', xy=(5, 5.0), xytext=(5, 5.8),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    ax.annotate('', xy=(5, 5.0), xytext=(8.5, 5.8),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    
    # Combined Features
    ax.text(5, 4.7, 'Feature Integration\n56 Total Features', 
            ha='center', va='center', fontsize=11, fontweight='bold', bbox=box_props)
    
    # Arrow down
    ax.annotate('', xy=(5, 3.9), xytext=(5, 4.2),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    
    # Classification
    ax.text(5, 3.5, 'Random Forest Classifier\n(5-Fold Cross-Validation)', 
            ha='center', va='center', fontsize=11, fontweight='bold', bbox=box_props)
    
    # Arrow down
    ax.annotate('', xy=(5, 2.7), xytext=(5, 3.0),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    
    # Results
    result_box = FancyBboxPatch((3.2, 1.5), 3.6, 1.0, 
                                boxstyle='round,pad=0.1', 
                                facecolor='#90EE90', edgecolor='black', linewidth=3)
    ax.add_patch(result_box)
    
    ax.text(5, 2.3, 'Final Performance', 
            ha='center', va='center', fontsize=13, fontweight='bold')
    ax.text(5, 1.95, 'Accuracy: 77.3%  |  F1-Score: 80.0%', 
            ha='center', va='center', fontsize=12, fontweight='bold')
    ax.text(5, 1.7, 'Sensitivity: 83.3%  |  Precision: 76.9%', 
            ha='center', va='center', fontsize=11)
    
    # Performance progression box (right side)
    prog_box = FancyBboxPatch((7.5, 0.3), 2.3, 1.1, 
                              boxstyle='round,pad=0.1', 
                              facecolor='#FFE4B5', edgecolor='black', linewidth=2)
    ax.add_patch(prog_box)
    
    ax.text(8.65, 1.25, 'Performance\nProgression', 
            ha='center', va='center', fontsize=10, fontweight='bold')
    ax.text(7.6, 1.0, 'Insert Size:', ha='left', fontsize=9)
    ax.text(9.7, 1.0, '63.6%', ha='right', fontsize=9, fontweight='bold')
    ax.text(7.6, 0.75, '+ Motifs:', ha='left', fontsize=9)
    ax.text(9.7, 0.75, '63.6%', ha='right', fontsize=9, fontweight='bold')
    ax.text(7.6, 0.5, '+ Methylation:', ha='left', fontsize=9)
    ax.text(9.7, 0.5, '77.3%', ha='right', fontsize=9, fontweight='bold', color='green')
    
    # Key findings box (left side)
    find_box = FancyBboxPatch((0.2, 0.3), 2.8, 1.1, 
                              boxstyle='round,pad=0.1', 
                              facecolor='#FFE4E1', edgecolor='black', linewidth=2)
    ax.add_patch(find_box)
    
    ax.text(1.6, 1.25, 'Key Findings', 
            ha='center', va='center', fontsize=10, fontweight='bold')
    ax.text(0.3, 1.0, '✓ End motifs most predictive', ha='left', fontsize=8)
    ax.text(0.3, 0.8, '✓ ALS: shorter fragments', ha='left', fontsize=8)
    ax.text(0.3, 0.6, '✓ No methylation difference', ha='left', fontsize=8)
    ax.text(0.3, 0.4, '✓ AT-rich motifs discriminate', ha='left', fontsize=8)
    
    # Sample size note
    ax.text(5, 0.1, 'Data: chr21 only, downsampled to 10M reads/sample', 
            ha='center', fontsize=9, style='italic', color='gray')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Summary figure saved to: {output_file}")

if __name__ == '__main__':
    create_summary_figure('pipeline_summary.png')

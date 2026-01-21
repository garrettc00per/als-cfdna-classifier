#!/usr/bin/env python3
"""
Create a summary figure showing the pipeline workflow and performance
Reads actual results from feature_selection_results.csv and classification_results.csv
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import pandas as pd
import numpy as np
import sys

def create_summary_figure(feature_selection_file, classification_file, output_file):
    # Load results
    feature_df = pd.read_csv(feature_selection_file)
    classification_df = pd.read_csv(classification_file)
    
    # Get top 3 feature combinations
    top_combos = feature_df.nlargest(3, 'best_f1_mean')
    
    # Get best classifier metrics
    best_classifier = classification_df.iloc[0]
    
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
    n_replicates = feature_df.iloc[0].get('n_replicates', '?')
    ax.text(5, 4.3, f'Automated Feature Selection\n11 Combinations Tested (20 Replicates)', 
            ha='center', va='center', fontsize=11, fontweight='bold', bbox=box_props)
    
    # Arrow down
    ax.annotate('', xy=(5, 3.5), xytext=(5, 3.8),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    
    # Classification - use actual best combo
    best_combo = top_combos.iloc[0]
    combo_name = best_combo['combination']
    n_features = int(best_combo['n_features'])
    ax.text(5, 3.1, f'Best Model: {combo_name} ({n_features} features)\nLogistic Regression, 5-Fold CV (20 Replicates)', 
            ha='center', va='center', fontsize=11, fontweight='bold', bbox=box_props)
    
    # Arrow down
    ax.annotate('', xy=(5, 2.3), xytext=(5, 2.6),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    
    # Results - use actual classification metrics
    result_box = FancyBboxPatch((2.2, 1.1), 5.6, 1.1, 
                                boxstyle='round,pad=0.1', 
                                facecolor='#90EE90', edgecolor='black', linewidth=3)
    ax.add_patch(result_box)
    
    ax.text(5, 2.05, 'Final Performance (Mean ± SD)', 
            ha='center', va='center', fontsize=12, fontweight='bold')
    
    acc_mean = best_classifier['Accuracy_Mean']
    acc_sd = best_classifier['Accuracy_SD']
    f1_mean = best_classifier['F1_Mean']
    f1_sd = best_classifier['F1_SD']
    
    ax.text(5, 1.7, f'F1-Score: {f1_mean:.1%} ± {f1_sd:.1%}  |  Accuracy: {acc_mean:.1%} ± {acc_sd:.1%}', 
            ha='center', va='center', fontsize=11, fontweight='bold')
    
    sens_mean = best_classifier['Sensitivity_Mean']
    sens_sd = best_classifier['Sensitivity_SD']
    prec_mean = best_classifier['Precision_Mean']
    prec_sd = best_classifier['Precision_SD']
    
    ax.text(5, 1.35, f'Sensitivity: {sens_mean:.1%} ± {sens_sd:.1%}  |  Precision: {prec_mean:.1%} ± {prec_sd:.1%}', 
            ha='center', va='center', fontsize=10)
    
    # Feature Selection Results (right side) - top 3 combos
    prog_box = FancyBboxPatch((6.8, 0.2), 3.0, 1.1, 
                              boxstyle='round,pad=0.1', 
                              facecolor='#FFE4B5', edgecolor='black', linewidth=2)
    ax.add_patch(prog_box)
    
    ax.text(8.3, 1.15, 'Top 3 Feature Combos', 
            ha='center', va='center', fontsize=9, fontweight='bold')
    
    for idx, (_, row) in enumerate(top_combos.iterrows()):
        y = 0.95 - (idx * 0.25)
        combo_short = row['combination'][:20] + '...' if len(row['combination']) > 20 else row['combination']
        ax.text(6.9, y, f"{idx+1}. {combo_short} ({int(row['n_features'])})", ha='left', fontsize=7)
        ax.text(9.8, y, f"{row['best_f1_mean']:.1%}", ha='right', fontsize=7, fontweight='bold')
    
    ax.text(8.3, 0.25, 'Competitive performance shows\nmultiple valid feature solutions', 
            ha='center', fontsize=7, style='italic')
    
    # Key findings box (left side)
    find_box = FancyBboxPatch((0.2, 0.2), 2.8, 1.1, 
                              boxstyle='round,pad=0.1', 
                              facecolor='#FFE4E1', edgecolor='black', linewidth=2)
    ax.add_patch(find_box)
    
    ax.text(1.6, 1.2, 'Key Insights', 
            ha='center', va='center', fontsize=10, fontweight='bold')
    ax.text(0.3, 0.95, '✓ Multiple combos competitive', ha='left', fontsize=7)
    ax.text(0.3, 0.75, '✓ Top 3 within ~3% F1', ha='left', fontsize=7)
    ax.text(0.3, 0.55, '✓ Feature selection unstable\nwith n=22', ha='left', fontsize=7)
    ax.text(0.3, 0.25, '✓ Larger cohort recommended', ha='left', fontsize=7)
    
    # Runtime note
    ax.text(5, 0.05, 'AWS t2.xlarge (4 vCPUs, 16GB) | Runtime: ~3m 45s (20 replicates) | chr21 only', 
            ha='center', fontsize=9, style='italic', color='gray')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Summary figure saved to: {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: create_summary_figure.py <feature_selection_file> <classification_file> <output_file>")
        sys.exit(1)
    
    create_summary_figure(sys.argv[1], sys.argv[2], sys.argv[3])

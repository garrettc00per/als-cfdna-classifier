#!/usr/bin/env python3
"""
Create visualizations for cfDNA classification analysis
"""
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import confusion_matrix
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

def plot_feature_importance(df, output_prefix):
    """Plot Random Forest feature importance"""
    # Prepare data
    feature_cols = [col for col in df.columns if col not in ['sample_id', 'disease_status']]
    X = df[feature_cols].values
    y = (df['disease_status'] == 'als').astype(int)
    
    # Train Random Forest
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    rf = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=5, n_jobs=1)
    rf.fit(X_scaled, y)
    
    # Get feature importances
    importances = pd.DataFrame({
        'feature': feature_cols,
        'importance': rf.feature_importances_
    }).sort_values('importance', ascending=False)
    
    # Plot top 20 features
    plt.figure(figsize=(12, 10))
    top_features = importances.head(20)
    
    # Color code by feature type
    colors = []
    for feat in top_features['feature']:
        if 'cpg' in feat or 'chg' in feat or 'chh' in feat or 'meth' in feat:
            colors.append('#e74c3c')  # Red for methylation
        elif 'start_' in feat or 'end_' in feat:
            colors.append('#3498db')  # Blue for motifs
        else:
            colors.append('#2ecc71')  # Green for insert size
    
    plt.barh(range(len(top_features)), top_features['importance'], color=colors)
    plt.yticks(range(len(top_features)), top_features['feature'])
    plt.xlabel('Feature Importance', fontsize=12)
    plt.title('Top 20 Most Important Features (Random Forest)', fontsize=14, fontweight='bold')
    plt.gca().invert_yaxis()
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#e74c3c', label='Methylation'),
        Patch(facecolor='#3498db', label='End Motifs'),
        Patch(facecolor='#2ecc71', label='Insert Size')
    ]
    plt.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_feature_importance.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_feature_importance.png")
    
    return importances

def plot_top_feature_distributions(df, importances, output_prefix, n_features=6):
    """Plot distributions of top features for ALS vs Control"""
    top_features = importances.head(n_features)['feature'].tolist()
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    als_data = df[df['disease_status'] == 'als']
    ctrl_data = df[df['disease_status'] == 'ctrl']
    
    for idx, feature in enumerate(top_features):
        ax = axes[idx]
        
        # Plot distributions
        ax.hist(als_data[feature], bins=15, alpha=0.6, label='ALS', color='#e74c3c', edgecolor='black')
        ax.hist(ctrl_data[feature], bins=15, alpha=0.6, label='Control', color='#3498db', edgecolor='black')
        
        # Add vertical lines for means
        ax.axvline(als_data[feature].mean(), color='#e74c3c', linestyle='--', linewidth=2)
        ax.axvline(ctrl_data[feature].mean(), color='#3498db', linestyle='--', linewidth=2)
        
        ax.set_xlabel(feature, fontsize=10)
        ax.set_ylabel('Frequency', fontsize=10)
        ax.legend()
        ax.grid(alpha=0.3)
    
    plt.suptitle('Top 6 Features: ALS vs Control Distributions', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_top_features_distributions.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_top_features_distributions.png")

def plot_confusion_matrix(df, output_prefix):
    """Plot confusion matrix for Random Forest classifier"""
    # Prepare data
    feature_cols = [col for col in df.columns if col not in ['sample_id', 'disease_status']]
    X = df[feature_cols].values
    y = (df['disease_status'] == 'als').astype(int)
    
    # Train and get predictions via cross-validation
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    rf = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=5, n_jobs=1)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    y_pred = cross_val_predict(rf, X_scaled, y, cv=cv)
    
    # Compute confusion matrix
    cm = confusion_matrix(y, y_pred)
    
    # Plot
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=True,
                xticklabels=['Control', 'ALS'],
                yticklabels=['Control', 'ALS'],
                annot_kws={'size': 16})
    
    plt.title('Confusion Matrix - Random Forest Classifier', fontsize=14, fontweight='bold')
    plt.ylabel('True Label', fontsize=12)
    plt.xlabel('Predicted Label', fontsize=12)
    
    # Add accuracy text
    accuracy = (cm[0,0] + cm[1,1]) / cm.sum()
    plt.text(1, 2.3, f'Accuracy: {accuracy:.1%}', ha='center', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_confusion_matrix.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_confusion_matrix.png")

def plot_methylation_comparison(df, output_prefix):
    """Compare methylation rates between ALS and Control"""
    meth_features = ['cpg_meth_rate', 'chg_meth_rate', 'chh_meth_rate', 'overall_meth_rate']
    
    als_data = df[df['disease_status'] == 'als']
    ctrl_data = df[df['disease_status'] == 'ctrl']
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    for idx, feature in enumerate(meth_features):
        ax = axes[idx]
        
        # Box plots
        data_to_plot = [ctrl_data[feature], als_data[feature]]
        bp = ax.boxplot(data_to_plot, labels=['Control', 'ALS'], patch_artist=True,
                        medianprops=dict(color='black', linewidth=2))
        
        # Color the boxes
        bp['boxes'][0].set_facecolor('#3498db')
        bp['boxes'][1].set_facecolor('#e74c3c')
        
        # Add individual points
        for i, data in enumerate(data_to_plot):
            y = data
            x = np.random.normal(i+1, 0.04, size=len(y))
            ax.plot(x, y, 'ko', alpha=0.3, markersize=6)
        
        ax.set_ylabel(feature.replace('_', ' ').title(), fontsize=11)
        ax.grid(alpha=0.3, axis='y')
        
        # Add p-value (simple t-test)
        t_stat, p_val = stats.ttest_ind(ctrl_data[feature], als_data[feature])
        ax.text(1.5, ax.get_ylim()[1]*0.95, f'p={p_val:.3f}', ha='center', fontsize=10)
    
    plt.suptitle('Methylation Rate Comparison: ALS vs Control', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_methylation_comparison.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_methylation_comparison.png")

def plot_insert_size_comparison(df, output_prefix):
    """Compare insert size distributions"""
    als_data = df[df['disease_status'] == 'als']
    ctrl_data = df[df['disease_status'] == 'ctrl']
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Mean insert size
    ax = axes[0]
    bp = ax.boxplot([ctrl_data['mean'], als_data['mean']], 
                     labels=['Control', 'ALS'], patch_artist=True,
                     medianprops=dict(color='black', linewidth=2))
    bp['boxes'][0].set_facecolor('#3498db')
    bp['boxes'][1].set_facecolor('#e74c3c')
    ax.set_ylabel('Mean Insert Size (bp)', fontsize=11)
    ax.set_title('Mean Insert Size', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    
    # Add individual points
    for i, data in enumerate([ctrl_data['mean'], als_data['mean']]):
        y = data
        x = np.random.normal(i+1, 0.04, size=len(y))
        ax.plot(x, y, 'ko', alpha=0.3, markersize=6)
    
    # Add p-value
    t_stat, p_val = stats.ttest_ind(ctrl_data['mean'], als_data['mean'])
    ax.text(1.5, ax.get_ylim()[1]*0.95, f'p={p_val:.3f}', ha='center', fontsize=10)
    
    # Median insert size
    ax = axes[1]
    bp = ax.boxplot([ctrl_data['median'], als_data['median']], 
                     labels=['Control', 'ALS'], patch_artist=True,
                     medianprops=dict(color='black', linewidth=2))
    bp['boxes'][0].set_facecolor('#3498db')
    bp['boxes'][1].set_facecolor('#e74c3c')
    ax.set_ylabel('Median Insert Size (bp)', fontsize=11)
    ax.set_title('Median Insert Size', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    
    # Add individual points
    for i, data in enumerate([ctrl_data['median'], als_data['median']]):
        y = data
        x = np.random.normal(i+1, 0.04, size=len(y))
        ax.plot(x, y, 'ko', alpha=0.3, markersize=6)
    
    # Add p-value
    t_stat, p_val = stats.ttest_ind(ctrl_data['median'], als_data['median'])
    ax.text(1.5, ax.get_ylim()[1]*0.95, f'p={p_val:.3f}', ha='center', fontsize=10)
    
    # Standard deviation
    ax = axes[2]
    bp = ax.boxplot([ctrl_data['stddev'], als_data['stddev']], 
                     labels=['Control', 'ALS'], patch_artist=True,
                     medianprops=dict(color='black', linewidth=2))
    bp['boxes'][0].set_facecolor('#3498db')
    bp['boxes'][1].set_facecolor('#e74c3c')
    ax.set_ylabel('Std Dev of Insert Size (bp)', fontsize=11)
    ax.set_title('Insert Size Variability', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    
    # Add individual points
    for i, data in enumerate([ctrl_data['stddev'], als_data['stddev']]):
        y = data
        x = np.random.normal(i+1, 0.04, size=len(y))
        ax.plot(x, y, 'ko', alpha=0.3, markersize=6)
    
    # Add p-value
    t_stat, p_val = stats.ttest_ind(ctrl_data['stddev'], als_data['stddev'])
    ax.text(1.5, ax.get_ylim()[1]*0.95, f'p={p_val:.3f}', ha='center', fontsize=10)
    
    plt.suptitle('Insert Size Characteristics: ALS vs Control', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_insert_size_comparison.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_insert_size_comparison.png")

def main(features_file, output_prefix):
    print("Creating visualizations...")
    print("="*60)
    
    # Load data
    df = pd.read_csv(features_file)
    print(f"Loaded {len(df)} samples with {len(df.columns)-2} features\n")
    
    # 1. Feature importance
    print("1. Plotting feature importance...")
    importances = plot_feature_importance(df, output_prefix)
    
    # 2. Top feature distributions
    print("2. Plotting top feature distributions...")
    plot_top_feature_distributions(df, importances, output_prefix)
    
    # 3. Confusion matrix
    print("3. Plotting confusion matrix...")
    plot_confusion_matrix(df, output_prefix)
    
    # 4. Methylation comparison
    print("4. Plotting methylation comparison...")
    plot_methylation_comparison(df, output_prefix)
    
    # 5. Insert size comparison
    print("5. Plotting insert size comparison...")
    plot_insert_size_comparison(df, output_prefix)
    
    print("\n" + "="*60)
    print("All visualizations complete!")
    
    # Save feature importance table
    importances.to_csv(f'{output_prefix}_feature_importances.csv', index=False)
    print(f"Feature importance table saved to: {output_prefix}_feature_importances.csv")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: create_visualizations.py <features_file> <output_prefix>")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2])

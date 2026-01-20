#!/usr/bin/env python3
"""
Test different feature combinations and select the best
"""
import sys
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold
import json

def main(features_file, output_file, best_features_file):
    # Load data
    df = pd.read_csv(features_file)
    
    # Define feature groups
    all_cols = df.columns.tolist()
    feature_groups = {
        'insert_size': ['mean', 'median', 'stddev', 'min', 'max', 'q25', 'q75', 'n_fragments'],
        'methylation': ['cpg_meth_rate', 'chg_meth_rate', 'chh_meth_rate', 'overall_meth_rate', 
                        'total_cpg', 'total_chg', 'total_chh', 'n_reads'],
        'motifs': [col for col in all_cols if 'start_' in col or 'end_' in col],
        'positions': ['region1_count', 'region2_count', 'region3_count'],
    }
    
    # Labels
    y = (df['disease_status'] == 'als').astype(int)
    
    # Test combinations
    print("=" * 70)
    print("FEATURE COMBINATION TESTING")
    print("=" * 70)
    print(f"Dataset: {len(df)} samples ({sum(y)} ALS, {len(y)-sum(y)} Control)")
    print()
    
    combinations = [
        ('Insert Size Only', ['insert_size']),
        ('Motifs Only', ['motifs']),
        ('Methylation Only', ['methylation']),
        ('Positions Only', ['positions']),
        ('Insert + Motifs', ['insert_size', 'motifs']),
        ('Insert + Methylation', ['insert_size', 'methylation']),
        ('Insert + Positions', ['insert_size', 'positions']),
        ('Motifs + Methylation', ['motifs', 'methylation']),
        ('Insert + Motifs + Methylation', ['insert_size', 'motifs', 'methylation']),
        ('All Features (No Positions)', ['insert_size', 'motifs', 'methylation']),
        ('All Features (With Positions)', ['insert_size', 'motifs', 'methylation', 'positions']),
    ]
    
    results = []
    best_score = 0
    best_features = None
    best_combination = None
    
    for name, groups in combinations:
        # Combine features
        features = []
        for group in groups:
            features.extend(feature_groups[group])
        
        X = df[features].values
        
        # Standardize
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        # Test both models
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        
        # Random Forest
        rf = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=5)
        rf_scores = cross_val_score(rf, X_scaled, y, cv=cv, scoring='f1')
        rf_f1 = rf_scores.mean()
        
        # Logistic Regression
        lr = LogisticRegression(random_state=42, max_iter=1000)
        lr_scores = cross_val_score(lr, X_scaled, y, cv=cv, scoring='f1')
        lr_f1 = lr_scores.mean()
        
        # Use best of both
        best_f1 = max(rf_f1, lr_f1)
        best_model = "RF" if rf_f1 > lr_f1 else "LR"
        
        print(f"{name:35s} | N={len(features):2d} | RF:{rf_f1:.3f} LR:{lr_f1:.3f} | Best:{best_model} {best_f1:.3f}")
        
        results.append({
            'combination': name,
            'n_features': len(features),
            'rf_f1': rf_f1,
            'lr_f1': lr_f1,
            'best_f1': best_f1,
            'best_model': best_model
        })
        
        # Track overall best
        if best_f1 > best_score:
            best_score = best_f1
            best_features = features
            best_combination = name
    
    print()
    print("=" * 70)
    print(f"BEST COMBINATION: {best_combination}")
    print(f"F1-Score: {best_score:.3f}")
    print(f"Number of features: {len(best_features)}")
    print("=" * 70)
    
    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")
    
    # Save best feature list
    with open(best_features_file, 'w') as f:
        json.dump({
            'combination': best_combination,
            'features': best_features,
            'f1_score': best_score
        }, f, indent=2)
    print(f"Best features saved to: {best_features_file}")

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: feature_selection.py <features_file> <output_file> <best_features_file>")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2], sys.argv[3])

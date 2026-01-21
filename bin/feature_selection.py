#!/usr/bin/env python3
"""
Test different feature combinations and select the best
Runs multiple replicates and reports mean ± SD
"""
import sys
import pandas as pd
import numpy as np
import random
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold
import json

N_REPLICATES = 5  # Change this to run more/fewer replicates

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
    
    # Store results for each combination across replicates
    all_results = {combo[0]: {'rf': [], 'lr': [], 'best': []} for combo in combinations}
    
    print("=" * 70)
    print(f"FEATURE COMBINATION TESTING ({N_REPLICATES} REPLICATES)")
    print("=" * 70)
    print(f"Dataset: {len(df)} samples ({sum(y)} ALS, {len(y)-sum(y)} Control)")
    print()
    
    # Run multiple replicates
    for rep in range(N_REPLICATES):
        random_seed = 42 + rep  # Different seed for each replicate
        np.random.seed(random_seed)
        random.seed(random_seed)
        
        print(f"--- Replicate {rep+1}/{N_REPLICATES} (Random Seed: {random_seed}) ---")
        
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
            cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_seed)
            
            # Random Forest
            rf = RandomForestClassifier(n_estimators=100, random_state=random_seed, max_depth=5, n_jobs=1)
            rf_scores = cross_val_score(rf, X_scaled, y, cv=cv, scoring='f1')
            rf_f1 = rf_scores.mean()
            
            # Logistic Regression
            lr = LogisticRegression(random_state=random_seed, max_iter=1000)
            lr_scores = cross_val_score(lr, X_scaled, y, cv=cv, scoring='f1')
            lr_f1 = lr_scores.mean()
            
            best_f1 = max(rf_f1, lr_f1)
            
            # Store results
            all_results[name]['rf'].append(rf_f1)
            all_results[name]['lr'].append(lr_f1)
            all_results[name]['best'].append(best_f1)
    
    print()
    print("=" * 70)
    print("SUMMARY: Mean ± SD across all replicates")
    print("=" * 70)
    print()
    
    results = []
    best_avg_score = 0
    best_combination_name = None
    
    for name, groups in combinations:
        # Compute actual feature count for this combination
        features = []
        for group in groups:
            features.extend(feature_groups[group])
        n_features = len(features)
        
        rf_scores = all_results[name]['rf']
        lr_scores = all_results[name]['lr']
        best_scores = all_results[name]['best']
        
        rf_mean, rf_sd = np.mean(rf_scores), np.std(rf_scores)
        lr_mean, lr_sd = np.mean(lr_scores), np.std(lr_scores)
        best_mean, best_sd = np.mean(best_scores), np.std(best_scores)
        
        # Determine which model is better on average
        best_model = "RF" if rf_mean > lr_mean else "LR"
        
        print(f"{name:35s} | N={n_features:2d}")
        print(f"  RF: {rf_mean:.3f} ± {rf_sd:.3f}")
        print(f"  LR: {lr_mean:.3f} ± {lr_sd:.3f}")
        print(f"  Best: {best_mean:.3f} ± {best_sd:.3f} ({best_model})")
        print()
        
        results.append({
            'combination': name,
            'n_features': n_features,
            'rf_f1_mean': rf_mean,
            'rf_f1_sd': rf_sd,
            'lr_f1_mean': lr_mean,
            'lr_f1_sd': lr_sd,
            'best_f1_mean': best_mean,
            'best_f1_sd': best_sd,
            'best_model': best_model
        })
        
        # Track overall best
        if best_mean > best_avg_score:
            best_avg_score = best_mean
            best_combination_name = name
    
    print("=" * 70)
    print(f"BEST COMBINATION (avg): {best_combination_name}")
    print(f"Average F1-Score: {best_avg_score:.3f}")
    print("=" * 70)
    
    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
    print(f"\nDetailed results saved to: {output_file}")
    
    # Save best feature list
    for combo_name, combo_groups in combinations:
        if combo_name == best_combination_name:
            best_features = []
            for group in combo_groups:
                best_features.extend(feature_groups[group])
            break
    
    with open(best_features_file, 'w') as f:
        json.dump({
            'combination': best_combination_name,
            'features': best_features,
            'f1_score_mean': best_avg_score,
            'n_replicates': N_REPLICATES
        }, f, indent=2)
    print(f"Best features saved to: {best_features_file}")

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: feature_selection.py <features_file> <output_file> <best_features_file>")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2], sys.argv[3])

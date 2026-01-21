#!/usr/bin/env python3
"""
Train final classifier using only the best feature combination
Runs multiple replicates and reports mean ± SD
"""
import sys
import pandas as pd
import numpy as np
import random
import json
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import confusion_matrix, classification_report

N_REPLICATES = 5  # Match feature_selection.py

def main(features_file, best_features_file, output_file, output_txt):
    # Load data
    df = pd.read_csv(features_file)
    
    # Load best features
    with open(best_features_file, 'r') as f:
        best_config = json.load(f)
    
    features = best_config['features']
    combination = best_config['combination']
    
    print(f"Using best feature combination: {combination}")
    print(f"Number of features: {len(features)}")
    print(f"Running {N_REPLICATES} replicates...\n")
    
    # Prepare data
    X = df[features].values
    y = (df['disease_status'] == 'als').astype(int)
    
    # Store results for each classifier across replicates
    results_by_classifier = {
        'Logistic Regression': {'accuracy': [], 'precision': [], 'recall': [], 'f1': []},
        'Random Forest': {'accuracy': [], 'precision': [], 'recall': [], 'f1': []}
    }
    
    output_lines = []
    output_lines.append(f"Dataset: {len(df)} samples")
    output_lines.append(f"  ALS: {sum(y)} samples")
    output_lines.append(f"  Control: {len(y) - sum(y)} samples")
    output_lines.append(f"\nFeature Combination: {combination}")
    output_lines.append(f"Total features: {len(features)}")
    output_lines.append(f"Replicates: {N_REPLICATES}\n")
    
    # Run multiple replicates
    for rep in range(N_REPLICATES):
        random_seed = 42 + rep  # Match feature_selection.py seeds
        np.random.seed(random_seed)
        random.seed(random_seed)
        
        # Standardize (fresh scaler each replicate)
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        # Train both models
        classifiers = {
            'Logistic Regression': LogisticRegression(random_state=random_seed, max_iter=1000),
            'Random Forest': RandomForestClassifier(n_estimators=100, random_state=random_seed, max_depth=5, n_jobs=1)
        }
        
        for name, clf in classifiers.items():
            cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_seed)
            y_pred = cross_val_predict(clf, X_scaled, y, cv=cv)
            
            accuracy = accuracy_score(y, y_pred)
            precision = precision_score(y, y_pred, zero_division=0)
            recall = recall_score(y, y_pred, zero_division=0)
            f1 = f1_score(y, y_pred, zero_division=0)
            
            results_by_classifier[name]['accuracy'].append(accuracy)
            results_by_classifier[name]['precision'].append(precision)
            results_by_classifier[name]['recall'].append(recall)
            results_by_classifier[name]['f1'].append(f1)
    
    # Compute means and standard deviations
    results = []
    for name in classifiers.keys():
        acc_mean, acc_sd = np.mean(results_by_classifier[name]['accuracy']), np.std(results_by_classifier[name]['accuracy'])
        prec_mean, prec_sd = np.mean(results_by_classifier[name]['precision']), np.std(results_by_classifier[name]['precision'])
        recall_mean, recall_sd = np.mean(results_by_classifier[name]['recall']), np.std(results_by_classifier[name]['recall'])
        f1_mean, f1_sd = np.mean(results_by_classifier[name]['f1']), np.std(results_by_classifier[name]['f1'])
        
        output_lines.append(f"{'='*60}")
        output_lines.append(f"{name}")
        output_lines.append(f"{'='*60}")
        output_lines.append(f"Accuracy:    {acc_mean:.3f} ± {acc_sd:.3f}")
        output_lines.append(f"Precision:   {prec_mean:.3f} ± {prec_sd:.3f}")
        output_lines.append(f"Sensitivity: {recall_mean:.3f} ± {recall_sd:.3f}")
        output_lines.append(f"F1-Score:    {f1_mean:.3f} ± {f1_sd:.3f}")
        output_lines.append("")
        
        results.append({
            'Classifier': name,
            'Feature_Combination': combination,
            'N_Features': len(features),
            'Accuracy_Mean': acc_mean,
            'Accuracy_SD': acc_sd,
            'Precision_Mean': prec_mean,
            'Precision_SD': prec_sd,
            'Sensitivity_Mean': recall_mean,
            'Sensitivity_SD': recall_sd,
            'F1_Mean': f1_mean,
            'F1_SD': f1_sd,
            'N_Replicates': N_REPLICATES
        })
    
    # Print and save
    output_text = '\n'.join(output_lines)
    print(output_text)
    
    with open(output_txt, 'w') as f:
        f.write(output_text)
    
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: classify_best_features.py <features_file> <best_features_file> <output_file> <output_txt>")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

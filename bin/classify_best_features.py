#!/usr/bin/env python3
"""
Train final classifier using only the best feature combination
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

# Set random seeds for reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
random.seed(RANDOM_SEED)

def main(features_file, best_features_file, output_file, output_txt):
    # Load data
    df = pd.read_csv(features_file)
    
    # Load best features
    with open(best_features_file, 'r') as f:
        best_config = json.load(f)
    
    features = best_config['features']
    combination = best_config['combination']
    
    print(f"Random Seed: {RANDOM_SEED}")
    print(f"Using best feature combination: {combination}")
    print(f"Number of features: {len(features)}")
    print()
    
    # Prepare data
    X = df[features].values
    y = (df['disease_status'] == 'als').astype(int)
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Train both models with explicit random seed
    classifiers = {
        'Logistic Regression': LogisticRegression(random_state=RANDOM_SEED, max_iter=1000),
        'Random Forest': RandomForestClassifier(n_estimators=100, random_state=RANDOM_SEED, max_depth=5)
    }
    
    results = []
    output_lines = []
    output_lines.append(f"Random Seed: {RANDOM_SEED}")
    output_lines.append(f"Dataset: {len(df)} samples")
    output_lines.append(f"  ALS: {sum(y)} samples")
    output_lines.append(f"  Control: {len(y) - sum(y)} samples")
    output_lines.append(f"\nFeature Combination: {combination}")
    output_lines.append(f"Total features: {len(features)}\n")
    
    for name, clf in classifiers.items():
        output_lines.append(f"{'='*60}")
        output_lines.append(f"{name}")
        output_lines.append(f"{'='*60}")
        
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)
        y_pred = cross_val_predict(clf, X_scaled, y, cv=cv)
        
        accuracy = accuracy_score(y, y_pred)
        precision = precision_score(y, y_pred, zero_division=0)
        recall = recall_score(y, y_pred, zero_division=0)
        f1 = f1_score(y, y_pred, zero_division=0)
        
        output_lines.append(f"Accuracy:    {accuracy:.3f}")
        output_lines.append(f"Precision:   {precision:.3f}")
        output_lines.append(f"Sensitivity: {recall:.3f}")
        output_lines.append(f"F1-Score:    {f1:.3f}")
        output_lines.append("")
        
        cm = confusion_matrix(y, y_pred)
        output_lines.append("Confusion Matrix:")
        output_lines.append(f"              Predicted")
        output_lines.append(f"              Control  ALS")
        output_lines.append(f"Actual Control   {cm[0,0]:3d}    {cm[0,1]:3d}")
        output_lines.append(f"       ALS       {cm[1,0]:3d}    {cm[1,1]:3d}")
        output_lines.append("")
        
        output_lines.append("Classification Report:")
        output_lines.append(classification_report(y, y_pred, target_names=['Control', 'ALS'], zero_division=0))
        
        results.append({
            'Classifier': name,
            'Feature_Combination': combination,
            'N_Features': len(features),
            'Accuracy': accuracy,
            'Precision': precision,
            'Sensitivity': recall,
            'F1-Score': f1
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

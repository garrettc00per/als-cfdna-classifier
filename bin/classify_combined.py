#!/usr/bin/env python3
"""
Binary classifier for ALS vs Control using combined features
"""
import sys
import pandas as pd
import numpy as np
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import confusion_matrix, classification_report

def main(features_file, output_file):
    # Load data
    df = pd.read_csv(features_file)
    
    # Separate features and labels
    # Exclude sample_id and disease_status columns
    feature_cols = [col for col in df.columns if col not in ['sample_id', 'disease_status']]
    X = df[feature_cols].values
    
    # Labels: als=1, ctrl=0
    y = (df['disease_status'] == 'als').astype(int)
    
    print(f"Dataset: {len(df)} samples")
    print(f"  ALS: {sum(y)} samples")
    print(f"  Control: {len(y) - sum(y)} samples")
    print(f"\nTotal features: {len(feature_cols)}")
    print()
    
    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Try multiple classifiers
    classifiers = {
        'Logistic Regression': LogisticRegression(random_state=42, max_iter=1000),
        'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42, max_depth=5)
    }
    
    results = []
    
    for name, clf in classifiers.items():
        print(f"{'='*60}")
        print(f"{name}")
        print(f"{'='*60}")
        
        # Use stratified K-fold cross-validation
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        
        # Get predictions via cross-validation
        y_pred = cross_val_predict(clf, X_scaled, y, cv=cv)
        
        # Calculate metrics
        accuracy = accuracy_score(y, y_pred)
        precision = precision_score(y, y_pred, zero_division=0)
        recall = recall_score(y, y_pred, zero_division=0)
        f1 = f1_score(y, y_pred, zero_division=0)
        
        print(f"Accuracy:    {accuracy:.3f}")
        print(f"Precision:   {precision:.3f}")
        print(f"Sensitivity: {recall:.3f}")
        print(f"F1-Score:    {f1:.3f}")
        print()
        
        # Confusion matrix
        cm = confusion_matrix(y, y_pred)
        print("Confusion Matrix:")
        print(f"              Predicted")
        print(f"              Control  ALS")
        print(f"Actual Control   {cm[0,0]:3d}    {cm[0,1]:3d}")
        print(f"       ALS       {cm[1,0]:3d}    {cm[1,1]:3d}")
        print()
        
        # Detailed classification report
        print("Classification Report:")
        print(classification_report(y, y_pred, target_names=['Control', 'ALS'], zero_division=0))
        
        results.append({
            'Classifier': name,
            'Accuracy': accuracy,
            'Precision': precision,
            'Sensitivity': recall,
            'F1-Score': f1
        })
    
    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: classify_combined.py <features_file> <output_file>")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2])

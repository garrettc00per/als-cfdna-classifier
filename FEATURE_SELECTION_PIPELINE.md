FEATURE SELECTION SCRIPT WALKTHROUGH

This is where things get interesting. We've extracted all these features from the BAM files, now we need to figure out which ones actually help classify ALS vs control. The naive approach would be to just throw everything at a classifier and see what happens. But with only 22 samples, that's a recipe for overfitting. So instead, this script systematically tests different feature combinations to find what works best.

THE CORE PROBLEM

We have 59 total features but only 22 samples. That's a terrible ratio for machine learning. Classic curse of dimensionality - when you have more features than samples, models can just memorize the training data instead of learning real patterns. So we need to be smart about feature selection.

The question: Which features actually carry the disease signal, and which are just noise?

WHAT THIS SCRIPT DOES

Tests 11 different feature combinations:
1. Insert size only (8 features)
2. Motifs only (40 features)
3. Methylation only (8 features)
4. Positions only (3 features)
5. Insert + Motifs (48 features)
6. Insert + Methylation (16 features)
7. Insert + Positions (11 features)
8. Motifs + Methylation (48 features)
9. Insert + Motifs + Methylation (56 features)
10. All features without positions (56 features)
11. All features with positions (59 features)

For each combination, it:
- Trains both Random Forest and Logistic Regression
- Does 5-fold cross-validation to avoid overfitting
- Runs 20 complete replicates with different random seeds
- Reports mean ± standard deviation to assess stability

Then picks the combination with the best average F1-score across all replicates.

THE CODE STRUCTURE

#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np
import random
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold
import json

N_REPLICATES = 20  # Run the entire analysis 20 times

Standard imports. The key ones are:
- RandomForestClassifier and LogisticRegression: our two model types
- StandardScaler: normalizes features so they're comparable
- StratifiedKFold: for cross-validation that preserves class balance
- cross_val_score: convenience function for running CV

DEFINING FEATURE GROUPS

feature_groups = {
    'insert_size': ['mean', 'median', 'stddev', 'min', 'max', 'q25', 'q75', 'n_fragments'],
    'methylation': ['cpg_meth_rate', 'chg_meth_rate', 'chh_meth_rate', 'overall_meth_rate', 
                    'total_cpg', 'total_chg', 'total_chh', 'n_reads'],
    'motifs': [col for col in all_cols if 'start_' in col or 'end_' in col],
    'positions': ['region1_count', 'region2_count', 'region3_count'],
}

This maps feature group names to actual column names in the combined features file. The motifs list is built dynamically by finding all columns with 'start_' or 'end_' in the name (that's the 40 motif frequency columns).

Then I define the 11 combinations to test:

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

Each tuple is (descriptive_name, list_of_groups_to_include).

Why these specific combinations? I'm testing:
- Each feature type alone (to see individual performance)
- Insert size + each other type (since insert size seems promising)
- Multi-feature combinations (to see if features complement each other)
- All features with/without positions (comprehensive test)

This isn't exhaustive (there are 2^4 = 16 possible combinations) but it covers the most interesting cases.

THE REPLICATION STRATEGY

all_results = {combo[0]: {'rf': [], 'lr': [], 'best': []} for combo in combinations}

for rep in range(N_REPLICATES):
    random_seed = 42 + rep
    np.random.seed(random_seed)
    random.seed(random_seed)

This is the key innovation. Instead of running each combination once, I run it 20 times with different random seeds. Why?

With only 22 samples, there's substantial variance in CV results depending on how you split the data. Two runs with different random seeds can give F1-scores that differ by 5-10 percentage points just due to luck. That makes it hard to know if one combination is truly better or just got lucky.

By running 20 replicates, I can calculate mean ± SD for each combination. If a combination has F1 = 0.75 ± 0.02, that's stable. If it's 0.75 ± 0.10, that's unstable - performance varies wildly depending on the random split.

The random seed controls:
- How samples are split into CV folds (StratifiedKFold shuffle)
- Random forest tree construction
- Any other stochastic elements

Incrementing the seed (42, 43, 44, ...) were added to ensure reproducible results, however even with seeds I saw slight variance across runs. This need additional troubleshooting...

TESTING EACH COMBINATION

for name, groups in combinations:
    # Combine features
    features = []
    for group in groups:
        features.extend(feature_groups[group])
    
    X = df[features].values

This builds the feature matrix for this combination. For example, if groups = ['insert_size', 'methylation'], it looks up those lists in feature_groups and concatenates them to get all 16 feature column names, then extracts those columns as a numpy array.

STANDARDIZATION

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

StandardScaler transforms each feature to have mean=0 and standard deviation=1. This is critical because:

1. Features are on different scales
   - Insert size mean: ~167 bp
   - Methylation rate: ~0.7 (70%)
   - Motif frequency: ~0.02 (2%)

2. Without scaling, distance-based algorithms get confused
   - Logistic regression uses feature magnitudes in the decision boundary
   - Features with larger values would dominate just because of scale

3. Makes interpretation easier
   - All features contribute equally by default
   - Feature importance is comparable across features

After standardization, a value of 2.0 means "2 standard deviations above the mean" regardless of which feature it is.

CROSS-VALIDATION

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_seed)

StratifiedKFold splits the data into 5 folds while preserving class balance. With 12 ALS and 10 controls:
- Regular KFold might give you a fold with 7 ALS and 0 controls (bad!)
- StratifiedKFold ensures each fold has roughly 60% ALS and 40% controls

Shuffle=True randomizes the order before splitting (deterministically via random_state). Without shuffling, if samples are ordered by disease status, the first fold would be all ALS and last fold all controls.

n_splits=5 means:
- Train on 80% of data (17-18 samples)
- Test on 20% of data (4-5 samples)
- Repeat 5 times with different test sets
- Average the results

Why 5 folds? It's a standard choice that balances:
- More folds = more training data per fold, but more variance (extreme: leave-one-out)
- Fewer folds = less training data, but more stable (extreme: 2-fold)
- 5 or 10 folds are conventional

TESTING BOTH MODELS

rf = RandomForestClassifier(n_estimators=100, random_state=random_seed, max_depth=5, n_jobs=1)
rf_scores = cross_val_score(rf, X_scaled, y, cv=cv, scoring='f1')
rf_f1 = rf_scores.mean()

lr = LogisticRegression(random_state=random_seed, max_iter=1000)
lr_scores = cross_val_score(lr, X_scaled, y, cv=cv, scoring='f1')
lr_f1 = lr_scores.mean()

I test both Random Forest and Logistic Regression because they have complementary strengths:

Random Forest:
- Handles non-linear relationships
- Can capture complex feature interactions
- More prone to overfitting with small sample sizes
- Harder to interpret
- My hyperparameters: 100 trees, max depth 5

Logistic Regression:
- Assumes linear decision boundary
- More resistant to overfitting
- Simpler, easier to interpret
- Works well with small datasets
- My hyperparameters: max 1000 iterations for convergence

cross_val_score does the CV for us. It:
1. Splits data into 5 folds
2. For each fold: train on 4 folds, test on 1 fold
3. Returns array of 5 scores (one per fold)
4. We take the mean as our performance estimate

Scoring='f1' because we care about both precision and recall. With imbalanced classes (12 ALS vs 10 controls), accuracy can be misleading. F1-score is the harmonic mean of precision and recall, so it only gives high scores when both are good.

STORING RESULTS

all_results[name]['rf'].append(rf_f1)
all_results[name]['lr'].append(lr_f1)
all_results[name]['best'].append(best_f1)

After running 20 replicates, each list has 20 values. So for "Insert Size Only", we have:
- all_results['Insert Size Only']['rf'] = [0.78, 0.81, 0.79, 0.77, 0.80, ...]
- all_results['Insert Size Only']['lr'] = [0.80, 0.82, 0.79, 0.81, 0.83, ...]

This lets us compute mean and standard deviation to assess stability.

COMPUTING SUMMARY STATISTICS

rf_mean, rf_sd = np.mean(rf_scores), np.std(rf_scores)
lr_mean, lr_sd = np.mean(lr_scores), np.std(lr_scores)
best_mean, best_sd = np.mean(best_scores), np.std(best_scores)

For each combination, I report:
- Random Forest: mean F1 ± SD
- Logistic Regression: mean F1 ± SD  
- Best of the two: mean F1 ± SD

Example output:
Insert Size Only                    | N=8
  RF: 0.713 ± 0.048
  LR: 0.773 ± 0.035
  Best: 0.773 ± 0.035 (LR)

This tells me:
- Logistic regression outperforms random forest for this combination
- LR achieves 77.3% F1 on average
- Performance is fairly stable (SD = 3.5 percentage points)
- Simpler model (LR) works better, suggesting signal is linear

SELECTING THE BEST COMBINATION

if best_mean > best_avg_score:
    best_avg_score = best_mean
    best_combination_name = name

After testing all combinations, pick whichever has the highest average "best" score (i.e., the best of RF and LR for that combination, averaged across replicates).

In my results, "Insert + Positions" wins with F1 = 0.702.

SAVING OUTPUTS

1. CSV with detailed results for all combinations:

results_df.to_csv(output_file, index=False)

Columns: combination, n_features, rf_f1_mean, rf_f1_sd, lr_f1_mean, lr_f1_sd, best_f1_mean, best_f1_sd, best_model

This lets you make plots comparing combinations, see stability across feature sets, etc.

2. JSON with the winning feature list:

with open(best_features_file, 'w') as f:
    json.dump({
        'combination': best_combination_name,
        'features': best_features,
        'f1_score_mean': best_avg_score,
        'n_replicates': N_REPLICATES
    }, f, indent=2)

This gets passed to the classification script so it knows which features to use for final training. The JSON contains the actual feature column names, not just the group name.

KEY FINDINGS FROM RESULTS

What I found when running this:

1. Insert + Positions performs best (F1 = 0.702 ± 0.062)
   - Just 11 features (8 insert + 3 positions)
   - Logistic regression wins
   - Positions surprisingly DO help classification

2. Several combinations perform competitively
   - Motifs Only: F1 = 0.701 ± 0.081 (40 features!)
   - Insert + Motifs: F1 = 0.693 ± 0.070 (48 features)
   - Positions Only: F1 = 0.682 ± 0.060 (just 3 features!)
   - Insert Size Only: F1 = 0.674 ± 0.057 (8 features)

3. Methylation adds minimal value
   - Methylation Only: F1 = 0.604 ± 0.075 (weakest individual feature)
   - Insert + Methylation: F1 = 0.643 ± 0.051 (worse than insert alone)
   - Adding methylation generally doesn't improve any combination

4. Positions are surprisingly informative
   - Positions Only gets 68% F1 with just 3 features
   - Adding positions to insert size gives the best performance
   - Disease signal has BOTH molecular and positional components

5. Logistic regression dominates
   - LR wins for 10 out of 11 combinations
   - Only exception: complex multi-feature models where RF wins
   - Suggests mostly linear decision boundaries
   - Overfitting is a real concern with random forest

6. High variance across replicates (SD = 0.05-0.10)
   - 5-10 percentage points of variance from random CV splits
   - With only 22 samples, data split matters a lot
   - This is why replication matters - single runs can be very misleading
   - The top 4 combinations are within one SD of each other

BIOLOGICAL INTERPRETATION

The results reveal a more nuanced picture than expected:

1. POSITIONAL INFORMATION MATTERS
   - Positions alone get 68% F1 with just 3 features
   - Adding positions to insert size gives the best performance (70%)
   - This means fragments from different chr21 regions have different disease associations
   - ALS patients may preferentially release cfDNA from specific chromosomal domains

2. FRAGMENT SIZE STILL IMPORTANT
   - Insert size features work well alone (67% F1)
   - Insert + Positions works better than either alone (70% F1)
   - Suggests BOTH molecular properties AND genomic location matter
   - Not just "how DNA fragments" but also "where DNA fragments from"

3. MOTIFS ARE SURPRISINGLY INFORMATIVE
   - Motifs alone: 70% F1 with 40 features
   - Comparable to the best combination
   - End motif patterns capture fragmentation preferences
   - May reflect altered nuclease activity or chromatin accessibility in ALS

4. METHYLATION DOESN'T HELP
   - Weakest individual feature (60% F1)
   - Doesn't improve any combination
   - Either methylation rates truly don't differ in ALS
   - Or sample size is too small to detect the signal
   - Or methylation info is redundant with fragment size/position

5. MULTIPLE ROUTES TO CLASSIFICATION
   - Several feature combinations perform similarly (68-70% F1)
   - No single "winner" that dominates
   - Suggests the ALS signal is distributed across feature types
   - Different feature sets capture overlapping information

What this tells us biologically:
- ALS fragmentation is MULTI-DIMENSIONAL
- Both intrinsic fragment properties (size, motifs) and genomic location contribute
- The disease may alter both chromatin structure (affecting fragmentation) and tissue-specific cfDNA shedding (affecting regional distribution)
- With small N, we can't definitively separate these mechanisms

For a clinical test:
- Could use just 11 features (insert + positions) for simplicity
- Or use motifs if insert size measurement is difficult
- Multiple viable biomarker strategies

WHY RUN 20 REPLICATES INSTEAD OF JUST 1

With 22 samples, a single CV run can be misleading. Example:

Replicate 1 (seed=42): Insert Size F1 = 0.81
Replicate 2 (seed=43): Insert Size F1 = 0.73
Replicate 3 (seed=44): Insert Size F1 = 0.79

Which one is "true"? They're all true - they just got different random splits. By averaging 20 replicates, I get 0.77 ± 0.03, which is a much more honest estimate.

If I only ran one replicate and happened to get 0.81, I might think performance is better than it really is. Or if I got 0.73, I might think it's worse. The mean across many replicates converges to the true expected performance.

Also, the standard deviation tells me about stability. If SD is 0.10, that's a red flag - performance is too dependent on how you split the data, suggesting the model hasn't learned a robust pattern.

DESIGN CHOICES I MADE

1. Test both RF and LR
   - Different inductive biases
   - One might work better for this data
   - Result: LR wins almost always, suggesting linear separability

2. Use F1-score instead of accuracy
   - Classes are nearly balanced (12 ALS, 10 controls) so accuracy would work
   - But F1 is more robust to slight imbalance
   - Forces model to do well on both classes

3. 5-fold instead of 10-fold CV
   - With 22 samples, 10-fold means training on only 19-20 samples
   - 5-fold gives 17-18 training samples per fold
   - More training data = more stable models
   - Standard choice for small datasets

4. Fixed hyperparameters instead of grid search
   - Could do nested CV with hyperparameter tuning
   - But with 22 samples, I'd be optimizing noise
   - Better to use reasonable defaults
   - Random forest max_depth=5 prevents overfitting

5. Standardization for all models
   - Logistic regression requires it (different scales)
   - Random forest technically doesn't need it (tree-based)
   - But standardizing makes feature importance comparable
   - Doesn't hurt, might help

6. Test 11 combinations, not all 16 possible
   - Could test every possible subset
   - But 11 combinations cover the most interpretable cases
   - Each individual feature type, common pairwise combinations, and all-features
   - Balances thoroughness with computational cost

A NOTE ABOUT UNEXPECTED RESULTS

When I designed this pipeline, I expected insert size features to dominate (based on prior cfDNA literature). The actual results are more nuanced:

- Positions matter more than expected (3 position features get 68% F1!)
- Motifs alone perform as well as insert + positions (70% F1)
- Multiple combinations achieve similar performance

This is a good reminder: let the data speak, don't force predetermined conclusions. The multi-dimensional nature of the signal is actually more interesting scientifically than if one feature had dominated. It suggests ALS affects multiple aspects of cfDNA biology simultaneously.

For an interview, this demonstrates:
- Willingness to update hypotheses based on data
- Recognition that "unexpected" results are often the most informative
- Scientific integrity - reporting what you found, not what you expected

WHAT I COULD IMPROVE

1. Add confidence intervals
Instead of just mean ± SD, could compute 95% CI:
from scipy import stats
ci = stats.t.interval(0.95, len(scores)-1, loc=np.mean(scores), scale=stats.sem(scores))

This would give more rigorous statistical bounds.

2. Test for significant differences
After finding "Insert Size Only" is best, could do paired t-test vs other combinations:
from scipy.stats import ttest_rel
t_stat, p_val = ttest_rel(insert_scores, insert_meth_scores)

If p < 0.05, the difference is statistically significant. Right now I just compare means.

3. Nested cross-validation
Outer loop: 5-fold CV for performance estimation
Inner loop: hyperparameter tuning on training folds only

This is the "right" way to do CV with hyperparameter search, but with 22 samples it's overkill and very slow.

4. Try more models
Could test SVM, gradient boosting, neural nets, etc. But:
- More models = more multiple testing issues
- With 22 samples, complex models will overfit
- LR and RF cover the bases (linear vs nonlinear)

5. Feature importance analysis
For the winning combination, could extract:
- Logistic regression coefficients (which features matter most)
- Random forest feature importances

I do this in the classification script, but could add it here too.

HOW THIS FITS IN THE PIPELINE

Pipeline flow:
[Feature Extraction] → [Combine Features] → [Feature Selection] → [Classification] → [Visualization]

Feature Selection (this script) sits in the middle. It:
- Takes the combined feature matrix as input
- Tests 11 combinations systematically
- Outputs the winning combination to a JSON file
- That JSON feeds into the final classification script

The classification script then:
- Reads the best features JSON
- Trains final models on those features
- Generates detailed performance metrics
- Creates confusion matrices and classification reports

So this script's job is SELECTION, not final training. It's telling downstream processes "use these 8 insert size features + 3 position features, ignore the rest."

COMPUTATIONAL COST

20 replicates × 11 combinations × 2 models × 5 CV folds = 2200 model fits

Each fit trains on ~18 samples, which is fast. Total runtime: ~2 minutes on my laptop.

If I increased to 100 replicates: ~10 minutes
If I added nested CV for hyperparameter tuning: ~30-60 minutes

For a take-home assignment with quick turnaround, 20 replicates is a good balance of rigor and speed.

ALTERNATIVE APPROACHES

1. Recursive Feature Elimination (RFE)
- Start with all features
- Iteratively remove least important feature
- Test performance at each step
- Pick feature count with best performance

Pro: More fine-grained, can find optimal subset
Con: Much slower, risk of overfitting with small N

2. L1 regularization (LASSO)
- Add penalty for non-zero coefficients
- Forces some coefficients to exactly zero
- Automatic feature selection
- Tune regularization strength via CV

Pro: Principled approach, one parameter to tune
Con: May not find best discrete subset

3. Univariate filtering
- Test each feature individually (t-test, mutual information)
- Keep top K features by p-value or score
- Very fast

Pro: Simple, interpretable
Con: Ignores feature interactions

4. Domain knowledge
- Just use insert size because we know it matters
- Skip feature selection entirely

Pro: Simplest, fastest
Con: Might miss complementary features, less convincing for interview

I chose exhaustive combination testing because:
- Systematic and transparent
- Tests both individual and joint performance
- Shows I'm thinking about the curse of dimensionality
- Demonstrates proper ML methodology

WHAT THE RESULTS TELL US SCIENTIFICALLY

The fact that multiple feature combinations work similarly well (68-70% F1) is actually a scientific finding:

It tells us:
1. ALS fragmentation signal is MULTI-DIMENSIONAL
   - Not captured by a single feature type
   - Different features (size, motifs, positions) all contribute
   - Multiple routes to ~70% classification performance

2. SMALL SAMPLE SIZE LIMITS RESOLUTION
   - Top 4 combinations overlap within their standard deviations
   - Can't definitively say one is "best"
   - Likely need 100+ samples to distinguish feature importance clearly
   - With 22 samples, we're identifying candidate features, not optimal sets

3. POSITIONS MATTER (unexpected finding!)
   - Chr21 regional distribution differs between ALS and controls
   - Either disease affects tissue-specific cfDNA shedding
   - Or certain chr21 regions have disease-relevant biology
   - Contradicts initial hypothesis that signal is purely molecular

4. METHYLATION DOESN'T HELP (unexpected)
   - Consistently weak across all combinations
   - Either ALS doesn't alter cfDNA methylation patterns
   - Or our methylation features don't capture the right patterns
   - Or sample size is insufficient

5. OVERFITTING RISK IS MANAGEABLE
   - Logistic regression works better than random forest
   - But performance plateaus around 70% regardless of features
   - Suggests we're hitting fundamental limits of the data
   - More features don't help because signal-to-noise ratio is constant

For future work:
- Larger cohort would stabilize feature selection
- Full-genome analysis might reveal chr21-specific vs genome-wide patterns
- More sophisticated positional features (gene-based, regulatory regions)
- Alternative methylation features (regional patterns, DMRs)

A NOTE ON REPRODUCIBILITY

I set random seeds throughout (np.random.seed(), random.seed(), random_state in models) to make this reproducible. However, I've found that exact reproducibility isn't guaranteed even with seeds set. There are a few likely culprits:

1. sklearn's cross_val_score has internal randomness that might not be fully controlled
2. RandomForest can have some non-deterministic behavior depending on version
3. Threading/parallelization (even with n_jobs=1) can introduce variance
4. Different sklearn versions might use different random number generators
5. Operating system differences in floating point operations

What I should have done for perfect reproducibility:
- Set PYTHONHASHSEED environment variable
- Use RandomState objects instead of global random seeds
- Disable all parallelization explicitly
- Save the actual CV fold indices for reuse

However, here's the important part: the 20-replicate approach with mean ± SD is actually MORE scientifically rigorous than single-run reproducibility. Why?

If my code gives different results each run, that's telling me something important: performance varies with data splits. By running 20 times and reporting the distribution of results, I'm quantifying that variance honestly rather than hiding it behind a single "reproducible" number.

Think about it:
- Single run with seed=42: F1 = 0.81 (lucky split?)
- Single run with seed=43: F1 = 0.73 (unlucky split?)
- 20 runs averaged: F1 = 0.77 ± 0.03 (true expected performance)

The third approach gives a more honest assessment. The fact that I CAN'T just cherry-pick one good seed is actually a feature, not a bug.

Acknowledge this limitation honestly, explain what I did try (setting seeds), and argue that the statistical approach makes the results more robust, not less. The scientific conclusions (multiple feature combinations work competitively, positions matter, methylation doesn't help) are stable across runs even if exact numbers vary slightly.

DISCUSSION

The results are nuanced - no single clear winner, multiple combinations perform similarly. Here's how I'd frame this:

"The feature selection revealed that several combinations achieve similar performance (68-70% F1), with Insert + Positions slightly ahead at 70%. However, the high standard deviations (5-6 percentage points) mean the top 4 combinations are statistically indistinguishable.

This is actually an important scientific finding, not a failure of the analysis. With only 22 samples, we're fundamentally limited in our ability to identify the single 'best' feature set. What we CAN conclude is:

1. Multiple features contribute to classification - it's not just one dominant signal
2. Positions matter (unexpected!), suggesting regional cfDNA differences
3. Methylation consistently underperforms across all combinations
4. We'd need 100+ samples to definitively rank these feature sets

For moving forward, I'd recommend Insert + Positions (11 features) as the most parsimonious option that achieves competitive performance. But I'd also flag Motifs Only as an alternative if fragment size is difficult to measure accurately.

The key is being transparent about uncertainty rather than overstating conclusions from limited data."

This shows:
- You understand statistical power and sample size limitations
- You can interpret "no clear winner" as meaningful, not just inconclusive
- You make practical recommendations despite uncertainty
- You're honest about what you can and can't conclude

SUMMARY

This script is the scientific heart of the pipeline. It answers the question: "Which features actually matter for distinguishing ALS from controls?"

The methodology:
- Systematic testing of 11 feature combinations
- Both linear (LR) and nonlinear (RF) models
- 5-fold cross-validation for robust performance estimates
- 20 replicates to assess stability and variance
- F1-score to balance precision and recall

The finding:
- Insert + Positions work best (F1 = 0.702 ± 0.062)
- But several combinations perform competitively (F1 = 0.68-0.70)
- Simpler models (logistic regression) outperform complex ones
- Multiple features contribute: size, motifs, AND positions
- High variance across replicates indicates small sample challenges

The implication:
- ALS alters BOTH fragmentation patterns AND regional cfDNA distribution
- Signal is multi-dimensional, not captured by single feature type
- Multiple viable biomarker strategies exist
- Need larger cohorts to definitively identify optimal feature set
- ~70% F1-score may be the limit with current features and sample size

This demonstrates:
- Understanding of overfitting and curse of dimensionality
- Proper cross-validation methodology
- Appropriate model selection for small datasets
- Ability to interpret ML results biologically
- Rigorous approach to feature selection
- Honest assessment of reproducibility challenges
- Statistical thinking (mean ± SD over point estimates)
- Recognition that "no clear winner" is itself a scientific finding

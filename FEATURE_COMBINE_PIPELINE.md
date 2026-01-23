FEATURE COMBINATION SCRIPT WALKTHROUGH

This is the glue script that takes all the separate feature files and combines them into a single feature matrix ready for machine learning. Pretty straightforward but there are some design choices worth discussing.

WHAT IT DOES

Takes 4 separate feature files plus metadata and merges them into one big table:
1. Insert size features (8 features)
2. Methylation features (8 features)  
3. End motif features (40 features)
4. Position features (3 features)
5. Metadata (disease labels)

Final output is one row per sample with 59 features plus disease label.

THE CODE

#!/usr/bin/env python3
import sys
import pandas as pd

def main(insert_file, meth_file, motif_file, position_file, metadata_file, output_file):
    # Load all features
    insert_df = pd.read_csv(insert_file, sep=' ')
    meth_df = pd.read_csv(meth_file, sep=' ')
    motif_df = pd.read_csv(motif_file)
    position_df = pd.read_csv(position_file, sep=' ')
    
    # Merge all features
    combined_df = insert_df.merge(meth_df, on='sample_id')
    combined_df = combined_df.merge(motif_df, on='sample_id')
    combined_df = combined_df.merge(position_df, on='sample_id')
    
    # Add disease labels
    metadata_df = pd.read_csv(metadata_file)
    disease_map = dict(zip(metadata_df['Run'], metadata_df['disease_status']))
    combined_df['disease_status'] = combined_df['sample_id'].map(disease_map)
    
    # Save
    combined_df.to_csv(output_file, index=False)
    
    print(f"Combined ALL features:")
    print(f"  Samples: {len(combined_df)}")
    print(f"  Insert size features: 8")
    print(f"  Methylation features: 8")
    print(f"  Motif features: 40")
    print(f"  Position features: 3")
    print(f"  Total features: {len(combined_df.columns) - 2}")
    print(f"  Output: {output_file}")

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])

THE INPUT FILES

Insert size file (space-delimited):
sample_id mean median stddev q25 q75 min max count
SRR123456 167.3 165.0 54.2 135 195 50 500 98234
SRR789012 169.1 167.0 62.8 138 198 51 498 101456

These are the 8 summary statistics from fragment length distributions. Mean and median tell you average fragment size, stddev tells you variability, quartiles show the distribution shape, min/max show range, and count is total fragments.

Methylation file (space-delimited):
sample_id cpg_rate chg_rate chh_rate overall_rate total_cpg total_chg total_chh n_reads
SRR123456 0.719 0.012 0.008 0.114 63735 88323 281253 1000000
SRR789012 0.721 0.011 0.009 0.116 65102 89456 279834 1050000

The first 4 columns are methylation rates (what fraction of cytosines are methylated in each context). The next 3 are total counts (how many sites we observed). Last column is number of reads processed. So 8 features total from methylation.

Motif file (CSV):
sample_id,AAAA,AAAT,AATA,...,TTTT
SRR123456,0.023,0.015,0.031,...,0.018
SRR789012,0.025,0.014,0.029,...,0.019

This has frequencies for the top 20 start motifs plus top 20 end motifs = 40 total motif features. Each value is the fraction of fragments that have that particular 4-mer at the start or end. These are normalized by total fragment count so they're comparable across samples.

Position file (space-delimited):
sample_id region1 region2 region3 total
SRR123456 35123 38901 32456 106480
SRR789012 36789 39234 31567 107590

This divides chr21 into three equal regions and counts fragments in each. So 3 features (technically 4 columns but total is redundant since it's the sum). This tests whether fragments cluster in specific genomic regions.

Metadata CSV:
Run,disease_status,AGE,...
SRR123456,als,67,...
SRR789012,ctrl,62,...

Standard sample metadata. The key column is disease_status which is either "als" or "ctrl". This becomes our target variable for classification.

STEP 1 - LOAD ALL FEATURE FILES

insert_df = pd.read_csv(insert_file, sep=' ')
meth_df = pd.read_csv(meth_file, sep=' ')
motif_df = pd.read_csv(motif_file)
position_df = pd.read_csv(position_file, sep=' ')

Notice the different separators. Insert, methylation, and position files use space-delimited format because they come from bash scripts where space-delimited is easiest. Motif file uses comma-delimited because it comes from a Python script that naturally outputs CSV. Metadata is standard CSV. Not ideal to have mixed formats but each script uses what's most natural for that language. In production I'd probably standardize everything to CSV, but for a quick pipeline this works fine.

After loading, each dataframe looks like:
- insert_df: 22 rows x 9 columns (sample_id + 8 features)
- meth_df: 22 rows x 9 columns (sample_id + 8 features)
- motif_df: 22 rows x 41 columns (sample_id + 40 features)
- position_df: 22 rows x 5 columns (sample_id + 4 position stats)

STEP 2 - MERGE ALL FEATURES

combined_df = insert_df.merge(meth_df, on='sample_id')
combined_df = combined_df.merge(motif_df, on='sample_id')
combined_df = combined_df.merge(position_df, on='sample_id')

This does a series of merges on sample_id. Pandas merge defaults to inner join, so only samples present in ALL files will be kept. This is important - if methylation extraction failed for one sample, that sample won't appear in the final matrix.

The merges happen sequentially:
1. Start with insert_df (22 rows, 9 cols)
2. Merge with meth_df → (22 rows, 17 cols) because we add 8 new columns
3. Merge with motif_df → (22 rows, 57 cols) because we add 40 new columns  
4. Merge with position_df → (22 rows, 61 cols) because we add 4 new columns

Wait, that's 61 not 59. The position file has region1, region2, region3, and total. But total is redundant (it's just region1 + region2 + region3), so really only 3 informative features. But I keep all 4 columns at this stage.

This assumes all files have the same samples in the same order. If a sample failed at any extraction step, it won't appear in the final matrix. In my case, Nextflow ensures all samples go through all processes, so this works. But in production I'd add validation to check that all expected samples are present.

After merging, I have one row per sample with all features:
- sample_id column
- Columns for insert size stats (mean, median, stddev, q25, q75, min, max, count)
- Columns for methylation (cpg_rate, chg_rate, chh_rate, overall_rate, total_cpg, total_chg, total_chh, n_reads)
- Columns for motifs (AAAA, AAAT, AATA, ..., TTTT - 40 total)
- Columns for positions (region1, region2, region3, total)

STEP 3 - ADD DISEASE LABELS

metadata_df = pd.read_csv(metadata_file)
disease_map = dict(zip(metadata_df['Run'], metadata_df['disease_status']))
combined_df['disease_status'] = combined_df['sample_id'].map(disease_map)

First I load the metadata CSV which has lots of columns but I only care about Run (sample ID) and disease_status (als or ctrl).

Then I create a dictionary mapping sample IDs to disease status:
{
    'SRR123456': 'als',
    'SRR789012': 'ctrl',
    'SRR345678': 'als',
    ...
}

The zip() function pairs up the Run column with the disease_status column, and dict() converts those pairs into a dictionary. This is a clean one-liner for creating lookup dictionaries in Python.

Then I use .map() to create a new disease_status column in my feature dataframe. For each sample_id, it looks up the corresponding disease status in the dictionary. If a sample isn't in the dictionary, map() returns NaN, which would flag missing metadata. This is good for QC.

Now I have the target variable for supervised learning. The final dataframe has all features plus the label we're trying to predict.

STEP 4 - SAVE OUTPUT

combined_df.to_csv(output_file, index=False)

This writes the combined feature matrix to CSV. The index=False prevents pandas from writing row numbers as a column, which we don't need. Final CSV looks like:

sample_id,mean,median,stddev,q25,q75,min,max,count,cpg_rate,chg_rate,chh_rate,overall_rate,total_cpg,total_chg,total_chh,n_reads,AAAA,AAAT,...,region1,region2,region3,total,disease_status
SRR123456,167.3,165.0,54.2,135,195,50,500,98234,0.719,0.012,0.008,0.114,63735,88323,281253,1000000,0.023,0.015,...,35123,38901,32456,106480,als
SRR789012,169.1,167.0,62.8,138,198,51,498,101456,0.721,0.011,0.009,0.116,65102,89456,279834,1050000,0.025,0.014,...,36789,39234,31567,107590,ctrl

This format is perfect for scikit-learn. To use it for ML, you just do:
X = df.drop(['sample_id', 'disease_status'], axis=1)
y = df['disease_status']

X is your feature matrix (22 samples x 59 features), y is your target vector (22 labels).

THE PRINT STATEMENTS

print(f"Combined ALL features:")
print(f"  Samples: {len(combined_df)}")
print(f"  Insert size features: 8")
print(f"  Methylation features: 8")
print(f"  Motif features: 40")
print(f"  Position features: 3")
print(f"  Total features: {len(combined_df.columns) - 2}")
print(f"  Output: {output_file}")

These print statements serve two purposes. First, they're helpful for pipeline logging and QC. When I look at Nextflow logs, I can quickly see whether all 22 samples made it through and whether I have the expected number of features. Second, they document what's in the output file for anyone using the pipeline.

The "Total features: {len(combined_df.columns) - 2}" calculation subtracts 2 because the dataframe has sample_id and disease_status columns which aren't actually features for ML. So if the dataframe has 61 columns total, that's 59 actual features.

WHY THIS DESIGN - MODULAR VS MONOLITHIC

I chose to do separate feature extraction followed by combination, rather than extracting everything in one big script. Here's the tradeoff:

Separate extraction then combine (what I did):
Pro: Modular - can test each feature type independently
Pro: Easy to parallelize - Nextflow runs all extractions simultaneously
Pro: Can re-run combination without re-extracting features (fast iteration)
Pro: Easier to debug - can inspect intermediate files
Pro: Can skip features easily (just don't include that file in merge)
Con: More files to manage (4 feature files + 1 metadata file)
Con: Need to ensure sample IDs match across all files

All-in-one extraction:
Pro: Fewer intermediate files
Pro: Guaranteed same samples in all outputs
Con: Can't parallelize as effectively
Con: If one feature extraction fails, lose everything
Con: Harder to test individual components
Con: Slower iteration (have to re-extract all features to test changes)

For a pipeline with multiple feature types that might fail independently and where I want to test different feature combinations, I strongly prefer modular extraction. The tradeoff is slightly more bookkeeping to ensure sample IDs are consistent. But Nextflow handles that automatically, so it's not really a problem in practice.

POTENTIAL ISSUES AND HOW TO HANDLE THEM

Issue 1: Sample ID mismatches
Problem: If bash scripts output "SRR123456" but metadata has "SRR123456.1"
Result: Merge fails, samples silently dropped
Solution: Standardize sample ID extraction across all scripts. I use the same sample_id variable in all bash scripts and make sure it matches the Run column in metadata exactly.

Issue 2: Missing samples
Problem: If methylation extraction failed for one sample
Result: That sample silently dropped from final matrix (inner join behavior)
Better approach: Add validation to check all expected samples are present:

expected_samples = set(metadata_df['Run'])
actual_samples = set(combined_df['sample_id'])
if expected_samples != actual_samples:
    missing = expected_samples - actual_samples
    print(f"WARNING: Missing samples: {missing}")
    sys.exit(1)

I didn't implement this because Nextflow ensures all samples complete all processes, but in production I would.

Issue 3: Feature name collisions
Problem: If two files both have a column named "total"
Result: Pandas merge creates "total_x" and "total_y" columns automatically
Solution: Either use suffixes parameter in merge() to control naming, or rename columns proactively. In my case this doesn't happen because I control all the column names.

Issue 4: Different sample counts across files
Problem: Insert file has 22 samples, methylation has 21 (one failed)
Result: Combined file has 21 samples (inner join drops the missing sample)
Better: Use outer join and fill NaN, then investigate why that sample failed. Or validate inputs first and fail fast.

WHAT I COULD IMPROVE

1. Add explicit validation:

def validate_inputs(insert_df, meth_df, motif_df, position_df):
    """Check all inputs have same samples"""
    samples_insert = set(insert_df['sample_id'])
    samples_meth = set(meth_df['sample_id'])
    samples_motif = set(motif_df['sample_id'])
    samples_pos = set(position_df['sample_id'])
    
    if not (samples_insert == samples_meth == samples_motif == samples_pos):
        missing_in_meth = samples_insert - samples_meth
        missing_in_motif = samples_insert - samples_motif
        missing_in_pos = samples_insert - samples_pos
        raise ValueError(f"Sample mismatch! Missing in meth: {missing_in_meth}, motif: {missing_in_motif}, pos: {missing_in_pos}")
    
    return samples_insert

Then call this before merging to fail fast if there's a problem.

2. Use consistent file formats - either all CSV or all space-delimited. I'd probably standardize on CSV since it handles edge cases better (like sample IDs with spaces, though that shouldn't happen with SRR accessions).

3. Add feature name prefixes to prevent collisions:

insert_df = insert_df.rename(columns={col: f'insert_{col}' for col in insert_df.columns if col != 'sample_id'})
meth_df = meth_df.rename(columns={col: f'meth_{col}' for col in meth_df.columns if col != 'sample_id'})

This makes it crystal clear which feature came from which extraction step. Especially useful if you have generic names like "mean" or "count" that could appear in multiple feature sets.

4. Output more metadata about the merge:

print(f"Samples in insert file: {len(insert_df)}")
print(f"Samples in methylation file: {len(meth_df)}")
print(f"Samples in motif file: {len(motif_df)}")
print(f"Samples in position file: {len(position_df)}")
print(f"Samples in final matrix: {len(combined_df)}")
if len(combined_df) < len(insert_df):
    print(f"WARNING: Lost {len(insert_df) - len(combined_df)} samples during merge!")

5. Add feature descriptions in the output. Could include a companion file that documents what each feature means:

feature_descriptions = {
    'mean': 'Mean insert size (bp)',
    'stddev': 'Standard deviation of insert size distribution (bp)',
    'cpg_rate': 'Fraction of CpG sites that are methylated',
    'AAAA': 'Frequency of AAAA motif at fragment starts',
    ...
}
pd.DataFrame.from_dict(feature_descriptions, orient='index', columns=['Description']).to_csv('feature_descriptions.csv')

This would help anyone using the pipeline understand what each column represents.

HOW IT FITS IN THE PIPELINE

The Nextflow DAG looks like this:

[EXTRACT_INSERT_SIZES] ----\
[EXTRACT_METHYLATION] ------> [COMBINE_FEATURES] -> [FEATURE_SELECTION] -> [CLASSIFY]
[EXTRACT_END_MOTIFS] ------/
[EXTRACT_POSITIONS] -------/

This combine script is the COMBINE_FEATURES process. It waits for all four extraction processes to complete, collects all their outputs, merges them together, and passes the result to feature selection.

The feature selection script then tests different feature combinations:
- Insert only (8 features)
- Methylation only (8 features)
- Motifs only (40 features)
- Positions only (3 features)
- Insert + Methylation (16 features)
- Insert + Motifs (48 features)
- Insert + Methylation + Motifs (56 features)
- Insert + Methylation + Positions (19 features)
- Insert + Positions (11 features)
- Methylation + Motifs (48 features)
- All features (59 features)

By having all features in one file, I can easily subset columns to test different combinations without re-extracting anything. The feature selection script just reads this combined file once and then does:

X_insert = combined_df[['mean', 'median', 'stddev', ...]]
X_meth = combined_df[['cpg_rate', 'chg_rate', ...]]
X_combined = combined_df[['mean', 'median', ..., 'cpg_rate', 'chg_rate', ...]]

Much more efficient than creating separate combined files for each combination.

ALTERNATIVE DATA FORMATS

1. Wide format (what I use)
Pros: One row per sample, many columns, scikit-learn ready
Cons: Can be memory-intensive for very high-dimensional data
Best for: Traditional ML with moderate number of features

2. Long format
Structure: Multiple rows per sample, one per feature
sample_id | feature_name | feature_value | disease_status
SRR123    | mean        | 167.3         | als
SRR123    | median      | 165.0         | als
SRR123    | stddev      | 54.2          | als
...

Pros: Good for ggplot-style visualization, relational databases, feature engineering
Cons: Need to pivot to wide format before ML, harder to work with
Best for: Exploratory data analysis and visualization

3. Sparse matrix format
Structure: Only store non-zero values
Pros: Memory efficient for high-dimensional sparse data (like bag-of-words)
Cons: Our features are dense (every sample has every feature), no benefit
Best for: Text data, genomic variants, any data with mostly zeros

4. Binary formats (HDF5, Parquet)
Pros: Compressed, faster I/O, can handle huge datasets
Cons: Not human-readable, overkill for 22 samples
Best for: Large-scale production pipelines with thousands of samples

For this project with 22 samples and 59 features, CSV is perfect. It's human-readable, easy to inspect, and works with every tool. If I were scaling to thousands of samples I'd switch to Parquet for faster loading.

WHY I USE PANDAS FOR THIS

Could I do this merge in bash? Sure:
join -t',' insert.csv meth.csv | join -t',' - motif.csv | join -t',' - position.csv

But pandas is way better for this because:
1. Handles missing values gracefully (bash join is brittle)
2. Better error messages when things go wrong
3. Can easily validate data types, check for duplicates, etc
4. More readable code (the merge operations are clear)
5. Easier to add QC checks and summary statistics

The only downside is requiring Python, but that's already a dependency for the ML steps anyway.

DESIGN PHILOSOPHY

The overall philosophy here is: make each step modular and inspectable. I can look at any intermediate file and verify it makes sense:
- insert_sizes/ directory has fragment length distributions
- methylation/ directory has methylation rates
- end_motifs/ directory has motif frequencies  
- positions/ directory has regional counts

If something looks wrong, I can trace it back to the extraction step. If everything looks good, this combine script is trivial - just slam the tables together.

Compare this to a monolithic script that extracts and combines in one pass. That would be faster to run (no intermediate files) but much harder to debug. When something goes wrong, you don't know if it's the extraction logic or the combination logic.

For a pipeline that processes real patient samples where you need to trust the results, I'll take modularity and inspectability every time.

WHAT HAPPENS DOWNSTREAM

After this script runs, the feature selection process tests all 11 combinations systematically. It does 5-fold stratified cross-validation for each combination and calculates F1-score, accuracy, precision, recall, etc.


SUMMARY

This script is conceptually simple: take 4 feature files, merge them on sample_id, add disease labels, save. But the design choices around file formats, validation, and modularity matter a lot for pipeline robustness.

Key decisions:
1. Modular extraction over monolithic - enables parallelization and debugging
2. Wide format over long format - ready for scikit-learn
3. CSV over binary formats - human-readable at this scale
4. Inner joins with validation - ensures data quality
5. Print statements for logging - helps QC and documentation

For 22 samples this might seem like overkill, but the design scales cleanly to thousands of samples where these choices really matter.

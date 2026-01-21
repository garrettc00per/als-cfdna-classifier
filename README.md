# cfDNA Fragmentomics Analysis Pipeline

Nextflow pipeline for multi-feature analysis of cell-free DNA with machine learning-based disease classification.

## Pipeline Overview
```
BAMs → [1] Extract Features → [2] Select Features → [3] Classify → [4] Visualize
       ✓ Complete             ✓ Complete           ✓ Complete   ✓ Complete
```

**Implemented:**

✅ Insert Size Distribution Analysis  
✅ End Motif Profiling (4-mer frequencies)  
✅ Methylation Pattern Extraction (XM tags)  
✅ Positional Distribution Analysis  
✅ Automated Feature Selection (11 combinations tested)  
✅ Binary Classification (ALS vs Control)  
✅ Publication-Quality Visualizations

---

## Quick Start
```bash
# Clone repository
git clone https://github.com/garrettc00per/als-cfdna-classifier.git
cd als_cfdna_classifier

# Install dependencies
mamba install -c bioconda samtools nextflow
mamba install pandas scikit-learn matplotlib seaborn scipy

#Install raw data
#Can be found here: https://drive.google.com/file/d/1Xt4SPbL8sYa9WAnUP_j0rZ_T9AJl7oTr/view

# Run full pipeline
nextflow run main.nf

# Results in results/ directory
ls results/
```

---

## Setup & Installation (Amazon Linux 2)

If you're setting up on a fresh AWS EC2 instance, follow these steps:

### 1. Install Git
```bash
sudo yum install -y git
git --version  # Verify
```

### 2. Install Java 17 (Required)
```bash
sudo yum install java-17-amazon-corretto -y
java -version  # Should show Java 17
```

### 3. Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
sudo mv ~/nextflow /usr/local/bin/
```

### 4. Install Miniforge
```bash
cd ~
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
source ~/.bashrc
```

### 5. Add bioconda channel and install tools
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

mamba install samtools -y


mamba install pandas scikit-learn matplotlib seaborn scipy
```

### 6. Download Raw Data
```bash
#Data can be found here: https://drive.google.com/file/d/1Xt4SPbL8sYa9WAnUP_j0rZ_T9AJl7oTr/view
#Download to your computer, upload to S3 and pull to your EC2 instance
```


### 7. Clone and Run
```bash
git clone https://github.com/garrettc00per/als-cfdna-classifier.git
cd als-cfdna-classifier
nextflow run main.nf
```

---

## Summary of Results (20 Replicates)

**Important Note on Reproducibility:** Due to the small sample size (n=22) and stochastic nature of cross-validation, results will vary between pipeline runs even with the same data. The rankings and absolute performance values below represent one execution; repeat runs may show different feature combinations as "best." This variability is expected and scientifically meaningful. It reflects the fundamental challenge of feature selection in low-sample datasets.

**Top Performing Combinations (Mean ± SD across 20 replicates):**
1. **Insert + Positions** (11 features): 70.2% F1 ± 6.2%
2. **Motifs Only** (40 features): 70.1% F1 ± 8.1%
3. **Insert + Motifs** (48 features): 69.3% F1 ± 7.0%
4. **Positions Only** (3 features): 68.2% F1 ± 6.0%

**Best Classifier Performance (Insert + Positions with 11 features):**
- **F1-Score: 69.4% ± 8.3%**
- **Accuracy: 68.6% ± 6.6%**
- **Sensitivity: 66.7% ± 12.1%**
- **Precision: 73.9% ± 7.0%**
- **Model: Logistic Regression**

**Key Insight:** Top 4 combinations are highly competitive (all ~68-70% F1), with overlapping confidence intervals. This demonstrates that **multiple cfDNA characteristics are informative for ALS detection**, and no single feature set dominates with n=22 samples. Different runs may rank these combinations differently. This is the expected behavior of feature selection with limited data, not a pipeline failure.

---

## Feature Extraction Modules

### 1. Insert Size Analysis (✅ Complete)
- Fragment length distribution from paired-end reads
- Statistical summaries (mean, median, stddev, quartiles)
- Per-sample fragment counts
- **Output:** `results/insert_sizes/`

### 2. End Motif Profiling (✅ Complete)
- 4-mer frequencies at fragment starts and ends
- Top 20 start motifs + top 20 end motifs
- Normalized by total fragment count
- Accounts for bisulfite treatment (C→T conversion)
- **Output:** `results/end_motifs/`

### 3. Methylation Extraction (✅ Complete)
- Parse XM tags from Bismark alignment
- CpG, CHG, CHH methylation rates
- Context-specific methylation counts
- Overall methylation statistics
- **Output:** `results/methylation/`

### 4. Positional Distribution (✅ Complete)
- Fragment start positions across chr21
- Regional density (3 genomic regions)
- 1Mb binned coverage tracks
- **Output:** `results/positions/`

### 5. Feature Selection (✅ Complete)
- Tests 11 feature combinations systematically
- Cross-validated performance (5-fold stratified)
- Automatic selection of optimal feature set
- **Output:** `results/feature_selection_results.csv`

### 6. Classification (✅ Complete)
- Logistic Regression + Random Forest models
- Stratified K-fold cross-validation
- Comprehensive performance metrics
- Confusion matrices and classification reports
- **Output:** `results/classification_*.csv`

### 7. Visualizations (✅ Complete)
- Feature importance plots
- Distribution comparisons (ALS vs Control)
- Confusion matrices
- Methylation and insert size boxplots
- Positional distribution plots
- Pipeline summary diagram
- **Output:** `results/plots/`

---

## Input

### BAM Files
Pre-processed bisulfite-treated cfDNA sequencing data:
- Aligned with Bismark aligner (provides XM methylation tags)
- Paired-end sequencing
- Chromosome 21 only (downsampled to ~10M reads/sample)
- Indexed (.bam.bai files present)

**Location:** `/home/ec2-user/celfie/bam/`  
**Naming:** `{SRRXXXXXX}.deduplicated.sorted_ds10mill_chr21.bam`

**Dataset:**
- 22 samples total
- 12 ALS patients
- 10 Healthy controls

### Metadata File
`celfie_cfDNA_ss.csv` - Sample metadata with disease status

Required columns:
- `Run` - Sample ID (SRR accession)
- `disease_status` - "als" or "ctrl"
- `AGE`, `tissue`, etc. (optional)

---

## Output
```
results/
├── insert_sizes/                          # Per-sample fragment lengths
│   └── {sample}_insert_sizes.txt
├── end_motifs/                            # Per-sample motif counts
│   └── {sample}_end_motifs.txt
├── methylation/                           # Per-sample methylation stats
│   └── {sample}_methylation.txt
├── positions/                             # Per-sample position distributions
│   ├── {sample}_position_stats.txt
│   └── {sample}_position_bins.txt
├── plots/                                 # All visualizations
│   ├── visualizations_feature_importance.png
│   ├── visualizations_confusion_matrix.png
│   ├── visualizations_methylation_comparison.png
│   ├── visualizations_insert_size_comparison.png
│   ├── visualizations_position_distributions.png
│   └── pipeline_summary.png
├── all_features_combined.csv              # Complete feature matrix
├── feature_selection_results.csv          # Performance of 11 combinations (20 replicates)
├── best_features.json                     # Optimal feature set
├── classification_results.csv             # Final model performance (20 replicates)
├── classification_output.txt              # Detailed metrics
├── insert_size_summary_labeled.txt        # Summary statistics
├── methylation_summary.txt                # Methylation rates
└── position_summary.txt                   # Regional distributions
```

---

## Performance

### Benchmarking
Tested on **AWS EC2 t2.xlarge** (4 vCPUs, 16 GB RAM):

| Task | Replicates | Wall Time | Notes |
|------|------------|-----------|-------|
| Full Feature Extraction | - | ~45 sec | 22 samples, parallel |
| Feature Selection | 20 | ~2 min 30 sec | 11 combinations × 20 replicates |
| Classification | 20 | ~45 sec | Final model training |
| Visualization | - | ~30 sec | 7 publication plots |
| **Full Pipeline** | **20 replicates** | **~3 min 45 sec** | **End-to-end** |

*Note: 5-replicate runs complete in ~2 min; 20-replicate runs recommended for robust feature selection*

### Resource Usage
- **CPU:** 2 cores per sample (parallelized via Nextflow)
- **Memory:** ~4 GB per sample
- **Disk:** ~500 MB for all outputs
- **Parallelization:** Up to 22 samples simultaneously

---

## Requirements

### Software
- Nextflow (v25.10.2+)
- Samtools (v1.23+)
- Python (v3.12+)

### Python Packages
```
pandas>=2.3.3
scikit-learn>=1.8.0
matplotlib>=3.10.8
seaborn>=0.13.2
scipy>=1.17.0
numpy>=2.4.1
```

Install via conda/mamba:
```bash
mamba install -c bioconda samtools nextflow
mamba install pandas scikit-learn matplotlib seaborn scipy
```

---

## Parameters

### Main Pipeline
- `--bam_dir` - BAM file directory (default: `/home/ec2-user/celfie/bam`)
- `--outdir` - Output directory (default: `results`)
- `--metadata` - Sample metadata file (default: `/home/ec2-user/celfie/celfie_cfDNA_ss.csv`)

### Nextflow Options
```bash
# Resume failed pipeline (skips completed steps)
nextflow run main.nf -resume

# Generate execution report
nextflow run main.nf -with-report report.html

# Profile execution timeline
nextflow run main.nf -with-timeline timeline.html
```

---

## Key Findings

### 1. Multiple Feature Combinations Show Competitive Performance
- Top 4 feature sets achieve 68-70% F1-score with overlapping error bars
- **Insert + Positions** (11 features), **Motifs Only** (40 features), and **Insert + Motifs** (48 features) perform nearly identically
- **Run-to-run variability is expected**: Different executions may rank these combinations differently due to CV fold randomization with n=22
- **Key Insight:** This is NOT a failure—it demonstrates that **with limited samples, multiple feature solutions are equally valid**
- **Biological Interpretation:** cfDNA changes in ALS manifest across multiple molecular scales (fragmentation patterns, sequence composition, genomic positioning)

### 2. Insert Size + Positional Features Are Most Informative (This Run)
- **Insert + Positions (11 features):** 70.2% F1 ± 6.2% (current best combination)
- Combining fragment size statistics with genomic position distributions improves discrimination
- Suggests ALS-specific fragmentation patterns are both quantitative (fragment length) and spatial (genomic localization)
- **However:** Motifs Only (70.1% F1) performs nearly identically—difference is within noise

### 3. End Motifs Show Surprisingly Strong Individual Performance  
- **Motifs Only**: 70.1% F1 ± 8.1% (competitive with best combination)
- 4-mer frequencies at fragment ends capture sequence-specific fragmentation preferences
- Likely reflects altered chromatin accessibility and nucleosome positioning in ALS
- High standalone performance suggests motif analysis deserves further investigation

### 4. Methylation Shows Weak Standalone Performance
- **Methylation Only**: 60.4% F1 ± 7.5% (lowest individual performer)
- Does not significantly improve performance when added to other features
- May provide complementary information in larger cohorts
- CpG/CHG/CHH patterns alone insufficient for ALS discrimination with chr21 data

### 5. Feature Selection Instability Reflects Small Sample Size
- Top 4 combinations within ~2% F1 (within statistical noise)
- **Rankings will change between runs** due to CV fold randomization
- This is **scientifically valid behavior**, not a bug
- **Biological Interpretation:** ALS involves **multifactorial cfDNA changes** rather than a single dominant biomarker
- **Recommendation:** Larger cohorts (n>100) needed for stable feature selection and clinical validation

---

## Key Features

### Automated Feature Selection
- Systematically tests all feature combinations
- Runs 20 replicates for robustness assessment
- Calculates cross-validated performance metrics
- Prevents manual bias in feature selection

### Modular Design
- Independent, reusable Nextflow processes
- Easy to customize and extend
- Can run individual modules separately
- Resume capability for failed runs

### Production Ready
- Comprehensive error handling
- Detailed classification reports
- Publication-quality visualizations
- Resource-aware parallelization

### Scientifically Rigorous
- Stratified K-fold cross-validation (5 folds × 20 replicates)
- Multiple classifiers tested (RF + Logistic Regression)
- Feature standardization (zero mean, unit variance)
- Proper handling of class imbalance

---

## Future Improvements

Given more time and resources, valuable extensions include:

### 1. Scalability Enhancements (High Priority)
- Support for full-genome analysis (beyond chr21)
- Streaming parsers for 50GB+ BAM files
- Distributed computing for large cohorts
- GPU acceleration for deep learning models

### 2. Advanced Feature Engineering
- **Regional methylation patterns:**
  - CpG island methylation
  - Promoter/enhancer-specific patterns
  - Differentially methylated regions (DMRs)
- **Fragment coverage profiles:**
  - Nucleosome positioning signals
  - Transcription factor binding site enrichment
- **Advanced motif analysis:**
  - k-mer analysis beyond 4-mers
  - Motif enrichment vs genomic background

### 3. Statistical Rigor
- Power analysis for sample size determination
- Multiple testing correction (Bonferroni/FDR)
- Permutation testing for significance
- External validation on independent cohort

### 4. Production Deployment
- Containerization (Docker/Singularity)
- Interactive dashboard (Shiny/Plotly)
- Automated PDF report generation

---

## Technical Details

### Software Versions
- Nextflow: v25.10.2
- Samtools: v1.23
- Python: v3.12.12
- scikit-learn: v1.8.0
- pandas: v2.3.3
- matplotlib: v3.10.8
- seaborn: v0.13.2
- scipy: v1.17.0
- numpy: v2.4.1

### Dataset Information
- **Source:** Cell-free DNA from ALS study: https://www.nature.com/articles/s41467-021-22901-x
- **Sequencing:** Paired-end, bisulfite-treated
- **Aligner:** Bismark (methylation-aware)
- **Genome:** hg38 (Chromosome 21 only)
- **Read depth:** ~10M reads per sample (downsampled)

### Model Details
**Best Performer (Insert + Positions, 11 features):**
- Algorithm: Random Forest (100 estimators, max_depth=5)
- Cross-validation: 5-fold stratified × 20 replicates
- Performance: 69.8% F1 ± 8.4%, 69.3% Accuracy ± 8.7%

**Recommendation:** Given the competitive performance of top 4 combinations, ensemble approaches or larger validation cohorts recommended for clinical deployment.

---

## Troubleshooting

### Common Issues

**Pipeline fails at feature extraction:**
```bash
# Check if samtools is installed
samtools --version

# Verify BAM files exist
ls /home/ec2-user/celfie/bam/*.bam
```

**Out of memory errors:**
```bash
# Reduce parallelization in nextflow.config
process.cpus = 1
process.memory = '2 GB'
```

**Python package not found:**
```bash
# Install all required packages
mamba install pandas scikit-learn matplotlib seaborn scipy
```

**Nextflow resumes from wrong point:**
```bash
# Clean work directory and start fresh
rm -rf work/ .nextflow*
nextflow run main.nf
```

---

## Project Structure
```
celfie-analysis/
├── main.nf                          # Main Nextflow workflow
├── nextflow.config                  # Pipeline configuration
├── README.md                        # This file
├── requirements.txt                 # Python dependencies
├── bin/                             # Analysis scripts
│   ├── extract_end_motifs.sh
│   ├── extract_methylation.sh
│   ├── extract_positions.sh
│   ├── compute_insert_stats.sh
│   ├── create_motif_features.py
│   ├── combine_all_features.py
│   ├── feature_selection.py
│   ├── classify_best_features.py
│   ├── create_visualizations.py
│   ├── plot_position_distributions.py
│   └── create_summary_figure.py
└── results/                         # Output directory
    ├── plots/
    ├── insert_sizes/
    ├── end_motifs/
    ├── methylation/
    ├── positions/
    └── *.csv, *.txt
```

---

## Author

Garrett Cooper, PhD
Genetics and Molecular Biology  
Department of Pediatrics

**Time Invested:** ~6 hours (implementation + visualization + documentation)  
**Lines of Code:** ~1,200 across all scripts  
**Pipeline Processes:** 14 automated processes

---

## References

### Data Source
Cell-free DNA sequencing data from ALS cohort study: https://www.nature.com/articles/s41467-021-22901-x

### Key Methods
- **Bismark aligner:** Methylation-aware alignment (Krueger & Andrews, 2011)
- **Random Forest:** Ensemble learning (Breiman, 2001)
- **Feature Selection:** Cross-validated comparison (Guyon & Elisseeff, 2003)

### Tools & Frameworks
- **Nextflow:** Workflow management (Di Tommaso et al., 2017)
- **Samtools:** BAM file processing (Li et al., 2009)
- **Scikit-learn:** Machine learning (Pedregosa et al., 2011)

---

## License

This code is provided for evaluation purposes as part of a take-home assignment.

---

**Last Updated:** January 2026

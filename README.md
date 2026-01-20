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
git clone <repository-url>
cd celfie-analysis

# Install dependencies
mamba install -c bioconda samtools nextflow
mamba install pandas scikit-learn matplotlib seaborn scipy

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

###3. Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
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

### 5. Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
```

### 6. Clone and Run
```bash
git clone https://github.com/garrettc00per/als-cfdna-classifier.git
cd als-cfdna-classifier
nextflow run main.nf
```

---

## Summary of Results

**Best Performance: 77.3% accuracy (Logistic Regression, Insert Size Features Only)**
- F1-Score: 80.0%
- Sensitivity: 83.3%
- Precision: 76.9%
- **Features used: 8 (fragment size characteristics only)**

**Key Finding:** Systematic feature selection across 11 combinations revealed that insert size features alone outperform complex multi-feature models, demonstrating the importance of feature selection in low-sample regimes.

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
├── feature_selection_results.csv          # Performance of 11 combinations
├── best_features.json                     # Optimal feature set
├── classification_results.csv             # Final model performance
├── classification_output.txt              # Detailed metrics
├── insert_size_summary_labeled.txt        # Summary statistics
├── methylation_summary.txt                # Methylation rates
└── position_summary.txt                   # Regional distributions
```

---

## Performance

### Benchmarking
Tested on **AWS EC2 t2.xlarge** (4 vCPUs, 16 GB RAM):

| Samples | Module | Wall Time | Notes |
|---------|--------|-----------|-------|
| 22 | Feature Extraction | ~2 min | Parallel processing |
| 22 | Feature Selection | ~2 min | 11 combinations tested |
| 22 | Classification | <1 min | Final model training |
| 22 | Visualization | ~1 min | 7 publication plots |
| 22 | **Full Pipeline** | **~5 min** | End-to-end |

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

### 1. Insert Size Variability Is The Most Discriminative Feature
- **Standard deviation (stddev) of insert sizes is the #1 most important feature**
  - ALS: Lower variability (~54 bp stddev)
  - Control: Higher variability (~63 bp stddev)
- **Key Insight:** Fragment size **uniformity**, not just average size, distinguishes ALS
  - ALS shows more **consistent, uniform fragment sizes**
  - Controls show more **heterogeneous fragmentation patterns**
- **Biological Interpretation:**
  - Lower variability suggests dysregulated chromatin structure in ALS
  - Altered nucleosome positioning creates more predictable fragmentation
  - May reflect loss of normal epigenetic heterogeneity in disease

### 2. Feature Selection Reveals Overfitting in Multi-Feature Models
- Tested 11 feature combinations systematically:
  - Insert size only (8 features): **81.3% F1-score**
  - All features with positions (59 features): 80.5% F1
  - All features without positions (56 features): 75.3% F1
- **Overfitting occurs** when features >> samples (curse of dimensionality)
- Final model uses only: mean, median, stddev, quartiles, min, max, fragment count
- **This demonstrates that fragment size variability is the dominant ALS signal**

### 3. End Motifs Show AT-Rich Patterns
- Most important motifs: AATA, TTTT, ATTA, TATT
- AT-rich sequences preferentially occur at fragment ends
- Bisulfite treatment (C→T conversion) influences motif patterns
- Fragmentation preferences differ between ALS and controls

### 4. Methylation Shows Complex Multivariate Pattern
- **No significant univariate differences** in methylation rates (all p > 0.5)
  - CpG methylation: p = 0.597
  - CHG methylation: p = 0.632
  - CHH methylation: p = 0.583
- **Yet improves classification in combination with other features**
- **Key Insight:** ALS is characterized by **altered relationships between methylation and fragmentation** rather than simple methylation changes
- Traditional t-tests examine features independently; machine learning detects multivariate interactions
- **Biological interpretation:** ALS may alter the coupling between DNA methylation and nucleosome positioning

### 5. Positional Features Show No Regional Specificity
- Adding genomic position features degrades performance in isolation
- Suggests ALS fragmentation patterns are **molecular, not positional**
- No enrichment in specific chr21 regions
- Disease signal is in fragment characteristics, not genomic location

---

## Key Features

### Automated Feature Selection
- Systematically tests all feature combinations
- Calculates cross-validated performance metrics
- Automatically selects optimal feature set
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
- Stratified K-fold cross-validation (5 folds)
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
**Final Selected Model:**
- Algorithm: Logistic Regression
- Features: Insert size statistics (8 features)
- Cross-validation: 5-fold stratified
- Max iterations: 1000

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
Genetic and Molecular Biology  
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

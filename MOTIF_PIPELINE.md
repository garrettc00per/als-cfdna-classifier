# Motif Feature Extraction Pipeline

## Overview

This pipeline extracts DNA sequence motifs from cell-free DNA fragments and converts them into normalized features for machine learning classification.

**Biological Rationale:** Cell-free DNA fragment ends have non-random sequence patterns influenced by nucleosome positioning, enzymatic cleavage preferences, and DNA accessibility. These patterns can serve as biomarkers for disease states like ALS.

---

## Pipeline Architecture
```
BAM files (aligned reads)
    ↓
[extract_end_motifs.sh]
    ↓
*_end_motifs.txt (raw motif counts)
    ↓
[create_motif_features.py]
    ↓
motif_features.csv (normalized ML features)
    ↓
[classification models]
```

---

## Step 1: Extract End Motifs (`extract_end_motifs.sh`)

### Purpose
Extract the first 4 bases from fragment START (5' end) and END (3' end) positions, then count occurrences of each 4-base motif.

### Input
- BAM file with paired-end reads
- Sample ID

### Output
Text file with format:
```
=== START MOTIFS ===
3336 AAAA
2100 GCTA
1500 TTTT
=== END MOTIFS ===
1500 TTTT
1200 CGCG
```

### How It Works

#### 1. Extract Properly Paired Reads
```bash
samtools view -f 0x2 "$bam_file"
```
- `-f 0x2` flag filters for properly paired reads only
- Ensures we have both ends of each fragment

#### 2. Extract Motifs with AWK
```bash
awk '{
    if ($9 > 0) {  # First in pair (READ1)
        seq = $10
        if (length(seq) >= 4) {
            print "START", substr(seq, 1, 4)
        }
    }
    if ($9 < 0) {  # Second in pair (READ2)
        seq = $10
        if (length(seq) >= 4) {
            print "END", substr(seq, 1, 4)
        }
    }
}'
```

**Key points:**
- `$9` (TLEN/template length) sign indicates read pair order:
  - Positive = READ1 (leftmost read, sequences fragment start)
  - Negative = READ2 (rightmost read, sequences fragment end)
- `$10` contains the DNA sequence
- `substr(seq, 1, 4)` extracts first 4 bases
- READ2 is stored reverse-complemented in BAM, so its first 4bp already represent the fragment end

**Why first 4bp of both reads?**
- READ1 sequences from 5' end → first 4bp = fragment START
- READ2 is reverse-complemented in BAM → first 4bp = fragment END

#### 3. Count and Sort Motifs
```bash
grep "^START" ${sample_id}_motifs_raw.txt | \
    awk '{print $2}' | \
    sort | uniq -c | sort -rn
```

**Pipeline breakdown:**
- `grep "^START"` - Get only START lines
- `awk '{print $2}'` - Extract just the motif (second column)
- `sort` - Group identical motifs together
- `uniq -c` - Count occurrences
- `sort -rn` - Sort by count (highest first)

---

## Step 2: Create Feature Matrix (`create_motif_features.py`)

### Purpose
Convert raw motif counts from all samples into normalized frequency features suitable for machine learning.

### Input
- Directory containing `*_end_motifs.txt` files
- Output CSV filename

### Output
CSV with normalized motif frequencies:
```csv
sample_id,start_AAAA,start_TTTT,...,end_TTTT,end_AAAA,...
SRR13404367,0.3336,0.1500,...,0.1875,0.1250,...
SRR13404368,0.3200,0.1600,...,0.1800,0.1300,...
```

---

### Function 1: `parse_motif_file(filepath)`

**Purpose:** Parse a single motif count file and return structured data

**Input:** Path to `*_end_motifs.txt` file

**Output:** Dictionary with structure:
```python
{
    'START': {'AAAA': 3336, 'GCTA': 2100, 'TTTT': 1500, ...},
    'END': {'TTTT': 1500, 'CGCG': 1200, 'AAAA': 800, ...}
}
```

**How it works:**

1. **Initialize storage:**
```python
motifs = {'START': {}, 'END': {}}
current_section = None
```
- Creates separate dictionaries for START and END motifs
- `current_section` tracks which section we're currently reading

2. **Read file line by line:**
```python
with open(filepath, 'r') as f:
    for line in f:
        line = line.strip()
```
- Opens file in read mode
- Strips whitespace from each line

3. **Identify section headers:**
```python
if '=== START MOTIFS ===' in line:
    current_section = 'START'
elif '=== END MOTIFS ===' in line:
    current_section = 'END'
```
- Sets `current_section` when encountering headers
- Subsequent lines are parsed into that section

4. **Parse data lines:**
```python
elif line and current_section:
    parts = line.split()  # "3336 AAAA" → ["3336", "AAAA"]
    if len(parts) == 2:
        count = int(parts[0])  # 3336
        motif = parts[1]       # "AAAA"
        motifs[current_section][motif] = count
```
- Only processes non-empty lines within a section
- Splits line by whitespace
- Validates format (must have exactly 2 parts)
- Stores count with motif as key

5. **Return parsed data:**
```python
return motifs
```

**Key Design Decisions:**
- **Separate START/END storage:** These represent different biological signals (5' vs 3' fragmentation patterns)
- **Dictionary structure:** Allows O(1) lookup by motif name
- **Robust parsing:** Validates line format and handles malformed data

---

### Function 2: `main(motif_dir, output_file)`

**Purpose:** Orchestrate the entire feature extraction pipeline

**Process Overview:**

#### Phase 1: Data Loading

**Find all motif files:**
```python
motif_files = [f for f in os.listdir(motif_dir) if f.endswith('_end_motifs.txt')]
```

**Initialize storage:**
```python
all_start_motifs = set()  # Collect unique START motifs
all_end_motifs = set()    # Collect unique END motifs
sample_data = {}          # Store parsed data for each sample
```

**Parse all files:**
```python
for filename in motif_files:
    sample_id = filename.replace('_end_motifs.txt', '')
    filepath = os.path.join(motif_dir, filename)
    
    motifs = parse_motif_file(filepath)  # Call parsing function
    sample_data[sample_id] = motifs      # Store results
    
    # Collect all unique motifs across samples
    all_start_motifs.update(motifs['START'].keys())
    all_end_motifs.update(motifs['END'].keys())
```

**Result:**
```python
sample_data = {
    'SRR13404367': {
        'START': {'AAAA': 3336, 'GCTA': 2100, ...},
        'END': {'TTTT': 1500, 'CGCG': 1200, ...}
    },
    'SRR13404368': {...},
    ...
}
```

#### Phase 2: Feature Selection

**Why feature selection?**
- Datasets may have 100+ unique motifs
- Many are rare and uninformative
- Too many features → overfitting
- Solution: Select top 20 most common motifs

**Sum motif counts across all samples:**
```python
start_totals = defaultdict(int)
end_totals = defaultdict(int)

for sample_id, motifs in sample_data.items():
    for motif, count in motifs['START'].items():
        start_totals[motif] += count  # Accumulate counts
    for motif, count in motifs['END'].items():
        end_totals[motif] += count
```

**Example accumulation:**
```python
# Sample 1: {'AAAA': 3336, 'GCTA': 2100}
# Sample 2: {'AAAA': 3200, 'GCTA': 2050}
# Sample 3: {'AAAA': 3100, 'TTTT': 1500}

# Result:
start_totals = {
    'AAAA': 3336 + 3200 + 3100 = 9636,
    'GCTA': 2100 + 2050 = 4150,
    'TTTT': 1500
}
```

**Select top 20 motifs:**
```python
# Sort by count (descending) and take first 20
top_start = sorted(start_totals.items(), key=lambda x: x[1], reverse=True)[:20]
top_end = sorted(end_totals.items(), key=lambda x: x[1], reverse=True)[:20]

# Extract just motif names
top_start_motifs = [m for m, c in top_start]
top_end_motifs = [m for m, c in top_end]
```

**Result:**
```python
top_start_motifs = ['AAAA', 'TTTT', 'GCTA', 'AAAT', ...]  # 20 motif names
top_end_motifs = ['TTTT', 'AAAA', 'CGCG', 'TATA', ...]    # 20 motif names
```

#### Phase 3: Feature Normalization

**Why normalize?**
- Different samples have different sequencing depths
- Sample 1: 10,000 total fragments
- Sample 2: 15,000 total fragments
- Raw counts aren't comparable
- Solution: Convert to frequencies (proportions)

**Build feature table:**
```python
rows = []
for sample_id in sorted(sample_data.keys()):
    motifs = sample_data[sample_id]
    
    # Calculate totals for normalization
    total_start = sum(motifs['START'].values())
    total_end = sum(motifs['END'].values())
    
    row = {'sample_id': sample_id}
    
    # Add normalized START motif frequencies
    for motif in top_start_motifs:
        count = motifs['START'].get(motif, 0)  # 0 if motif not present
        freq = count / total_start if total_start > 0 else 0
        row[f'start_{motif}'] = freq
    
    # Add normalized END motif frequencies
    for motif in top_end_motifs:
        count = motifs['END'].get(motif, 0)
        freq = count / total_end if total_end > 0 else 0
        row[f'end_{motif}'] = freq
    
    rows.append(row)
```

**Example normalization:**
```python
# Sample has 10,000 total START motifs
# AAAA appears 3336 times
# Frequency = 3336 / 10000 = 0.3336 (33.36%)

row = {
    'sample_id': 'SRR13404367',
    'start_AAAA': 0.3336,  # 33.36% of fragments start with AAAA
    'start_TTTT': 0.1500,  # 15.00%
    'start_GCTA': 0.2100,  # 21.00%
    ...
    'end_TTTT': 0.1875,    # 18.75% of fragments end with TTTT
    'end_AAAA': 0.1250,    # 12.50%
    ...
}
```

#### Phase 4: Output Generation

**Create DataFrame:**
```python
df = pd.DataFrame(rows)
```

**Save to CSV:**
```python
df.to_csv(output_file, index=False)
```

**Final output:**
```csv
sample_id,start_AAAA,start_TTTT,start_GCTA,...,end_TTTT,end_AAAA,end_CGCG,...
SRR13404367,0.3336,0.1500,0.2100,...,0.1875,0.1250,0.1500,...
SRR13404368,0.3200,0.1600,0.2050,...,0.1800,0.1300,0.1450,...
SRR13404369,0.3100,0.1550,0.2000,...,0.1850,0.1200,0.1400,...
```

---

## Key Design Principles

### 1. Separation of START and END Motifs
**Why:** Fragment ends represent different biological processes
- START: Reflects where DNA fragmentation initiated
- END: Reflects where fragmentation terminated
- These patterns are influenced by different factors (nucleosome positioning, enzymatic preferences)
- Combining them would lose information

**Example:**
```python
# Healthy sample
start_AAAA = 0.35, end_AAAA = 0.15

# Disease sample  
start_AAAA = 0.35, end_AAAA = 0.30  # Increased at END!

# Model can learn: "High end_AAAA = disease signature"
```

### 2. Feature Selection (Top 20)
**Why:** Balance information vs. overfitting
- More features = more information (good)
- More features = more overfitting risk (bad)
- Top 20 captures most signal while keeping model manageable
- Rare motifs (< top 20) are typically noise

### 3. Frequency Normalization
**Why:** Make samples comparable
- Sequencing depth varies between samples
- Raw counts reflect both biology and technical variation
- Frequencies represent true proportions
- Enables fair comparison across samples

**Formula:**
```
frequency = motif_count / total_motif_count
```

**Example:**
```
Sample 1: AAAA count = 3336, total = 10000 → freq = 0.3336
Sample 2: AAAA count = 6000, total = 15000 → freq = 0.4000

Sample 2 has higher frequency despite different depths!
```

---

## Usage

### Command Line

**Extract motifs from BAM:**
```bash
./bin/extract_end_motifs.sh input.bam sample_001
```

**Create feature matrix:**
```bash
python bin/create_motif_features.py results/end_motifs/ motif_features.csv
```

### Nextflow Pipeline

The pipeline is automated using Nextflow:
```groovy
// Extract motifs from all BAM files
EXTRACT_END_MOTIFS(bam_channel)

// Collect all motif files
all_motifs = EXTRACT_END_MOTIFS.out.collect()

// Create feature matrix
CREATE_MOTIF_FEATURES(all_motifs)
```

---

## Output Interpretation

### Feature Names
- `start_AAAA`: Frequency of AAAA motif at fragment starts (5' end)
- `end_TTTT`: Frequency of TTTT motif at fragment ends (3' end)

### Values
- Range: 0.0 to 1.0 (proportions)
- Interpretation: 0.3336 = 33.36% of fragments have this motif
- Sum across all motifs in a category ≈ 1.0 (accounts for top 20 only)

### Biological Significance
- Poly-A/Poly-T motifs typically most common (AT-rich regions)
- Disease may show altered fragmentation patterns
- Different motif frequencies at START vs END indicate directional preferences

---

## Implementation Notes

### Python Libraries Used
- `pandas`: DataFrame manipulation and CSV I/O
- `numpy`: Statistical calculations
- `collections.defaultdict`: Automatic zero initialization for counting

### Bash Tools Used
- `samtools`: BAM file processing
- `awk`: Pattern matching and text extraction
- `grep`, `sort`, `uniq`: Text processing pipeline

### Performance Considerations
- File I/O: Reads each BAM once, writes text output
- Memory: Stores motif counts for all samples in memory (manageable for typical datasets)
- Runtime: Dominated by BAM processing (samtools), feature extraction is fast

---

## Future Enhancements

**Potential improvements:**
1. **Configurable feature count:** Allow user to specify number of top motifs (currently hardcoded at 20)
2. **Motif length flexibility:** Support variable-length motifs (currently fixed at 4bp)
3. **Position-specific features:** Capture motif positions along fragments
4. **Dinucleotide frequencies:** Add GC content, CpG frequency as features
5. **Quality filtering:** Filter motifs by base quality scores

---

## References

- Cell-free DNA fragmentation patterns: doi:10.1038/nature23879
- Nucleosome positioning and cfDNA: doi:10.1038/s41586-019-1426-6

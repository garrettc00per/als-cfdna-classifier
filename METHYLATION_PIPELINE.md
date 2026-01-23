METHYLATION EXTRACTION CODE WALKTHROUGH

Okay so let me walk you through how this methylation extraction works. It's actually pretty straightforward once you understand what Bismark puts in the XM tag.

============================================================
THE BASIC IDEA
============================================================

Bismark adds this XM tag to every read that encodes methylation state at every cytosine position. The tag uses different letters for methylated vs unmethylated cytosines in different sequence contexts. So instead of storing raw methylation calls, we can just parse these strings and count up the different letter types.


============================================================
THE CODE
============================================================

samtools view -f 0x2 "$bam_file" | \
awk '{
    for(i=12; i<=NF; i++) {
        if($i ~ /^XM:Z:/) {
            xm = substr($i, 6)
            
            cpg_meth = gsub(/Z/, "", xm)      # Methylated CpG
            cpg_unmeth = gsub(/z/, "", xm)    # Unmethylated CpG
            chg_meth = gsub(/X/, "", xm)      # Methylated CHG
            chg_unmeth = gsub(/x/, "", xm)    # Unmethylated CHG
            chh_meth = gsub(/H/, "", xm)      # Methylated CHH
            chh_unmeth = gsub(/h/, "", xm)    # Unmethylated CHH
            
            total_cpg_meth += cpg_meth
            total_cpg_unmeth += cpg_unmeth
            total_chg_meth += chg_meth
            total_chg_unmeth += chg_unmeth
            total_chh_meth += chh_meth
            total_chh_unmeth += chh_unmeth
            n_reads++
        }
    }
}
END {
    total_cpg = total_cpg_meth + total_cpg_unmeth
    total_chg = total_chg_meth + total_chg_unmeth
    total_chh = total_chh_meth + total_chh_unmeth
    total_c = total_cpg + total_chg + total_chh
    
    cpg_rate = (total_cpg > 0) ? total_cpg_meth / total_cpg : 0
    chg_rate = (total_chg > 0) ? total_chg_meth / total_chg : 0
    chh_rate = (total_chh > 0) ? total_chh_meth / total_chh : 0
    overall_rate = (total_c > 0) ? (total_cpg_meth + total_chg_meth + total_chh_meth) / total_c : 0
    
    print "'$sample_id'", cpg_rate, chg_rate, chh_rate, overall_rate, total_cpg, total_chg, total_chh, n_reads
}'


============================================================
STEP 1: FILTER READS
============================================================

samtools view -f 0x2 "$bam_file"

This just grabs properly paired reads. Flag 0x2 means the read mapped in a proper pair - both mates aligned at the expected distance and orientation. I'm filtering for these because they're higher quality and less likely to have alignment artifacts that mess up methylation calls.


============================================================
STEP 2: FIND THE XM TAG
============================================================

for(i=12; i<=NF; i++) {
    if($i ~ /^XM:Z:/) {
        xm = substr($i, 6)

BAM format has 11 mandatory columns, then optional tags start at field 12. So I loop through everything from field 12 onwards looking for anything that starts with "XM:Z:". When I find it, substr($i, 6) just strips off the "XM:Z:" prefix and gives me the actual methylation string.

Here's what that looks like with a real read:

XM:Z:hh.h.....Zx........h........h...h............x.......h.............x.......x......h..........

After substr, I get:
hh.h.....Zx........h........h...h............x.......h.............x.......x......h..........


============================================================
STEP 3: COUNT METHYLATION STATES
============================================================

cpg_meth = gsub(/Z/, "", xm)      # Methylated CpG
cpg_unmeth = gsub(/z/, "", xm)    # Unmethylated CpG
chg_meth = gsub(/X/, "", xm)      # Methylated CHG
chg_unmeth = gsub(/x/, "", xm)    # Unmethylated CHG
chh_meth = gsub(/H/, "", xm)      # Methylated CHH
chh_unmeth = gsub(/h/, "", xm)    # Unmethylated CHH

The gsub() function returns the NUMBER of substitutions it makes. So when I do gsub(/Z/, "", xm), it finds all the Z's, replaces them with nothing, and returns how many it found.

The XM encoding is:
- Z = methylated CpG
- z = unmethylated CpG  
- X = methylated CHG
- x = unmethylated CHG (note: lowercase!)
- H = methylated CHH
- h = unmethylated CHH
- . = not a cytosine

So for that example string, I get:
- cpg_meth = 1 (one Z)
- cpg_unmeth = 0 (no z's)
- chg_unmeth = 3 (three x's)
- chh_unmeth = ~15 (about 15 h's)

The lowercase vs uppercase is important - lowercase always means unmethylated. Bisulfite converts unmethylated C to T, so Bismark can distinguish methylated (C stays C) from unmethylated (C becomes T).


============================================================
STEP 4: ACCUMULATE ACROSS READS
============================================================

total_cpg_meth += cpg_meth
total_cpg_unmeth += cpg_unmeth
...
n_reads++

These are just running totals. Every read adds its counts to the overall sample-level counts. So if read 1 has 3 methylated CpGs and read 2 has 5, total_cpg_meth becomes 8. Do this across a million reads and you get genome-wide methylation statistics.

This is memory efficient because I'm not storing individual read data - just maintaining counters. Can process 50GB BAM files using basically no RAM.


============================================================
STEP 5: CALCULATE RATES
============================================================

total_cpg = total_cpg_meth + total_cpg_unmeth
cpg_rate = (total_cpg > 0) ? total_cpg_meth / total_cpg : 0

Methylation rate is straightforward: methylated cytosines divided by total cytosines in that context. The ternary operator just protects against dividing by zero if we somehow didn't see any CpGs.

Example:
- Methylated CpGs: 45,823
- Unmethylated CpGs: 17,912  
- Total CpGs: 63,735
- Methylation rate: 45,823 / 63,735 = 0.719 = 71.9%

I calculate this separately for CpG, CHG, and CHH because they have totally different biology. CpG is the main regulatory methylation in mammals (usually 60-80% methylated). CHG and CHH are mostly noise in mammals, <5% methylated.


============================================================
STEP 6: OUTPUT
============================================================

print "'$sample_id'", cpg_rate, chg_rate, chh_rate, overall_rate, total_cpg, total_chg, total_chh, n_reads

Final output is 9 columns:
SRR123456 0.719 0.012 0.008 0.114 63735 88323 281253 1000000

That's: sample ID, then the four methylation rates (CpG, CHG, CHH, overall), then the total counts for each context, then number of reads processed.


============================================================
WHY THIS MATTERS FOR THE ML PIPELINE
============================================================

These methylation rates become features for classification. But here's the interesting finding - when I ran t-tests comparing ALS vs controls:

CpG methylation: p = 0.597 (not significant)
CHG methylation: p = 0.632 (not significant)  
CHH methylation: p = 0.583 (not significant)


============================================================
DESIGN CHOICES I MADE
============================================================

1. Stream processing with awk instead of loading into Python
   - Memory efficient, handles huge files
   - Fast - processes 10M reads in ~30 seconds
   - Could parallelize across chromosomes if needed

2. Per-sample processing in Nextflow
   - Each sample is independent
   - Can run all 22 samples simultaneously
   - Easy to scale to larger cohorts

3. Output both rates and counts
   - Rates are the ML features
   - Counts let me verify coverage is adequate
   - Can check for samples with unusually low cytosine counts

4. Division by zero protection
   - Edge case but important for robustness
   - Prevents crashes on weird edge cases


============================================================
POTENTIAL ISSUES I'VE SEEN
============================================================

Main failure mode: BAM not from Bismark
- Will just output all zeros
- Quick check: samtools view file.bam | head | grep "XM:Z:"

Suspiciously low CpG methylation (<40%)
- Usually means bisulfite over-conversion or sample degradation
- Check against expected ranges for your tissue type

High CHH methylation (>10%)  
- Incomplete bisulfite conversion
- Or sequencing errors being called as methylation
- Should have conversion efficiency controls

One thing I'd add with more time: parse the XR and XG tags to separate original top/bottom strand. Right now I'm pooling both strands which is fine for most applications but strand-specific methylation can be informative.


============================================================
HOW IT FITS THE BIGGER PIPELINE
============================================================

This is one of four feature extraction modules:
1. Insert sizes (fragment lengths)
2. End motifs (fragmentation preferences)
3. Methylation (this script)
4. Positional distributions

All outputs get combined into a single feature matrix, then I do systematic feature selection to find the optimal combination.

The whole pipeline runs in about 4 minutes on AWS for all 22 samples, which is fast enough for interactive development but would need optimization for clinical-scale deployment with hundreds of samples.


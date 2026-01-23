POSITION EXTRACTION CODE WALKTHROUGH

Alright, this script is extracting where fragments are located along chromosome 21. The goal is to see if ALS fragments cluster in specific genomic regions versus being uniformly distributed.

============================================================
WHAT IT'S DOING
============================================================

Two main things:
1. Binning the genome into 1Mb windows and counting fragments in each bin
2. Dividing chr21 into three equal regions (start, middle, end) to check for broad regional enrichment

This creates two types of positional features for the ML model.


============================================================
THE CODE BREAKDOWN
============================================================

#!/bin/bash
bam_file=$1
sample_id=$2
chr21_length=46709983  # hg38 chr21 length

Standard setup. Chr21 is about 47Mb long in hg38. This is hardcoded which is fine for this project but I'd pull it from the BAM header in a production pipeline.


============================================================
EXTRACT POSITIONS
============================================================

samtools view -f 0x2 "$bam_file" | \
awk -v chr_len="$chr21_length" '{
    if ($9 > 0) {  # Only first mate

The -f 0x2 is the same as before - properly paired reads only.

The if ($9 > 0) is important. Field 9 is TLEN (template length), which is the insert size. For paired-end reads:
- First mate in pair: TLEN is positive
- Second mate in pair: TLEN is negative (same magnitude, opposite sign)

By filtering for $9 > 0, I'm only processing each fragment ONCE. Otherwise I'd count every fragment twice (once for each mate). This is a common pattern in paired-end analysis.


============================================================
CALCULATE FRAGMENT BOUNDARIES
============================================================

    start_pos = $4
    end_pos = $4 + $9

Field 4 is the alignment position (where this read starts).
Field 9 is the insert size (distance to mate).

So if a read aligns at position 10,000,000 with insert size 150:
- Fragment starts at 10,000,000
- Fragment ends at 10,000,150

This gives me the full fragment coordinates, not just where each individual read mapped.


============================================================
BIN INTO 1MB WINDOWS
============================================================

    start_bin = int(start_pos / 1000000)
    end_bin = int(end_pos / 1000000)
    bins[start_bin]++

This divides chr21 into 1Mb bins. So:
- Positions 0-999,999 → bin 0
- Positions 1,000,000-1,999,999 → bin 1
- Positions 2,000,000-2,999,999 → bin 2
- etc.

Chr21 is ~47Mb so I get about 47 bins total.

I'm only incrementing the start_bin, not the end_bin. This counts where fragments BEGIN, not where they span. For fragments that cross bin boundaries, this avoids double-counting and keeps the interpretation clean: "how many fragments start in this region?"

The bins array is an associative array in awk - it auto-creates keys as needed. Very convenient for sparse data.


============================================================
THREE-REGION DISTRIBUTION
============================================================

    total++
    if (start_pos < chr_len/3) region1++
    else if (start_pos < 2*chr_len/3) region2++
    else region3++

This is a simpler, coarser-grained view. Divide chr21 into three equal chunks:
- Region 1: positions 0 to 15.6Mb (first third)
- Region 2: positions 15.6Mb to 31.2Mb (middle third)
- Region 3: positions 31.2Mb to 46.7Mb (last third)

Count how many fragments fall in each region.

Why three regions? It's simple and biologically meaningful. Chr21 has:
- Centromeric regions (often repetitive, lower fragment recovery)
- Gene-rich regions
- Gene-poor regions

If ALS preferentially sheds cfDNA from certain chromosomal domains, this would show up here.


============================================================
OUTPUT TWO FILES
============================================================

END {
    # File 1: Binned counts (1Mb resolution)
    print "# Position bins (1Mb windows)" > "'${sample_id}_position_bins.txt'"
    for (bin in bins) {
        print bin, bins[bin] > "'${sample_id}_position_bins.txt'"
    }

This creates a file like:

Position bins (1Mb windows)
0 1234
1 1189
2 1301
3 1276
...
46 891

Each line is: bin_number fragment_count

This gives fine-grained positional information. Could use this to plot fragment density across chr21.


    # File 2: Regional stats (coarse-grained)
    print "'$sample_id'", region1, region2, region3, total
}' > ${sample_id}_position_stats.txt

This creates a single-line summary:

SRR123456 35123 38901 32456 106480

That's: sample_id, count_in_region1, count_in_region2, count_in_region3, total_fragments

These four numbers (or technically three since total = sum of regions) become ML features.


============================================================
WHAT THE FEATURES MEAN
============================================================

Region-based features test the hypothesis: "Do ALS patients have different regional fragmentation patterns?"

For example:
- Maybe ALS preferentially releases cfDNA from telomeric regions
- Maybe certain chr21 genes are dysregulated in ALS, leading to more fragments from those loci
- Maybe chromatin accessibility differs by region in disease

The ML model can learn: "ALS samples have higher region1/total ratio" or "Controls have more uniform distribution across regions."


============================================================
WHY THIS APPROACH
============================================================

1. Two resolutions for different questions
   - 1Mb bins: detailed positional signal, good for visualization
   - Three regions: simpler feature space for ML, less overfitting risk

2. Only count fragment starts
   - Avoids complexity of fragments spanning bins
   - Clear interpretation: "initiation density"
   - Prevents double-counting at boundaries

3. Use only first mate
   - Each fragment counted exactly once
   - Simple $9 > 0 filter
   - Computationally efficient



============================================================
POTENTIAL IMPROVEMENTS
============================================================

1. Normalize by mappability
   - Some regions are easier to align to than others
   - Low-complexity regions get fewer reads
   - Could bias fragment counts

2. Account for CpG density
   - CpG islands affect fragmentation
   - Could normalize bins by CpG content

3. Gene-based binning instead of fixed 1Mb
   - Maybe fragments cluster around genes specifically
   - Could use gene annotations from GTF file

4. Strand-specific analysis
   - Forward vs reverse strand fragment origins
   - Might reveal transcription-coupled fragmentation


============================================================
BENEFITS
============================================================

1. Validates the approach
   - Systematically tested positions, found they don't matter
   - Better than just assuming they wouldn't help

2. Biological insight
   - Knowing positions DON'T matter is valuable
   - Tells us ALS signature is intrinsic to fragment properties

3. Quality control
   - Can check for weird coverage artifacts
   - Ensures fragments distribute reasonably across chr21
   - Would catch sample prep issues (like enrichment bias)

4. Future-proofing
   - If someone wants to try full-genome later
   - Already have the infrastructure to extract positions
   - Just need to modify chr_len and loop over chromosomes


============================================================
DESIGN TRADEOFFS
============================================================

Fixed 1Mb bins vs adaptive bins:
- Fixed: Simple, predictable, easy to compare across samples
- Adaptive (e.g., by gene density): More biologically meaningful but harder to interpret

Three equal regions vs quartiles vs deciles:
- Three: Minimal features, less overfitting with small N
- More: Higher resolution, but risk of overfitting with only 22 samples



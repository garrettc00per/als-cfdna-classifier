#!/bin/bash
# Extract fragment start and end positions

bam_file=$1
sample_id=$2

# Get chromosome 21 length (we'll bin positions)
chr21_length=46709983  # hg38 chr21 length

# Extract fragment positions
samtools view -f 0x2 "$bam_file" | \
awk -v chr_len="$chr21_length" '{
    if ($9 > 0) {  # Only first mate
        start_pos = $4
        end_pos = $4 + $9
        
        # Bin into 1Mb windows for density calculation
        start_bin = int(start_pos / 1000000)
        end_bin = int(end_pos / 1000000)
        
        # Count fragments in each bin
        bins[start_bin]++
        
        # Track overall distribution
        total++
        if (start_pos < chr_len/3) region1++
        else if (start_pos < 2*chr_len/3) region2++
        else region3++
    }
}
END {
    # Output binned counts
    print "# Position bins (1Mb windows)" > "'${sample_id}_position_bins.txt'"
    for (bin in bins) {
        print bin, bins[bin] > "'${sample_id}_position_bins.txt'"
    }
    
    # Output regional distribution
    print "'$sample_id'", region1, region2, region3, total
}' > ${sample_id}_position_stats.txt

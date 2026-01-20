#!/bin/bash
# Extract 4bp end motifs from fragment reads

bam_file=$1
sample_id=$2

# Extract sequences from properly paired reads
# Get 5' end (start) motifs from read1 and 3' end motifs from read2
samtools view -f 0x2 "$bam_file" | \
awk '{
    if ($9 > 0) {  # First in pair
        seq = $10
        if (length(seq) >= 4) {
            # Start motif: first 4bp of read1
            print "START", substr(seq, 1, 4)
        }
    }
    if ($9 < 0) {  # Second in pair  
        seq = $10
        if (length(seq) >= 4) {
            # End motif: first 4bp of read2 (reverse complement represents fragment end)
            print "END", substr(seq, 1, 4)
        }
    }
}' > ${sample_id}_motifs_raw.txt

# Count start motifs
echo "=== START MOTIFS ===" > ${sample_id}_end_motifs.txt
grep "^START" ${sample_id}_motifs_raw.txt | \
    awk '{print $2}' | \
    sort | uniq -c | sort -rn >> ${sample_id}_end_motifs.txt

# Count end motifs  
echo "=== END MOTIFS ===" >> ${sample_id}_end_motifs.txt
grep "^END" ${sample_id}_motifs_raw.txt | \
    awk '{print $2}' | \
    sort | uniq -c | sort -rn >> ${sample_id}_end_motifs.txt

# Clean up
rm ${sample_id}_motifs_raw.txt

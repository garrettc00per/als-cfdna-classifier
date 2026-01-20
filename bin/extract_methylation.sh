#!/bin/bash
# Extract methylation features from XM tags

bam_file=$1
sample_id=$2

# Extract XM tags and calculate methylation stats
samtools view -f 0x2 "$bam_file" | \
awk '{
    for(i=12; i<=NF; i++) {
        if($i ~ /^XM:Z:/) {
            # Extract the XM tag value
            xm = substr($i, 6)
            
            # Count different methylation types
            cpg_meth = gsub(/Z/, "", xm)      # Methylated CpG
            cpg_unmeth = gsub(/z/, "", xm)    # Unmethylated CpG
            chg_meth = gsub(/X/, "", xm)      # Methylated CHG
            chg_unmeth = gsub(/x/, "", xm)    # Unmethylated CHG
            chh_meth = gsub(/H/, "", xm)      # Methylated CHH
            chh_unmeth = gsub(/h/, "", xm)    # Unmethylated CHH
            
            # Accumulate counts
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
    # Calculate methylation rates
    total_cpg = total_cpg_meth + total_cpg_unmeth
    total_chg = total_chg_meth + total_chg_unmeth
    total_chh = total_chh_meth + total_chh_unmeth
    total_c = total_cpg + total_chg + total_chh
    
    cpg_rate = (total_cpg > 0) ? total_cpg_meth / total_cpg : 0
    chg_rate = (total_chg > 0) ? total_chg_meth / total_chg : 0
    chh_rate = (total_chh > 0) ? total_chh_meth / total_chh : 0
    overall_rate = (total_c > 0) ? (total_cpg_meth + total_chg_meth + total_chh_meth) / total_c : 0
    
    # Print results
    print "'$sample_id'", cpg_rate, chg_rate, chh_rate, overall_rate, total_cpg, total_chg, total_chh, n_reads
}' > ${sample_id}_methylation.txt

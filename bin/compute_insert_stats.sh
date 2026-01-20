#!/bin/bash
# Usage: compute_insert_stats.sh insert_sizes.txt sample_id

input_file=$1
sample_id=$2

awk -v sample="$sample_id" '{
    sizes[NR] = $1
    sum += $1
    sumsq += $1*$1
    if ($1 > max) max = $1
    if (min == "" || $1 < min) min = $1
}
END {
    n = NR
    mean = sum/n
    # Calculate median (sort array)
    asort(sizes)
    if (n % 2) {
        median = sizes[(n+1)/2]
    } else {
        median = (sizes[n/2] + sizes[n/2+1])/2
    }
    q25 = sizes[int(n*0.25)]
    q75 = sizes[int(n*0.75)]
    stddev = sqrt(sumsq/n - mean*mean)
    
    print sample, mean, median, stddev, min, max, q25, q75, n
}' "$input_file"

#!/bin/bash
# Add disease status to the summary file

stats_file=$1
metadata=$2

# Print header with disease_status column
echo "sample_id mean median stddev min max q25 q75 n_fragments disease_status"

# Add disease status for each sample
tail -n +2 "$stats_file" | while read line; do
    sample_id=$(echo "$line" | awk '{print $1}')
    # Column 21 when naively split by comma (because of quoted fields)
    disease=$(awk -F',' -v id="$sample_id" '$1 == id {print $21}' "$metadata")
    echo "$line $disease"
done

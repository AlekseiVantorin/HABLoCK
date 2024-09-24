#!/bin/bash

# Input file
triobin_file="./contigs_HG06808.txt"

input_dump="../../mnt/projects/assembly/project2/crossingover/counted_dump/HG06808/"
output_dump="../../mnt/projects/assembly/project2/crossingover/counted_dump/HG06808_sorted/"

# Read each line of the input file
while IFS= read -r line; do
    # Extract the first word
    first_word=$(echo "$line" | awk '{print $1}')

    sort -o "${output_dump}HG06808_${first_word}_sorted.txt" "${input_dump}HG06808_${first_word}.txt"

done < "$triobin_file"

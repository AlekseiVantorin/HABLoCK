#!/bin/bash

# Input file
triobin_file="./contigs_HG06808.txt"

output_kmc="../../mnt/projects/assembly/project2/crossingover/counted_kmc/HG06808/"
output_dump="../../mnt/projects/assembly/project2/crossingover/counted_dump/HG06808/"

# Read each line of the input file
while IFS= read -r line; do
    # Extract the first word
    first_word=$(echo "$line" | awk '{print $1}')

    child_input_file="../../mnt/projects/assembly/project2/crossingover/counted_contigs/HG06808/HG06808_${first_word}.fasta"
    ./kmc -k27 -m32 -ci1 -fa "${child_input_file}" "${output_kmc}HG06808_${first_word}" .
    ./kmc_dump "${output_kmc}HG06808_${first_word}" "${output_dump}HG06808_${first_word}.txt"

done < "$triobin_file"

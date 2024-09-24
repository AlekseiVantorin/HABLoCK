#!/bin/bash

# Input file
triobin_file="./contigs_HG06808.txt"

child_dump="../../mnt/projects/assembly/project2/crossingover/counted_dump/HG06808_sorted/HG06808_"
maternal_dump="../../mnt/projects/assembly/project2/crossingover/counted_dump/sorted_dump_pan010.txt"
paternal_dump="../../mnt/projects/assembly/project2/crossingover/counted_dump/sorted_dump_pan011.txt"

output_path="../../mnt/projects/assembly/project2/crossingover/counted_intersections/HG06808/"

# Read each line of the input file
while IFS= read -r line; do
    # Extract the first word
    first_word=$(echo "$line" | awk '{print $1}')

#   ./intersection paternal_dump maternal_dump child_dump paternal_origin maternal_origin
    ./intersection "${paternal_dump}" "${maternal_dump}" "${child_dump}${first_word}_sorted.txt" "${output_path}HG06808_pan011_${first_word}.txt" "${output_path}HG06808_pan010_${first_word}.txt"

done < "$triobin_file"

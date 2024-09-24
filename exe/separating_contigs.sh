#!/bin/bash
cd ../../mnt/

# Input file
triobin_file="../home/avantorin/contigs_HG06808.txt"

# Paternal FASTA file
#parental_path="./projects/assembly/project2/samples/HG06807/parents/"
#paternal_blood="pan011/blood/"
#paternal_cell_one="pan011/cell1/"
#paternal_cell_two="pan011/cell2/"

#maternal_blood="pan010/blood/"
#maternal_cell_one="pan010/cell1/"
#maternal_cell_two="pan010/cell2/"

#paternal_blood_one="CCATTCTGTA-TAACTTCACG_S36_L003_R1_001.fastq.gz"
#paternal_cell_one_one="HTMKJDSX2_CGGCGTTCTA-GAGTTGAGGC_L004_R1.fastq.gz"
#paternal_cell_two_one="HTTMGDSX2_CGGCGTTCTA-GAGTTGAGGC_L002_R1.fastq.gz"
#paternal_blood_two="CCATTCTGTA-TAACTTCACG_S36_L003_R2_001.fastq.gz"
#paternal_cell_one_two="HTMKJDSX2_CGGCGTTCTA-GAGTTGAGGC_L004_R2.fastq.gz"
#paternal_cell_two_two="HTTMGDSX2_CGGCGTTCTA-GAGTTGAGGC_L002_R2.fastq.gz"

#maternal_blood_one="AATGTCCGAG-CTTACGAGTA_S10_L003_R1_001.fastq.gz"
#maternal_cell_one_one="HTMKJDSX2_TCGCTCTGGA-TACTCAGCTG_L004_R1.fastq.gz"
#maternal_cell_two_one="HTTMGDSX2_TCGCTCTGGA-TACTCAGCTG_L002_R1.fastq.gz"
#maternal_blood_two="AATGTCCGAG-CTTACGAGTA_S10_L003_R2_001.fastq.gz"
#maternal_cell_one_two="HTMKJDSX2_TCGCTCTGGA-TACTCAGCTG_L004_R2.fastq.gz"
#maternal_cell_two_two="HTTMGDSX2_TCGCTCTGGA-TACTCAGCTG_L002_R2.fastq.gz"

child_fasta="./projects/dantipov/assemblies/HG06808/assembly.fasta"

output_path="./projects/assembly/project2/crossingover/counted_contigs/"

# Read each line of the input file
while IFS= read -r line; do
    # Extract the first word
    first_word=$(echo "$line" | awk '{print $1}')

    # Create temporary file tmp.reg and put the first word in it
    echo "$first_word" > "../home/avantorin/tmp.reg"

    # Run seqtk subseq with the temporary file and paternal FASTA
#    paternal_output_file_blood_one="${output_path}${paternal_blood}pan011_${first_word}_R1.fastq"
#    paternal_output_file_blood_two="${output_path}${paternal_blood}pan011_${first_word}_R2.fastq"
#    paternal_output_file_cell_one_one="${output_path}${paternal_cell_one}pan011_${first_word}_R1.fastq"
#    paternal_output_file_cell_one_two="${output_path}${paternal_cell_one}pan011_${first_word}_R2.fastq"
#    paternal_output_file_cell_two_one="${output_path}${paternal_cell_two}pan011_${first_word}_R1.fastq"
#    paternal_output_file_cell_two_two="${output_path}${paternal_cell_two}pan011_${first_word}_R2.fastq"

#    ../home/avantorin/seqtk subseq "${parental_path}${paternal_blood}${paternal_blood_one}" tmp.reg > "${paternal_output_file_blood_one}"
#    ../home/avantorin/seqtk subseq "${parental_path}${paternal_blood}${paternal_blood_two}" tmp.reg > "${paternal_output_file_blood_two}"
#    ../home/avantorin/seqtk subseq "${parental_path}${paternal_blood}${paternal_cell_one_one}" tmp.reg > "${paternal_output_file_cell_one_one}"
#    ../home/avantorin/seqtk subseq "${parental_path}${paternal_blood}${paternal_cell_one_two}" tmp.reg > "${paternal_output_file_cell_one_two}"
#    ../home/avantorin/seqtk subseq "${parental_path}${paternal_blood}${paternal_cell_two_one}" tmp.reg > "${paternal_output_file_cell_two_one}"
#    ../home/avantorin/seqtk subseq "${parental_path}${paternal_blood}${paternal_cell_two_two}" tmp.reg > "${paternal_output_file_cell_two_two}"

#    maternal_output_file_blood_one="${output_path}${maternal_blood}pan010_${first_word}_R1.fastq"
#    maternal_output_file_blood_two="${output_path}${maternal_blood}pan010_${first_word}_R2.fastq"
#    maternal_output_file_cell_one_one="${output_path}${maternal_cell_one}pan010_${first_word}_R1.fastq"
#    maternal_output_file_cell_one_two="${output_path}${maternal_cell_one}pan010_${first_word}_R2.fastq"
#    maternal_output_file_cell_two_one="${output_path}${maternal_cell_two}pan010_${first_word}_R1.fastq"
#    maternal_output_file_cell_two_two="${output_path}${maternal_cell_two}pan010_${first_word}_R2.fastq"

#    ../home/avantorin/seqtk subseq "${parental_path}${maternal_blood}${maternal_blood_one}" tmp.reg > "${maternal_output_file_blood_one}"
#    ../home/avantorin/seqtk subseq "${parental_path}${maternal_blood}${maternal_blood_two}" tmp.reg > "${maternal_output_file_blood_two}"
#    ../home/avantorin/seqtk subseq "${parental_path}${maternal_blood}${maternal_cell_one_one}" tmp.reg > "${maternal_output_file_cell_one_one}"
#    ../home/avantorin/seqtk subseq "${parental_path}${maternal_blood}${maternal_cell_one_two}" tmp.reg > "${maternal_output_file_cell_one_two}"
#    ../home/avantorin/seqtk subseq "${parental_path}${maternal_blood}${maternal_cell_two_one}" tmp.reg > "${maternal_output_file_cell_two_one}"
#    ../home/avantorin/seqtk subseq "${parental_path}${maternal_blood}${maternal_cell_two_two}" tmp.reg > "${maternal_output_file_cell_two_two}"

    child_output_file="./projects/assembly/project2/crossingover/counted_contigs/HG06808/HG06808_${first_word}.fasta"
    ../home/avantorin/seqtk subseq "$child_fasta" "../home/avantorin/tmp.reg" > "$child_output_file"

    # Clear tmp.reg for the next iteration
    > "../home/avantorin/tmp.reg"
done < "$triobin_file"

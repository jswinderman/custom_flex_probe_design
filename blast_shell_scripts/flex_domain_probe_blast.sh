#!/bin/bash

# Define variables
database="hs.refseq.rna"
input_file="231112_fd_sub_lhs_rhs_probes.fasta"
output_file="231112_FD_sub_probe_blast.txt"
min_length=20  # Minimum alignment length

# Run BLAST and redirect output to a temporary file
blastn -db "$database" -query "$input_file" -task "blastn-short" -out "$output_file.tmp" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"

# Filter the results using AWK and redirect to the final output file
awk -v min="$min_length" '$4 >= min {print}' "$output_file.tmp" > "$output_file"

# Remove the temporary file
rm "$output_file.tmp"
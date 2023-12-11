#!/bin/bash


# download and create BLAST databases for human refseq genomic
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.*.rna.gbff.gz
for a in human.*.rna.gbff.gz; do gunzip $a; done


wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.*.rna.fna.gz
for a in human.*.rna.fna.gz; do gunzip $a; done

# Concatenate FNA files
cat human.*.rna.fna > concatenated_human_rna.fna

makeblastdb -dbtype nucl \
    -parse_seqids \
    -in concatenated_human_rna.fna \
    -out hs.refseq.rna \
    -title "hsapiens_refseq_"

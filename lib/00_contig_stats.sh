#!/bin/bash

###=====================###
### CONFIGURATION START ###
###=====================###

# FASTA file of preliminary assembled genome
genomeFastaFile="/path/to/prelimiary_genome_fasta_file.fa"

# Output directory
outputDirectory="/path/to/output_directory/"

###=====================###
### CONFIGURATION END   ###
###=====================###

contigs=`cat $genomeFastaFile | grep "^>" | awk -F " " '{print $1}' | sed 's/>//' | sort -n`
samtools faidx $genomeFastaFile $contigs | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | sort -k2n | awk '{sum+=$2} END {print sum}'
cat $genomeFastaFile | grep "^>" | awk -F " " '{print $1}' | sed 's/>//' | sort -n` | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | sort -k2n > $outputDirectory/contigLengths.txt

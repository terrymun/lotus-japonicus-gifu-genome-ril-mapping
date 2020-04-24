#!/bin/bash

###=====================###
### CONFIGURATION START ###
###=====================###

# Output directory
outputDirectory="/path/to/output_directory/"

# Individual folders for the output of each RIL population
gifuBurtiiOutputDirectory="${outputDirectory}/Gifu-Burtii"
gifuMG20OutputDirectory="${outputDirectory}/Gifu-MG20"

# FASTA file of preliminary assembled genome
genomeFastaFile="/path/to/prelimiary_genome_fasta_file.fa"

###=====================###
### CONFIGURATION END   ###
###=====================###


# Make folders for each RIL population
`mkdir -p $gifuBurtiiOutputDirectory`
`mkdir -p $gifuMG20OutputDirectory`

# Create bam file list containing the biparental bam files for each population
echo "${outputDirectory}/Gifu/Gifu.picardDeduped.bam" > "${gifuBurtiiOutputDirectory}/bamFileList.txt"
echo "${outputDirectory}/Burttii/Burttii.merged.bam" >> "${gifuBurtiiOutputDirectory}/bamFileList.txt"

echo "${outputDirectory}/Gifu/Gifu.picardDeduped.bam" > "${gifuMG20OutputDirectory}/bamFileList.txt"
echo "${outputDirectory}/MG20/MG20.merged.bam" >> "${gifuMG20OutputDirectory}/bamFileList.txt"

# Generate contig regions to perform mpileup in parallel
regions=`cat $genomeFastaFile | grep "^>" | sed 's/>//' |  cut -d ' ' -f1 | sort -n | uniq`

# Gifu x Burttii RILs start with C<...>
`find $outputDirectory -type f -name "C*.merged.bam" | sort >> $gifuBurtiiOutputDirectory/bamFileList.txt`

# Gifu x MG20 RILs start with R<...>
`find $outputDirectory -type f -name "R*.merged.bam" | sort >> $gifuMG20OutputDirectory/bamFileList.txt`

# Perform mpileup for each RIL population
for region in $regions; do
  for rilDir in $(find ${outputDirectory} -type d -name 'Gifu-*'); do
    ###======================###
    ### JOB SUBMISSION START ###
    ###======================###
    jobCommand="samtools mpileup \
      -f ${genomeFastaFile} \
      -v -I \
      -r $region \
      -b ${rilDir}/bamFileList.txt | 
    tee ${rilDir}/variants_${region}.vcf |
    bcftools call \
      -m -v \
      -Oz \
      --threads==1 \
      - > ${rilDir}/variantsCalled_${region}.vcf && \
    echo \"Contig: ${region}\""
    ###======================###
    ### JOB SUBMISSION END   ###
    ###======================###
  done
done

#!/bin/bash

###=====================###
### CONFIGURATION START ###
###=====================###

# Output directory
outputDirectory="/path/to/output_directory/"

# FASTA file of preliminary assembled genome
genomeFastaFile="/path/to/prelimiary_genome_fasta_file.fa"

###=====================###
### CONFIGURATION END   ###
###=====================###



# Make folders for each RIL population
`mkdir -p ${outputDirectory}/Gifu-Burtii`
`mkdir -p ${outputDirectory}/Gifu-MG20`

# Create bam file list containing the biparental bam files for each population
echo "${outputDirectory}/Gifu/Gifu.picardDeduped.bam" > "${outputDirectory}/Gifu-Burtii/bamFileList.txt"
echo "${outputDirectory}/Burtii/Burtii.merged.bam" >> "${outputDirectory}/Gifu-Burtii/bamFileList.txt"

echo "${outputDirectory}/Gifu/Gifu.picardDeduped.bam" > "${outputDirectory}/Gifu-MG20/bamFileList.txt"
echo "${outputDirectory}/MG20/MG20.merged.bam" >> "${outputDirectory}/Gifu-MG20/bamFileList.txt"

# Generate contig regions to perform mpileup in parallel
regions=`cat $genomeFastaFile | grep "^>" | sed 's/>//' |  cut -d ' ' -f1 | sort -n | uniq`

# Gifu x Burtii RILs start with C<...>
`find $outputDirectory -type f -name "C*.merged.bam" | sort >> ${outputDirectory}/Gifu-Burtii/bamFileList.txt`

# Gifu x MG20 RILs start with R<...>
`find $outputDirectory -type f -name "R*.merged.bam" | sort >> ${outputDirectory}/Gifu-MG20/bamFileList.txt`

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

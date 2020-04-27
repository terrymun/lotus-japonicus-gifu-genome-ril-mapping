#!/bin/bash

###=====================###
### CONFIGURATION START ###
###=====================###

# Output directory
outputDirectory="/path/to/output_directory/"

###=====================###
### CONFIGURATION END   ###
###=====================###



for rilDir in $(find ${outputDirectory} -type d -name 'Gifu-*'); do
  # Create list of VCF files in directory
  # Then run concatenation job
  find $rilDir -type f -name "refiltered_*.vcf" | sort -n > ${rilDir}/refiltered_vcfList.txt && \
  ###======================###
  ### JOB SUBMISSION START ###
  ###======================###
  jobCommand="bcftools concat \
    -f ${rilDir}/refiltered_vcfList.txt \
    -Ov > ${rilDir}/refiltered.merged.vcf && \
  vcftools \
    --vcf ${rilDir}/refiltered.merged.vcf \
    --out ${rilDir}/refiltered.merged \
    --012"
  ###======================###
  ### JOB SUBMISSION END ###
  ###======================###
done

for rilDir in $(find ${outputDirectory} -type d -name 'Gifu-*'); do
  # Create list of VCF files in directory
  # Then run concatenation job
  find $rilDir -type f -name "variants_*.vcf" | sort -n > ${rilDir}/vcfList.txt && \
  ###======================###
  ### JOB SUBMISSION START ###
  ###======================###
  jobCommand="bcftools concat \
    -f ${rilDir}/vcfList.txt \
    -Ov > ${rilDir}/merged.vcf && \
  vcftools \
    --vcf ${rilDir}/merged.vcf \
    --out ${rilDir}/merged \
    --012"
  ###======================###
  ### JOB SUBMISSION END   ###
  ###======================###
done

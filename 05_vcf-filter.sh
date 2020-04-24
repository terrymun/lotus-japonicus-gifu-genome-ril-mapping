#!/bin/bash

###=====================###
### CONFIGURATION START ###
###=====================###

# Output directory
outputDirectory="/path/to/output_directory/"

# Minimum quality
minimumQuality="30"

# Minimum depth
minimumDepth="50"

###=====================###
### CONFIGURATION END   ###
###=====================###



for vcf in $(find ${outputDirectory} -type f -name "variantsCalled*.vcf" | sort -n); do
  ###======================###
  ### JOB SUBMISSION START ###
  ###======================###
  jobCommand="bcftools filter \
      --threads=1 \
      -Ou \
      --include 'QUAL>=${minimumQuality} & DP>=${minimumDepth}' ${vcf} | \
    bcftools view \
      --threads=1 \
      -m2 -M2 -v snps \
      -Ov - > ${vcf/variantsCalled/filtered_Q30DP50}"
  ###======================###
  ### JOB SUBMISSION END   ###
  ###======================###
done

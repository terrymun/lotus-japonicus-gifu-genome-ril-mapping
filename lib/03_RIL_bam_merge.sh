#!/bin/bash

###=====================###
### CONFIGURATION START ###
###=====================###

# Output directory
outputDirectory="/path/to/output_directory/"

###=====================###
### CONFIGURATION END   ###
###=====================###



# Go through all RIL folders
for ril in $(find ${outputDirectory} -maxdepth 1 -mindepth 1 -type d | sort); do
  
  # Get basename
  baseRIL=`basename $ril`
  
  # 1. Get all the .sorted.bam files and merge them into one
  # 2. Index the merged bam file
  # 3. Run mpileup
  # 4. Pass to bcftools for variant calling
  # 5. Filter low quality variants
  bamFiles=`find ${ril} -type f -name "*.picardDeduped.bam" | sort | paste -d ' ' -s`

  ###======================###
  ### JOB SUBMISSION START ###
  ###======================###
  jobCommand="sambamba merge \
    -t 4 \
    ${ril}/${baseRIL}.merged.bam ${bamFiles}"
  ###======================###
  ### JOB SUBMISSION END   ###
  ###======================###
done

#!/bin/bash

###=====================###
### CONFIGURATION START ###
###=====================###

# Output directory
outputDirectory="/path/to/output_directory/"

###=====================###
### CONFIGURATION END   ###
###=====================###



# Find all sorted BAM files
for file in $(find ${outputDirectory} -type f -name "*.sorted.bam" | sort); do
  ###======================###
  ### JOB SUBMISSION START ###
  ###======================###
  jobCommand="picard MarkDuplicates \
  I=${file} \
  O=${file/.sorted.bam/.picardDeduped.bam/} \
  METRICS_FILE=${file/.sorted.bam/_dupMetrics.txt/} \
  VALIDATION_STRINGENCY=SILENT \
  REMOVE_DUPLICATES=TRUE"

  ###======================###
  ### JOB SUBMISSION END   ###
  ###======================###
done

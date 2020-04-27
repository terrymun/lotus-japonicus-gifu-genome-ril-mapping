#!/bin/bash

###=====================###
### CONFIGURATION START ###
###=====================###

# Output directory
outputDirectory="/path/to/output_directory/"

###=====================###
### CONFIGURATION END   ###
###=====================###



# Remove previous refiltering files
`find $outputDirectory -type f -name "refiltered_*.vcf" -exec rm {} \;`

# Run refiltering job
for rilPop in $(find ${outputDirectory} -type d -name "Gifu-*"); do
  ###======================###
  ### JOB SUBMISSION START ###
  ###======================###
  jobCommand="python ./homref_het_filter.py -l ${rilPop}/refilteringLog.csv ${rilPop}/filtered_*.vcf"
  ###======================###
  ### JOB SUBMISSION END   ###
  ###======================###
done

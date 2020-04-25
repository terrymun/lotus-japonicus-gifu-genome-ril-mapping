#!/bin/bash

###=====================###
### CONFIGURATION START ###
###=====================###

# Output directory
outputDirectory="/path/to/output_directory/"

###=====================###
### CONFIGURATION END   ###
###=====================###



# Gifu x Burtii population
###======================###
### JOB SUBMISSION START ###
###======================###
jobCommand="python process-genotype-calls.py -i ${outputDirectory}/Gifu-Burtii/refiltered.merged.012 --popType RIL8 --popName GifuBurtii --cherryPick 0 --window 5 --noMapSize 0"
###======================###
### JOB SUBMISSION END   ###
###======================###

# Gifu x MG20 population
###======================###
### JOB SUBMISSION START ###
###======================###
jobCommand="python process-genotype-calls.py -i ${outputDirectory}/Gifu-MG20/refiltered.merged.012 --popType RIL9 --popName GifuMG20 --cherryPick 0 --window 5 --noMapSize 0"
###======================###
### JOB SUBMISSION END   ###
###======================###

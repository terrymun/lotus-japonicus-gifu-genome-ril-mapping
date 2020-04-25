#!/bin/bash

###=====================###
### CONFIGURATION START ###
###=====================###

# NOTE: You will need to define the path to the directory that contains:
#       - paired-end reads for Gifu×Burtii RIL population
#       - paired-end reads for MG20×Gifu RIL population
gifuBurtiiDir="/path/to/Gifu_x_Burtii_RIL_reads/"
mg20GifuDir="/path/to/MG20_x_Gifu_RIL_reads/"

# Temporary working directory for BWA-MEM
tmpDir="/path/to/temp_working_directory/"

# FASTA file of preliminary assembled genome
genomeFastaFile="/path/to/prelimiary_genome_fasta_file.fa"

# FASTQ file for paired end reads of Gifu
gifuReadsPair1="/path/to/gifu_paired_end_reads_pair1.fq"
gifuReadsPair2="/path/to/gifu_paired_end_reads_pair1.fq"

# Output directory
outputDirectory="/path/to/output_directory/"

###=====================###
### CONFIGURATION END   ###
###=====================###



# Create output directory
`mkdir -p $outputDirectory`

# Get fastq files from both RILs
list1=`find $gifuBurtiiDir -type f -name "*.fastq" | grep "_R1_001.fastq" | sed 's/_R1_001.fastq//' | sort`
list2=`find $mg20GifuDir -type f -name "*fastq" | grep "_R1_001.fastq" | sed 's/_R1_001.fastq//' | sort`

for file in $list1 $list2; do

  # Get the RIL name
  RILname=$(echo `basename $(dirname ${file}_R1_001.fastq)` | sed 's/Sample_//' | sed 's/\-//')
  echo "RIL name is: ${RILname}"
  
  # Get basefile name
  Basefile=`basename ${file}`
  echo "Basefile name is: ${Basefile}"
  
  # Make directory for each RIL
  `mkdir -p ${outputDirectory}/${RILname}`
  
  # Store variables for readgroup entries, required by GATK
  RGPL="ILLUMINA"
  RGSM=`echo $Basefile | cut -d '_' -f 1`
  RGLN=`echo $Basefile | cut -d '_' -f 3`
  RGID="${RGSM}_${RGLN}"
  RGLB="foo"
  readgroup="@RG\tID:${RGID}\tPL:${RGPL}\tSM:${RGSM}\tLB:${RGLB}"
  
  # Job submission
  # 1. Use bwa-mem for alignment
  # 2. Use samtools view to convert to bam
  # 3. Use samtools sort to sort bam file
  ###======================###
  ### JOB SUMBISSION START ###
  ###======================###
  jobCommand="bwa mem \
  -t 8 \
  -R \"${readgroup}\" \
  ${genomeFastaFile} \
  ${file}_R1_001.fastq ${file}_R2_001.fastq | \
  sambamba view -f bam -S -t 4 /dev/stdin | \
  sambamba sort -m 16GB --tmpdir ${tmpDir} -t 4 -o ${outputDirectory}/${RILname}/${Basefile}.sorted.bam /dev/stdin"
  ###======================###
  ### JOB SUMBISSION END   ###
  ###======================###
done

# Genearte BAM file for Gifu
###======================###
### JOB SUMBISSION START ###
###======================###
jobCommand="bwa mem \
-t 16 \
-R \"@RG\tID:Gifu\tPL:ILLUMINA\tSM:Gifu\tLB:foo\" \
${genomeFastaFile} ${gifuReadsPair1} ${gifuReadsPair2} | \
sambamba view -f bam -S -t 16 /dev/stdin | \
sambamba sort -m 128GB --tmpdir ${tmpDir} -t 16 -o ${outputDirectory}/Gifu/Gifu.sorted.bam /dev/stdin"
###======================###
### JOB SUMBISSION END   ###
###======================###

# Generate BAM file for MG20
# TBD

# Geneate BAM file for Burtii
# TBD

# RIL mapping

This readme contains all the necessary scripts needed

# Pre-requisites

The following datasets must be available on your machine. They can be downloaded from NCBI's SRA website in the links provided below:

* Gifu paired end reads: https://submit.ncbi.nlm.nih.gov/subs/sra/SUB4699487/overview
* All the paired-end read `.fastq` files for *Lotus japonicus* RIL populations:
  * Gifu&times;Burtii: https://submit.ncbi.nlm.nih.gov/subs/sra/SUB4699503/overview, and
  * MG20&times;Gifu: https://submit.ncbi.nlm.nih.gov/subs/sra/SUB4699504/overview

The following tools will need to be installed and be available in your environment before running any of the bash scripts in this repository. The version number indicated for each tool used is the minimum version used when analysis was done, and verified to be working.

* [bcftools 1.3](https://github.com/samtools/bcftools)
* [bwa 0.7.15](https://github.com/lh3/bwa)
* [picard 2.7.1](https://github.com/broadinstitute/picard) (run in Java 1.8.0 environment)
* [sambamba 0.6.5](https://github.com/biod/sambamba)
* [samtools 1.3](https://github.com/samtools/samtools)
* [vcftools 0.1.15](https://github.com/vcftools/vcftools)

In addition, the final step, `08_variant-smoothing.sh` runs a python script that requires [numpy](https://numpy.org/) to be installed. At the time of writing, the Python version used is 2.7.13.

# How to use

The bash scripts should be run sequentially, starting from the lowest index, `01`, and up till `08`. The scripts have been intentionally split into several files to allow catching of errors or fine-tuning of parameters in between step, instead of combining them into a large file. It is, of course, possible to automate the entire process, but this is not part of the scope of this repository.

## Manual configuration required

Each bash script contains certain variables that must be configured manually on your end. They are marked within the blocks at the top of the files:

```bash
###=====================###
### CONFIGURATION START ###
###=====================###

# (configuration here)

###=====================###
### CONFIGURATION END   ###
###=====================###
```

## Job submission to your own computing cluster

The scripts contains markers that indicate certain codes that should be executed as part of a job submission script in the computing cluster that you are using. These operations are computationally expensive and should not run on the same thread as your terminal. These lines are marked with:

```bash
###======================###
### JOB SUBMISSION START ###
###======================###

jobCommand="(arbitrary command here)"

###======================###
### JOB SUBMISSION END   ###
###======================###
```

The `jobCommand` should be executed as in by your computing cluster's own job management and dispatch system. Please consult with your syster administrator or person-in-charge if you have any doubts on submitting jobs to a queue.

## What about the Python scripts?

TL;DR: Leave them alone.

The python files are simply dependencies used by the bash scripts and **should not be executed out of context or order, and should never be used individually**.
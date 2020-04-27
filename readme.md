# RIL mapping for *Lotus japonicus* Gifu genome assembly

This readme contains all the necessary scripts needed to repeat the mapping step in the assembly efforts behind *Lotus japonicus* Gifu v1.2 genome. It deals with how we have used paired-end reads obtained from two distinct recombinant inbred line (RIL) populations, namely Gifu&times;Burtii and MG20&times;Gifu, to polish our assembly.

# Pre-requisites

The following datasets must be available on your machine. They can be downloaded from NCBI's SRA website in the links provided below:

* Paired-end reads from the parent lines:
  * Gifu: https://www.ebi.ac.uk/ena/data/view/SAMEA4807266 (part of the [PRJEB27969](https://www.ebi.ac.uk/ena/data/view/PRJEB27969) study),
  * Burtii: [TBA], and
  * MG20: https://www.ebi.ac.uk/ena/data/view/SAMEA4807284 (part of the [PRJEB27969](https://www.ebi.ac.uk/ena/data/view/PRJEB27969) study)

* All the paired-end reads for the RIL populations:
  * Gifu&times;Burtii: https://submit.ncbi.nlm.nih.gov/subs/sra/SUB4699503/overview, and
  * MG20&times;Gifu: https://submit.ncbi.nlm.nih.gov/subs/sra/SUB4699504/overview

The following tools will need to be installed and be available in your environment before running any of the bash scripts in this repository. The version number indicated for each tool used is the minimum version used when analysis was done, and verified to be working.

* [bcftools 1.3](https://github.com/samtools/bcftools)
* [bwa 0.7.15](https://github.com/lh3/bwa)
* [picard 2.7.1](https://github.com/broadinstitute/picard) (run in Java 1.8.0 environment)
* [sambamba 0.6.5](https://github.com/biod/sambamba)
* [samtools 1.3](https://github.com/samtools/samtools)
* [vcftools 0.1.15](https://github.com/vcftools/vcftools)

The two final steps of the RIL mapping pipelines requires running python scripts from bash, namely `07_vcf-concat.sh` and `08_variant-smoothing.sh`. Their dependencies are:

* [Python 2.7.13](https://www.python.org/downloads/release/python-2713/)
* [NumPy](https://numpy.org/)

If you want to perform data visualization as the final step, you will need to have R installed and the following packages installed, too:

* ggplot2
* gplots
* reshape2
* optparse
* ape
* RColorBrewer

R packages can be installed by simply running `install.packages(PACKAGE_NAME)`.

# How to use

## Pipeline (see the `/lib` directory)

The pipeline consists of bash scripts, located in the `/lib` directory, that should be run sequentially, starting from the lowest index, `01`, and up till `08`. The scripts have been intentionally split into several files to allow catching of errors or fine-tuning of parameters in between step, instead of combining them into a large file. It is, of course, possible to automate the entire process, but this is not part of the scope of this repository.

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

The `jobCommand` should be executed as-is by your computing cluster's own job management and dispatch system. Please consult with your system administrator or person-in-charge if you have any doubts on submitting jobs to a queue.

## What about the Python scripts?

TL;DR: Leave them alone.

The python files are simply dependencies used by the bash scripts and **should not be executed out of context or order, and should never be used individually**.

## Visualization

The Rscript in the `/visualization` directory allows you use visualize data generated from the genetic map. Please note that all steps in the `/lib` directory need to be completed first: this Rscript requires certain files to be in place before you can proceed.

The Rscript will generate charts in PDF and PNG format in the respective RIL population folder, namely `Gifu-Burttii` and `Gifu-MG20` folders in the output directory. The Rscript **will need to be run twice**: and remember to swap the `pop1` and `pop2` variables around, so that both individual RIL population will be processed.
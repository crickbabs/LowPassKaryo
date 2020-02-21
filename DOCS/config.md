

## Pre-requisites 

You will need the following software tools to be installed on the system:

1. A job scheduler - the current code base includes a wrapper for ['slurm'](https://slurm.schedmd.com/documentation.html)
2. ['R'](https://cran.r-project.org/)
3. The Bioconductor package ['QDNASeq'](https://bioconductor.org/packages/release/bioc/html/QDNAseq.html)
4. ['FastQC'](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)
5. A quality and adapter trimming tool - the current code base includes a wrapper for [`Trimgalore`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
6. An alignment tool - the current code base includes a wrapper for ['bwa'](https://github.com/lh3/bwa) version 0.7.15-r1140
7. ['Samtools'](http://samtools.sourceforge.net/) - the current code base includes a wrapper for version 0.1.19-44428cd



## Set-up

Obtain a local copy of the pipeline by running the following command from the directory in which you would like to store the source code:

```bash
git clone https://github.com/crickbabs/LowPassKaryo
```


Next, edit the following scripts to customise to your own system:

1. ./LowPassKaryo.R - Set the variable "SITE.CONFIG.FILE" to point at your copy of config.R (which is to be found in the SRC directory)
2. ./SRC/config.R - This file should be edited to represent your own system. In particular the following variables should be set explicitly
CONTACT.ADDRESS : the email address of the first point of contact for people receiving the report files.
SOURCE.DIR: 	the location of the pipeline source files (almost certainly the directory that config.R is stored in)
GENOMES.LOOKUP.FILE:	the location of the genomes look up file -- see below for the appropriate format for this file
BINARY.CALLS: 	the system specific call for each binary used. If you use a module system you should include the module load statements here. 

If you choose to use an aligner or trimmer other than bwa/Trimmomatic respectively, then you will also need to provide an appropriate WRAPPED_do script and update the appropriate WRAPPED.FILES and VERSION.COMMANDS entries in config.R


If you use a scheduler other than slurm, you will also need to define appropriate build.submisison, get.id & check.for.id functions in the scheduler_submisison_functions.R script & update the assignment of 
build.scheduler.submisison, get.scheduler.id, check.queue.for.id, DEFAULT.SCHEDULER.OPTIONS and THREADS.REGEXP. 


##GENOMES.LOOKUP.FILE

This should be a tab delimited text file containing five columns with the following names
1. Tag:	This column should contain the "name" corresponding to this genome which will be included in the design file (e.g. Homo sapiens or mm10)
2. Description: The formal name of the genome including build and release if appropriate (e.g. GRCh38-r39 or hg19)
3. Species: The generic name for the species (e.g. Homo sapiens, Mus musculus)
4. Reference_Path: path to the reference sequence file to be passed to your aligner of choice.
5. QDNAseq_Annotated_Bins_Path: path to an RDat object for use by QDNASeq. The script build_qdnaseq_binfile.R can be used to generate objects in the appropriate format - see the script for more details

See genome_lookup_table_template.txt for an example of the format.









## Introduction

An ['R'](https://cran.r-project.org/) pipeline for automatic processing of low-pass whole genome sequencing data to detect copy number variation using the ['QDNASeq'](https://bioconductor.org/packages/release/bioc/html/QDNAseq.html) package.


## Pipeline summary

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc))
2. Adapter/Quality trimming ([`Trimgalore`][https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/])
3. Post trimming QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc))
4. Alignment (['bwa'][https://github.com/lh3/bwa] (version 0.7.15-r1140))
5. Sorting and indexing (['Samtools'][http://samtools.sourceforge.net/])
6. Copy number calling (['QDNASeq'][https://bioconductor.org/packages/release/bioc/html/QDNAseq.html])
7. Summary report generation (['R'][https://cran.r-project.org/])


##Typical workflow 

After low-pass whole genome sequencing of a number of samples, a typical workflow will involve 

1. Creating a design file associating each set of FastQ files with the appropriate sample, genome and annotation information. 
2. Passing this design file to the main LowPassKaryo_Wrapper.R script which will sanity check the parameters and then handle submission of procesing jobs to your HPC cluster/farm.
3. On sucessful completion, the pipeline will produce one pdf file containing QDNASeq copy number profiles for each species included in the processing run and an html report containing primary alignment QC metrics and recording the software versions used.


Details of the local configuration required to set up the pipeline and also instructions on how to subsequently run it may be found in the DOCS/ directory

1. [Local Config](DOCS/config.md)
2. [Usage](DOCS/usage.md)



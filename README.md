# LowPassKaryo
An R pipeline to process low pass whole genome sequencing and call copy number variation


Code will be made public in the near future.

## Pipeline summary

1. Raw read QC (FastQC)
2. (Optional) Adapter/Quality trimming (Trimgalore)
3. (Optional) Post trimming QC (FastQC)
4. Alignment (BWA)
5. Sorting and indexing (samtools)
6. Copy number calling (QDNASeq)
7. Summary report generation (R)

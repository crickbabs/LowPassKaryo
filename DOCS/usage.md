

## Useage


The pipeline may be started by passing  LowPassKaryo_Wrapper.R to Rscript with one or two command line parameters, e.g.

```bash
Rscript LowPassKaryo_Wrapper.R --help
```

The principal useage though will be to provide a design file:


```bash
Rscript LowPassKaryo_Wrapper.R --design lowpass_design_file.txt
```


The design file should be a tab separated text file consisting of two sections: 

A two column header containing
1. A name for the project/run
2. The directory containign the FastQ files
3. The base working directory - a sub-directory will be created using the project name and results and temporary files will be created here.
4. Optional Scheduler parameters - entries other than NA here will over-ride (i.e. completely replace) the existing defaults.
(Column 1 contains the variable tags Column 2 contains the values)


Following this is a four column table with one row per sample and the following columns
1. ##LIMS_ID - A unique identifier for this sample corresponding to the first underscore separated field of the fastq files names. e.g. for NGS1234_S26_L005_R1_001.fastq.gz the LIMS_ID would be NGS1234
2. ##SAMPLE_NAME - a descriptive name for the sample, e.g. WildType-Time0-Replicate1
3. ##GENOME_TAG - which genome should the sample be aligned to. The entries here should correspond to the "Tag" entries of the genome_lookup_table.txt
4. ##ANNO_BED - an optional bed file indicating regions of interest to be annotated in the logR plots. Set to NA to skip annotation.

A template design file indicating the expected format can be generated as follows

```bash
Rscript LowPassKaryo_Wrapper.R ----template-design
```


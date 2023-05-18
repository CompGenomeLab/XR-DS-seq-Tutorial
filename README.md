# XR-DS-seq-Tutorial
 GitHub repository tutorial for XR-Seq and Damage-Seq methodologies, providing code examples, datasets, and documentation for studying DNA damage and repair processes.

## Introduction to NGS

## NGS File Format

## Introduction to R, Python, Bash (whichever is necessary)

## Downloading reference genome and generating related files 

## Quality control

## Adaptor handling, mapping, quality trimming, and converting to bed

## Sorting, filtering, calculating length dist., filtering damage-seq samples by motif, producing bigwig files
### Sorting BED files
The BED files are sorted according to genomic coordinates of the reads.

> sort -k1,1 -k2,2n -k3,3n ${SAMPLE}.bed > ${SAMPLE}_sorted.bed

### Obtaining read length distribution and read count
Checking the 


## Simulating the sample reads

## Plotting length distribution, nucleotide enrichment, and bam correlations

## Running pipeline with snakemake

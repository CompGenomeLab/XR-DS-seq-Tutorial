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

    sort -k1,1 -k2,2n -k3,3n ${SAMPLE}.bed > ${SAMPLE}_sorted.bed

### Obtaining read length distribution and read count
In order see if our data contains the reads with lengths expexted from the XR-seq or. Damage-seq methods, read length distribution of the data should be extracted.

>

In addition, we count the reads in the BED files and use this count in the following steps. 

>

 ### Filtering Damage-seq data by motif

Damage-seq reads are expected to contain C and T nucleotides at certain position due to the nature of the NER and Damage-seq experimental protocol. Therefore, we use this as an additional filtering step to eliminate the reads that don’t fit this criteria.

>

### Generating Bigwig files

Bigwig file format includes representation of the distribution reads in each genomic window without respect to plus and minus strands. It is a compact file format that is required for many tools in the further analysis steps. 

The first step for generation of bigwig files from BED files is generating BedGraph files. Since the strands of the reads are not taken into account in the Bigwig and BedGraph files

>

Then, from these Bedgraph files, we can generate Bigwig files. 


## Simulating the sample reads
Because CPD and (6-4)PP damage types require certain nucleotides in certain positions, the genomic locations rich in adjacent CC, TC, CT or TT dinucleotides may be prone to receiving more UV damage while other regions that are poor in these dinucleotides receive less damage. Therefore, the sequence contents may bias our analysis results while comparing the damage formation or NER efficiency of two genomic regions. In order to eliminate the effect of  sequence content, we create synthetic sequencing data from the real Damage-seq and XR-seq data, which give us the expected damage counts and NER efficiencies, respectively, from the sequence content of the genomic areas of interest.

The read counts from the simulated Damage-seq and XR-seq data are then used to normalize our real Damage-seq and XR-seq data to eliminate the sequence content bias. 

[Boquila github link]

    boquila


## Plotting length distribution, nucleotide enrichment, and bam correlations

### Plotting read length distribution

We create a histogram of the read length distribution to clearly understand which read length is mostly found in our data. This information is used as a proof about the success of Damage-seq or XR-seq experiment in capturing the right reads. 

>


### Plotting nucleotide enrichment of the reads

Since the reads from Damage-seq and XR-seq data are expected to contain C and T nucleotides in certain positions, we plot the nucleotide enrichment in each position in reads. This also is used as a proof of data quality. 


### Bam correlations

Another way to assess the data quality is to compare samples and replicates. Replicates should be alike in terms of genomic distribution of their reads while, for example, CPD samples should be different than (6-4)PP samples. To check this, we plot bam correlations. 

Pearson corr??


### Plotting read distribution on genes

deepTools is an easy way to plot read distribution in certain genomic regions. This tool uses the Bigwig files to count the reads in genomic windows and requires the regions of interest in BED format. To plot a line graph of the read  distribution on genes:

    deeptools 

## Running pipeline with snakemake

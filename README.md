# XR-DS-seq-Tutorial
 GitHub repository tutorial for XR-Seq and Damage-Seq methodologies, providing code examples, datasets, and documentation for studying DNA damage and repair processes.

## Introduction to NGS

## NGS File Format

## Introduction to R, Python, Bash (whichever is necessary)

## Downloading reference genome and generating related files 

## Quality control

## Adaptor handling, mapping, quality trimming, and converting to bed

## Processing BED files and assessing data quality
### Sorting BED files
After obtaining the BED files, the reads in the BED files are sorted according to genomic coordinates.

    sort -k1,1 -k2,2n -k3,3n hela_xr_cpd.bed > hela_ds_cpd_sorted.bed
    sort -k1,1 -k2,2n -k3,3n hela_ds_cpd.bed > hela_xr_cpd_sorted.bed

### Filtering for chromosomes
Reads that are aligned to regions other than chromosomes 1-22 and X are filtered and the rest is kept for further analysis.

    grep "^chr" hela_ds_cpd_sorted.bed \ 
        grep -v -e "chrY" -e "chrM" > \
        hela_ds_cpd_sorted_chr.bed
    
    grep "^chr" hela_xr_cpd_sorted.bed \ 
        grep -v -e "chrY" -e "chrM" > \
        hela_xr_cpd_sorted_chr.bed

### Separating strands
The reads in the BED files are separated according to which strand they were mapped on.
```
    awk '{if($6=="+"){print}}' hela_ds_cpd_sorted_chr.bed > \
        hela_ds_cpd_sorted_chr_plus.bed
    awk '{if($6=="-"){print}}' hela_ds_cpd_sorted_chr.bed > \
        hela_ds_cpd_sorted_chr_minus.bed

    awk '{if($6=="+"){print}}' hela_xr_cpd_sorted_chr.bed > \
        hela_xr_cpd_sorted_chr_plus.bed
    awk '{if($6=="-"){print}}' hela_xr_cpd_sorted_chr.bed > \
        hela_xr_cpd_sorted_chr_minus.bed
```

### Centering the damage site in for Damage-seq reads
Due to the experimental protocol of Damage-seq, the damaged dipyrimidines are located at the two bases upstream of the reads. We use [bedtools flank](https://www.google.com/search?q=bedtools+flank&rlz=1C1GCEU_enTR1010TR1010&oq=bedtools&aqs=chrome.0.69i59j69i57j69i64j69i59l2j69i60l3.2002j0j7&sourceid=chrome&ie=UTF-8) and [bedtools slop](https://bedtools.readthedocs.io/en/latest/content/tools/slop.html) to obtain 10-nucleotide long read locations with the damaged nucleotides at 5th and 6th positions.
```
    bedtools flank -i  hela_ds_cpd_sorted_chr_plus.bed -l 6 -r 0 > \
    hela_ds_cpd_sorted_chr_plus_flank.bed
    bedtools flank -i hela_ds_cpd_sorted_chr_minus.bed -l 0 -r 6 > \
    hela_ds_cpd_sorted_chr_minus_flank.bed
```
```
    bedtools slop -i hela_ds_cpd_sorted_chr_plus_flank.bed \
        -l 0 -r 4 > \
        hela_ds_cpd_sorted_chr_plus_10.bed
    bedtools slop -i hela_ds_cpd_sorted_chr_minus_flank.bed \
        -l 4 -r 0 > \
        hela_ds_cpd_sorted_chr_minus_10.bed
```
```
    cat hela_ds_cpd_sorted_chr_plus_10.bed hela_ds_cpd_sorted_chr_minus_10.bed > \
        hela_ds_cpd_sorted_chr_10.bed
```

### Obtaining FASTA files from BED files
To convert BED files to FASTA format:
    
    bedtools getfasta \
        -fi GRCh38.fa \
        -bed hela_ds_cpd_sorted_chr_plus_10.bed \
        -fo hela_ds_cpd_sorted_plus_10.fa \
        -s

    bedtools getfasta \
        -fi GRCh38.fa \
        -bed hela_ds_cpd_sorted_chr_minus_10.bed \
        -fo hela_ds_cpd_sorted_minus_10.fa \
        -s

### Filtering Damage-seq data by motif
Damage-seq reads are expected to contain C and T nucleotides at certain position due to the nature of the NER and Damage-seq experimental protocol. Therefore, we use this as an additional filtering step to eliminate the reads that don’t fit this criteria.

    python3 scripts/fa2bedByChoosingReadMotifs.py hela_ds_cpd_sorted_plus_10.fa
    python3 scripts/fa2bedByChoosingReadMotifs.py hela_ds_cpd_sorted_minus_10.fa

    cat hela_ds_cpd_sorted_ds_dipyrimidines_plus.bed \   
        hela_ds_cpd_sorted_ds_dipyrimidines_minus.bed > \
        hela_ds_cpd_sorted_ds_dipyrimidines.bed

### Obtaining read length distribution and read count
In order see if our data contains the reads with lengths expected from the XR-seq or Damage-seq methods, read length distribution of the data should be extracted.

    awk '{{print $3-$2}}' hela_ds_cpd_sorted_ds_dipyrimidines.bed \
        sort -k1,1n \
        uniq -c \
        sed 's/\s\s*/ /g' \
        awk '{{print $2"\\t"$1}}' > \
        hela_ds_cpd_sorted_ds_dipyrimidines_ReadLengthDist.txt

    awk '{{print $3-$2}}' hela_xr_cpd_sorted_chr.bed \&
        sort -k1,1n \& 
        uniq -c \& 
        sed 's/\s\s*/ /g' \&
        awk '{{print $2"\\t"$1}}' > \
        hela_xr_cpd_sorted_chr_ReadLengthDist.txt

In addition, we count the reads in the BED files and use this count in the following steps. 

    grep -c "^" hela_ds_cpd_sorted_ds_dipyrimidines.bed > hela_ds_cpd_sorted_ds_dipyrimidines_readCount.txt
    grep -c "^" hela_xr_cpd_sorted_chr.bed > hela_xr_cpd_sorted_chr_readCount.txt

### Assessing dinucleotide content


    fa2kmerAbundanceTable.py 
### Generating BigWig files

BigWig file format includes representation of the distribution reads in each genomic window without respect to plus and minus strands. It is a compact file format that is required for many tools in the further analysis steps. 

The first step for generation of BigWig files from BED files is generating BedGraph files. Since the strands of the reads are not taken into account in the BigWig and BedGraph files

    bedtools genomecov -i hela_ds_cpd_sorted_ds_dipyrimidines_plus.bed \
        -g ref_genome/GRCh38.fa.fai \
        -bg -scale \
        $(cat hela_ds_cpd_sorted_ds_dipyrimidines_readCount.txt | \
        awk '{print 1000000/$1}') > \
        hela_ds_cpd_sorted_ds_dipyrimidines_plus.bdg
    
    bedtools genomecov -i hela_ds_cpd_sorted_ds_dipyrimidines_minus.bed \
        -g ref_genome/GRCh38.fa.fai \
        -bg -scale \
        $(cat hela_ds_cpd_sorted_ds_dipyrimidines_readCount.txt | \
        awk '{print 1000000/$1}') > \
        hela_ds_cpd_sorted_ds_dipyrimidines_minus.bdg

    bedtools genomecov -i hela_xr_cpd_sorted_chr_plus.bed \
        -g ref_genome/GRCh38.fa.fai \
        -bg -scale \
        $(cat hela_xr_cpd_sorted_chr_readCount.txt | \
        awk '{print 1000000/$1}') > \
        hela_xr_cpd_sorted_chr_plus.bdg        
    
    bedtools genomecov -i hela_xr_cpd_sorted_chr_minus.bed \
        -g ref_genome/GRCh38.fa.fai \
        -bg -scale \
        $(cat hela_xr_cpd_sorted_chr_readCount.txt | \
        awk '{print 1000000/$1}') > \
        hela_xr_cpd_sorted_chr_minus.bdg

Then, from these BedGraph files, we generate BigWig files. 

    sort -k1,1 -k2,2n hela_ds_cpd_sorted_ds_dipyrimidines_plus.bdg > \
        hela_ds_cpd_sorted_ds_dipyrimidines_plus_sorted.bdg
    sort -k1,1 -k2,2n hela_ds_cpd_sorted_ds_dipyrimidines_minus.bdg > \
        hela_ds_cpd_sorted_ds_dipyrimidines_minus_sorted.bdg
    
    sort -k1,1 -k2,2n hela_xr_cpd_sorted_chr_plus.bdg > \
        hela_xr_cpd_sorted_chr_plus_sorted.bdg
    sort -k1,1 -k2,2n hela_xr_cpd_sorted_chr_minus.bdg > \
        hela_xr_cpd_sorted_chr_minus_sorted.bdg
    
    bedGraphToBigWig hela_ds_cpd_sorted_ds_dipyrimidines_plus_sorted.bdg

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

deepTools is an easy way to plot read distribution in certain genomic regions. This tool uses the BigWig files to count the reads in genomic windows and requires the regions of interest in BED format. To plot a line graph of the read  distribution on genes:

    deeptools 

## Running pipeline with snakemake

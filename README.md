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

    sort -k1,1 -k2,2n -k3,3n results/hela_ds_cpd.bed > results/hela_ds_cpd_sorted.bed
    sort -k1,1 -k2,2n -k3,3n results/hela_xr_cpd.bed > results/hela_xr_cpd_sorted.bed

### Filtering for chromosomes
Reads that are aligned to regions other than chromosomes 1-22 and X are filtered and the rest is kept for further analysis.

    grep "^chr" results/hela_ds_cpd_sorted.bed | grep -v -e "chrY" -e "chrM" > results/hela_ds_cpd_sorted_chr.bed
    
    grep "^chr" results/hela_xr_cpd_sorted.bed | grep -v -e "chrY" -e "chrM" > results/hela_xr_cpd_sorted_chr.bed

### Separating strands
The reads in the BED files are separated according to which strand they were mapped on.
```
    awk '{if($6=="+"){print}}' results/hela_ds_cpd_sorted_chr.bed > results/hela_ds_cpd_sorted_chr_plus.bed
    awk '{if($6=="-"){print}}' results/hela_ds_cpd_sorted_chr.bed > results/hela_ds_cpd_sorted_chr_minus.bed

    awk '{if($6=="+"){print}}' results/hela_xr_cpd_sorted_chr.bed > results/hela_xr_cpd_sorted_chr_plus.bed
    awk '{if($6=="-"){print}}' results/hela_xr_cpd_sorted_chr.bed > results/hela_xr_cpd_sorted_chr_minus.bed
```

### Centering the damage site in for Damage-seq reads
Due to the experimental protocol of Damage-seq, the damaged dipyrimidines are located at the two bases upstream of the reads. We use [bedtools flank](https://www.google.com/search?q=bedtools+flank&rlz=1C1GCEU_enTR1010TR1010&oq=bedtools&aqs=chrome.0.69i59j69i57j69i64j69i59l2j69i60l3.2002j0j7&sourceid=chrome&ie=UTF-8) and [bedtools slop](https://bedtools.readthedocs.io/en/latest/content/tools/slop.html) to obtain 10-nucleotide long read locations with the damaged nucleotides at 5th and 6th positions.
```
    cut -f1,2 ref_genome/GRCh38.fa.fai > ref_genome/sizes.chrom

    bedtools flank -i  results/hela_ds_cpd_sorted_chr_plus.bed -g ref_genome/sizes.chrom -l 6 -r 0 > results/hela_ds_cpd_sorted_chr_plus_flank.bed
    bedtools flank -i results/hela_ds_cpd_sorted_chr_minus.bed -g ref_genome/sizes.chrom -l 0 -r 6 > results/hela_ds_cpd_sorted_chr_minus_flank.bed
```
```
    bedtools slop -i results/hela_ds_cpd_sorted_chr_plus_flank.bed -g ref_genome/sizes.chrom -l 0 -r 4 > results/hela_ds_cpd_sorted_chr_plus_10.bed
    bedtools slop -i results/hela_ds_cpd_sorted_chr_minus_flank.bed -g ref_genome/sizes.chrom -l 4 -r 0 > results/hela_ds_cpd_sorted_chr_minus_10.bed
```
```
    cat results/hela_ds_cpd_sorted_chr_plus_10.bed results/hela_ds_cpd_sorted_chr_minus_10.bed > results/hela_ds_cpd_sorted_chr_10.bed
```

### Obtaining FASTA files from BED files
To convert BED files to FASTA format:
    
    bedtools getfasta -fi ref_genome/GRCh38.fa -bed results/hela_ds_cpd_sorted_chr_plus_10.bed -fo results/hela_ds_cpd_sorted_plus_10.fa -s
    bedtools getfasta -fi ref_genome/GRCh38.fa -bed results/hela_ds_cpd_sorted_chr_minus_10.bed -fo results/hela_ds_cpd_sorted_minus_10.fa -s

### Filtering Damage-seq data by motif
Damage-seq reads are expected to contain C and T nucleotides at certain position due to the nature of the NER and Damage-seq experimental protocol. Therefore, we use this as an additional filtering step to eliminate the reads that don’t fit this criteria.

    python3 scripts/fa2bedByChoosingReadMotifs.py -i results/hela_ds_cpd_sorted_plus_10.fa -o results/hela_ds_cpd_sorted_ds_dipyrimidines_plus.bed -r '.{4}(c|t|C|T){2}.{4}'
    python3 scripts/fa2bedByChoosingReadMotifs.py -i results/hela_ds_cpd_sorted_minus_10.fa -o results/hela_ds_cpd_sorted_ds_dipyrimidines_minus.bed -r '.{4}(c|t|C|T){2}.{4}'

    cat results/hela_ds_cpd_sorted_ds_dipyrimidines_plus.bed results/hela_ds_cpd_sorted_ds_dipyrimidines_minus.bed > results/hela_ds_cpd_sorted_ds_dipyrimidines.bed

### Obtaining read length distribution and read count
In order see if our data contains the reads with lengths expected from the XR-seq or Damage-seq methods, read length distribution of the data should be extracted.

    awk '{print $3-$2}' results/hela_ds_cpd_sorted_chr.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}' > results/hela_ds_cpd_sorted_chr_ReadLengthDist.txt

    awk '{print $3-$2}' results/hela_xr_cpd_sorted_chr.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}' > results/hela_xr_cpd_sorted_chr_ReadLengthDist.txt

In addition, we count the reads in the BED files and use this count in the following steps. 

    grep -c "^" results/hela_ds_cpd_sorted_ds_dipyrimidines.bed > results/hela_ds_cpd_sorted_ds_dipyrimidines_readCount.txt
    grep -c "^" results/hela_xr_cpd_sorted_chr.bed > results/hela_xr_cpd_sorted_chr_readCount.txt

### Assessing nucleotide content
We check the reads for the enrichment of nucleotides and dinucleotides at specific positions.

    bedtools getfasta -fi ref_genome/GRCh38.fa -bed results/hela_ds_cpd_sorted_ds_dipyrimidines.bed -fo results/hela_ds_cpd_sorted_ds_dipyrimidines.fa -s

    bedtools getfasta -fi ref_genome/GRCh38.fa -bed results/hela_xr_cpd_sorted_chr.bed -fo results/hela_xr_cpd_sorted_chr.fa -s

    python3 scripts/fa2kmerAbundanceTable.py -i results/hela_ds_cpd_sorted_ds_dipyrimidines.fa -k 1 -o results/hela_ds_cpd_sorted_ds_dipyrimidines_nucleotideTable.txt

    python3 scripts/fa2kmerAbundanceTable.py -i results/hela_xr_cpd_sorted_chr.fa -k 1 -o results/hela_xr_cpd_sorted_chr_nucleotideTable.txt

### Generating BigWig files

BigWig file format includes representation of the distribution reads in each genomic window without respect to plus and minus strands. It is a compact file format that is required for many tools in the further analysis steps. 

The first step for generation of BigWig files from BED files is generating BedGraph files. Since the strands of the reads are not taken into account in the BigWig and BedGraph files

    conda install -c mvdbeek ucsc_tools

    bedtools genomecov -i results/hela_ds_cpd_sorted_ds_dipyrimidines_plus.bed -g ref_genome/GRCh38.fa.fai -bg -scale $(cat results/hela_ds_cpd_sorted_ds_dipyrimidines_readCount.txt | awk '{print 1000000/$1}') > results/hela_ds_cpd_sorted_ds_dipyrimidines_plus.bdg
    
    bedtools genomecov -i results/hela_ds_cpd_sorted_ds_dipyrimidines_minus.bed -g ref_genome/GRCh38.fa.fai -bg -scale $(cat results/hela_ds_cpd_sorted_ds_dipyrimidines_readCount.txt | awk '{print 1000000/$1}') > results/hela_ds_cpd_sorted_ds_dipyrimidines_minus.bdg

    bedtools genomecov -i results/hela_xr_cpd_sorted_chr_plus.bed -g ref_genome/GRCh38.fa.fai -bg -scale $(cat results/hela_xr_cpd_sorted_chr_readCount.txt | awk '{print 1000000/$1}') > results/hela_xr_cpd_sorted_chr_plus.bdg        
    
    bedtools genomecov -i results/hela_xr_cpd_sorted_chr_minus.bed -g ref_genome/GRCh38.fa.fai -bg -scale $(cat results/hela_xr_cpd_sorted_chr_readCount.txt | awk '{print 1000000/$1}') > results/hela_xr_cpd_sorted_chr_minus.bdg

Then, from these BedGraph files, we generate BigWig files. 

    sort -k1,1 -k2,2n results/hela_ds_cpd_sorted_ds_dipyrimidines_plus.bdg > results/hela_ds_cpd_sorted_ds_dipyrimidines_plus_sorted.bdg
    sort -k1,1 -k2,2n results/hela_ds_cpd_sorted_ds_dipyrimidines_minus.bdg > results/hela_ds_cpd_sorted_ds_dipyrimidines_minus_sorted.bdg
    
    sort -k1,1 -k2,2n results/hela_xr_cpd_sorted_chr_plus.bdg > results/hela_xr_cpd_sorted_chr_plus_sorted.bdg
    sort -k1,1 -k2,2n results/hela_xr_cpd_sorted_chr_minus.bdg > results/hela_xr_cpd_sorted_chr_minus_sorted.bdg
        
    bedGraphToBigWig results/hela_ds_cpd_sorted_ds_dipyrimidines_plus_sorted.bdg ref_genome/GRCh38.fa.fai results/hela_ds_cpd_sorted_ds_dipyrimidines_plus.bw
    bedGraphToBigWig results/hela_ds_cpd_sorted_ds_dipyrimidines_minus_sorted.bdg ref_genome/GRCh38.fa.fai results/hela_ds_cpd_sorted_ds_dipyrimidines_minus.bw

    bedGraphToBigWig results/hela_xr_cpd_sorted_chr_plus_sorted.bdg ref_genome/GRCh38.fa.fai results/hela_xr_cpd_sorted_chr_plus.bw
    bedGraphToBigWig results/hela_xr_cpd_sorted_chr_minus_sorted.bdg ref_genome/GRCh38.fa.fai results/hela_xr_cpd_sorted_chr_minus.bw

## Simulating the sample reads
Because CPD and (6-4)PP damage types require certain nucleotides in certain positions, the genomic locations rich in adjacent CC, TC, CT or TT dinucleotides may be prone to receiving more UV damage while other regions that are poor in these dinucleotides receive less damage. Therefore, the sequence contents may bias our analysis results while comparing the damage formation or NER efficiency of two genomic regions. In order to eliminate the effect of  sequence content, we create synthetic sequencing data from the real Damage-seq and XR-seq data, which give us the expected damage counts and NER efficiencies, respectively, from the sequence content of the genomic areas of interest. We use [Boquila](https://github.com/CompGenomeLab/boquila) to generate simulated data.

```
    conda install boquila -c bioconda 

    boquila --fasta results/hela_ds_cpd_sorted_ds_dipyrimidines.fa --bed results/hela_ds_cpd_sorted_ds_dipyrimidines_sim.bed --ref ref_genome/GRCh38.fa --regions ref_genome/GRCh38.ron --kmer 2 --seed 1 --sens 2 > results/hela_ds_cpd_sorted_ds_dipyrimidines_sim.fa

    boquila --fasta results/hela_xr_cpd_sorted_chr.fa --bed results/hela_xr_cpd_sorted_chr_sim.bed --ref ref_genome/GRCh38.fa --regions ref_genome/GRCh38.ron --kmer 2 --seed 1 --sens 2 > results/hela_xr_cpd_sorted_chr_sim.fa

```
    awk '{if($6=="+"){print}}' results/hela_ds_cpd_sorted_ds_dipyrimidines_sim.bed > results/hela_ds_cpd_sorted_ds_dipyrimidines_sim_plus.bed
    awk '{if($6=="-"){print}}' results/hela_ds_cpd_sorted_ds_dipyrimidines_sim.bed > results/hela_ds_cpd_sorted_ds_dipyrimidines_sim_minus.bed

    awk '{if($6=="+"){print}}' results/hela_xr_cpd_sorted_chr_sim.bed > results/hela_xr_cpd_sorted_chr_sim_plus.bed
    awk '{if($6=="-"){print}}' results/hela_xr_cpd_sorted_chr_sim.bed > results/hela_xr_cpd_sorted_chr_sim_minus.bed

The read counts from the simulated Damage-seq and XR-seq data are then used to normalize our real Damage-seq and XR-seq data to eliminate the sequence content bias. 
    
    grep -c '^' results/hela_xr_cpd_sorted_chr_sim.bed > results/hela_xr_cpd_sorted_chr_sim_readCount.txt

    bedtools genomecov -i results/hela_xr_cpd_sorted_chr_sim_plus.bed -g ref_genome/GRCh38.fa.fai -bg -scale $(cat results/hela_xr_cpd_sorted_chr_sim_readCount.txt | awk '{print 1000000/$1}') > results/hela_xr_cpd_sorted_chr_sim_plus.bdg
    
    bedtools genomecov -i results/hela_xr_cpd_sorted_chr_sim_minus.bed -g ref_genome/GRCh38.fa.fai -bg -scale $(cat results/hela_xr_cpd_sorted_chr_sim_readCount.txt | awk '{print 1000000/$1}') > results/hela_xr_cpd_sorted_chr_sim_minus.bdg

    sort -k1,1 -k2,2n results/hela_xr_cpd_sorted_chr_sim_plus.bdg > results/hela_xr_cpd_sorted_chr_sim_plus_sorted.bdg
    sort -k1,1 -k2,2n results/hela_xr_cpd_sorted_chr_sim_minus.bdg > results/hela_xr_cpd_sorted_chr_sim_minus_sorted.bdg

    bedGraphToBigWig results/hela_xr_cpd_sorted_chr_sim_plus_sorted.bdg ref_genome/GRCh38.fa.fai results/hela_xr_cpd_sorted_chr_sim_plus.bw
    bedGraphToBigWig results/hela_xr_cpd_sorted_chr_sim_minus_sorted.bdg ref_genome/GRCh38.fa.fai results/hela_xr_cpd_sorted_chr_sim_minus.bw

## Plotting length distribution, nucleotide enrichment, and bam correlations

### Plotting read length distribution

We create a histogram of the read length distribution to clearly understand which read length is mostly found in our data. This information is used as a proof about the success of Damage-seq or XR-seq experiment in capturing the right reads. 

    conda install r-rbokeh
    conda install -c conda-forge r-ggplot2
    conda install -c conda-forge r-dplyr

    R
    library(ggplot2)
    read_len <- read.table("results/hela_xr_cpd_sorted_chr_ReadLengthDist.txt")
    colnames(read_len) <- c("length", "counts")
    xrLenDistPlot <- ggplot(data = read_len, aes(x = length, y = counts)) +
        geom_bar(stat='identity') +
        xlab("Read length") +
        ylab("Read count")
    ggsave("results/xrLenDistPlot.png")

### Plotting nucleotide enrichment of the reads

Since the reads from Damage-seq and XR-seq data are expected to contain C and T nucleotides in certain positions, we plot the nucleotide enrichment in each position in reads. This also is used as a proof of data quality. 

    R
    library(ggplot2)
    nucl_content <- read.table("results/hela_xr_cpd_sorted_chr_nucleotideTable.txt")
    nucl_content_gathered <- gather(nucl_content, nucl, count, -kmer)
    nucl_content_gathered$kmer <- factor(nucl_content_gathered$kmer, levels = c("C", "T", "G", "A"))
    xrNuclPlot <- ggplot(data = nucl_content_gathered, aes(x = nucl, y = count, fill = kmer)) +
        geom_bar(position = "fill", stat='identity') +
        xlab("Position") +
        ylab("Count") +
        labs(fill = "") +
        scale_fill_manual(values = c("seagreen3","gray60", "steelblue2", "steelblue4"))
    ggsave("results/xrNuclPlot.png")


### Bam correlations

Another way to assess the data quality is to compare samples and replicates. Replicates should be alike in terms of genomic distribution of their reads while, for example, CPD samples should be different than (6-4)PP samples.

    conda install -c bioconda deeptools

    multiBamSummary bins --bamfiles results/hela_ds_cpd_cutadapt_sorted_dedup.bam results/hela_xr_cpd_cutadapt_sorted_dedup.bam --minMappingQuality 20 --outFileName results/bamSummary.npz
    
    plotCorrelation -in results/bamSummary.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o results/bamCorr.png \

### Plotting read distribution on genes

[deepTools](https://deeptools.readthedocs.io/en/develop/) is an easy way to plot read distribution in certain genomic regions. This tool uses the BigWig files to count the reads in genomic windows and requires the regions of interest in BED format. To plot a line graph of the read  distribution on genes:

    bigwigCompare --bigwig1 results/hela_xr_cpd_sorted_chr_plus.bw --bigwig2 results/hela_xr_cpd_sorted_chr_sim_plus.bw --operation ratio --outFileFormat bigwig --outFileName hela_xr_cpd_norm_sim_plus.bw

    computeMatrix scale-regions -S results/hela_xr_cpd_norm_sim_plus.bw -R ref_genome/GRCh38_genes.bed --outFileName results/hela_xr_cpd_norm_sim_plus_on_GRCh_genes_1kb_scaleregions_computeMatrix.out -b 1000 -a 1000 --smartLabels

    plotHeatmap --matrixFile results/hela_xr_cpd_norm_sim_plus_on_GRCh_genes_1kb_scaleregions_computeMatrix.out --outFileName hela_xr_cpd_norm_sim_plus_on_GRCh_genes_1kb_heatmap_k2.png --xAxisLabel "Position with respect to genes" --yAxisLabel "Genes" --kmeans 2 --heatmapWidth 10 --startLabel "Start" --endLabel "End"

## Running pipeline with snakemake


# Analysis Tutorial for XR-seq and Damage-seq

This tutorial focuses on the manual analysis of XR-seq and Damage-seq data from raw reads to genome-wide tracks and quality control plots. Use the sections below for step-by-step commands.

Prerequisites:
- A GitHub account (for obtaining materials and sharing results). Create one at [https://github.com/join](https://github.com/join).
- Linux environment is recommended. If you are not on Linux, please follow our WSL setup guide: [wsl_vscode_github.pdf](wsl_vscode_github.pdf).
- [Visual Studio Code](https://code.visualstudio.com/) is highly recommended (not mandatory) for the exercises and file navigation.
- Conda environments should be installed and ready-to-use (see [ENVIRONMENTS.md](ENVIRONMENTS.md)).

Outline:
- Concepts: XR vs DS signals, strand conventions, lesion/products.
- QC: FastQC/MultiQC.
- Adapter handling: cutadapt (XR trim-and-keep; DS discard-trim).
- Alignment and BAM ops: bowtie2, samtools, duplicate handling, MAPQ filtering.
- BED conversion; chromosome/strand filtering.
- Damage-seq motif filtering at lesion positions.
- Coverage and normalization: bedGraph/bigWig per strand; quick IGV review.
- QC and analyses: length distributions, nucleotide enrichment, replicate correlations, region-level profiles (deepTools).
- Advanced: simulation-based controls with boquila.

Requirements:
- Example FASTQs in `samples/`; outputs under `results/`.
- Reference genome in `ref_genome/` (instructions below use GRCh38 p14).

### Introduction to XR-seq

Excision Repair sequencing (XR-seq) captures the short oligomers excised by nucleotide excision repair (NER). Key features:

- Captures excision products produced by NER; signals reflect excision activity across the genome.
- Fragment lengths are protocol- and organism-dependent (e.g., often ~24–32 nt after trimming in many eukaryotes, shorter in some prokaryotes); signal is strand-aware.
- Interpreting XR-seq focuses on excision footprints and strand assignments rather than any single repair subpathway.
- Analysis goals: adapter trimming (retain trimmed reads), alignment with short-read-aware parameters, deduplication/quality filtering, convert to BED with correct strandedness, generate normalized coverage (bigWig), and compute QA metrics (length distribution, nucleotide enrichment, correlations).

### Introduction to Damage-seq

Damage-seq maps DNA lesions (e.g., CPD, (6-4)PP) by capturing polymerase-arrested sites. Key features:

- Reads report lesion positions; for UV photoproducts (CPD, (6-4)PP), lesions occur at dipyrimidines; motif filtering improves specificity.
- Lesion position is typically two bases upstream of the read’s 5′ end; downstream steps center windows accordingly before motif checks.
- Adapter handling is diagnostic: presence of adapter indicates no polymerase arrest (no lesion), so discard-trim is used to keep only lesion-containing reads.
- Analysis goals: adapter detection/removal (discard-trim), alignment, high-quality unique mappings, lesion-site coordinate definition, motif filtering, normalized coverage, and QA metrics (length distributions, motif/nucleotide enrichment).

### Further reading

- Hu, J., Li, W., Adebali, O., et al. Genome-wide mapping of nucleotide excision repair with XR-seq. Nature Protocols 14, 248–282 (2019). A detailed, step-by-step XR-seq laboratory and analysis protocol with practical considerations and preliminary analysis examples.  
  Link: https://www.nature.com/articles/s41596-018-0093-7

- Hu, J., Adar, S., Selby, C. P., Lieb, J. D. & Sancar, A. Genome-wide analysis of human global and transcription-coupled excision repair of UV damage at single-nucleotide resolution. Genes & Development 29, 948–960 (2015). Introduces single-nucleotide–resolution repair maps, separating global and transcription-coupled repair components in human cells.  
  Link: https://genesdev.cshlp.org/content/29/9/948

- Adar, S., Hu, J., Lieb, J. D. & Sancar, A. Genome-wide kinetics of DNA excision repair in relation to chromatin state and mutagenesis. PNAS (2016). Quantifies excision repair kinetics genome-wide and relates repair dynamics to chromatin features and mutational patterns.  
  Link: https://www.pnas.org/doi/10.1073/pnas.1614430113

- Li, W., Hu, J., et al. (2017) PNAS. Extends genome-wide mapping by integrating UV damage formation and repair dynamics to reveal determinants of repair heterogeneity.  
  Link: https://www.pnas.org/doi/10.1073/pnas.1706522114


## NGS File Format

This section explains the file formats that you will encounter troughout the tutorial. You can find more detail about each of these formats together with others from the [link](https://genome.ucsc.edu/FAQ/FAQformat.html).

- FASTQ:

  - The FASTQ format is used to store both the sequencing reads and their corresponding quality scores.
  - It consists of four lines per sequence read:
  - Line 1: Begins with a '@' symbol followed by a unique identifier for the read.
  - Line 2: Contains the actual nucleotide sequence of the read.
  - Line 3: Starts with a '+' symbol and may optionally contain the same unique identifier as Line 1.
  - Line 4: Contains the quality scores corresponding to the nucleotides in Line 2.
  - FASTQ files are commonly generated by sequencing platforms and are widely used as input for various NGS data analysis tools.

- FASTA:

  - The FASTA format is a simple and widely used text-based format for representing nucleotide or protein sequences.
  - Each sequence in a FASTA file consists of two parts:
    - A single-line description or identifier that starts with a '>' symbol, followed by a sequence identifier or description.
    - The sequence data itself, which can span one or more lines.
  - The sequence lines contain the actual nucleotide or amino acid symbols representing the sequence.
  - FASTA files can store multiple sequences, each with its own identifier and sequence data.
  - FASTA format is commonly used for storing and exchanging sequences, such as reference genomes, gene sequences, or protein sequences.
  - Many bioinformatics tools and databases accept FASTA files as input for various sequence analysis tasks, including alignment, motif discovery, and similarity searches.

- SAM/BAM:

  - SAM (Sequence Alignment/Map) is a text-based format used to store the alignment information of reads to a reference genome.
  - BAM (Binary Alignment/Map) is the binary equivalent of SAM, which is more compact and allows for faster data processing.
  - Both SAM and BAM files contain the same information, including read sequences, alignment positions, mapping qualities, and optional tags for additional metadata.
  - SAM/BAM files are crucial for downstream analysis tasks such as variant calling, differential expression analysis, and visualization in genome browsers.
  - SAM/BAM files can be generated using aligners like Bowtie2 or BWA.

- BED:

  - BED (Browser Extensible Data) is a widely used format for representing genomic intervals or features.
  - It consists of a tab-separated plain-text file with columns representing different attributes of each genomic feature.
  - The standard BED format includes columns for chromosome, start position, end position, feature name, and additional optional columns for annotations or scores.
  - BED files are commonly used for visualizing and analyzing genomic features such as gene coordinates, binding sites, peaks, or regulatory regions.
  - BED files can be generated from BAM files using tools like BEDTools or samtools.

- BEDGRAPH:

  - BEDGRAPH is a file format used to represent continuous numerical data across the genome, such as signal intensities, coverage, or scores.
  - Similar to BED files, BEDGRAPH files are plain-text files with columns representing genomic intervals and corresponding values.
  - The standard BEDGRAPH format includes columns for chromosome, start position, end position, and a numerical value representing the data for that interval.
  - BEDGRAPH files are typically used to visualize and analyze genome-wide data, such as ChIP-seq signal, DNA methylation levels, or RNA-seq coverage.
  - The values in a BEDGRAPH file can be positive or negative, representing different aspects of the data.
  - BEDGRAPH files can be generated from BAM files using tools like BEDTools or bedGraphToBigWig.

- BigWig:

  - BigWig is a binary file format designed for efficient storage and retrieval of large-scale numerical data across the genome.
  - BigWig files are created from BEDGRAPH files and provide a compressed and indexed representation of the data, allowing for fast random access and visualization.
  - BigWig files store the data in a binary format and include additional indexing information for quick retrieval of data within specific genomic regions.
  - They are commonly used for visualizing and analyzing genome-wide data in genome browsers or for performing quantitative analyses across different genomic regions.
  - BigWig files can be generated from BEDGRAPH files using tools like bedGraphToBigWig or through conversion from other formats like BAM or WIG.

## Create conda environments

See [ENVIRONMENTS.md](ENVIRONMENTS.md) for creating and switching between the step-specific environments used in this tutorial. After creating the environments, ensure working directories exist:

```bash
mkdir -p ref_genome/Bowtie2 results qc
```

## Downloading reference genome and generating related files

This section provides commands to download a reference genome (GRCh38) and generate related files.
For alignment and BAM operations, activate the mapping environment:

```bash
conda activate mapping
```

Next, you will download the reference genome.
Because you will analyze results/HeLa cells that are a type of human cancer cells (derived from cervical cancer cells), GRCh38 genome is chosen as the reference:

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.p14.genome.fa.gz -O ref_genome/GRCh38.fa.gz && gunzip -f ref_genome/GRCh38.fa.gz
```

After downloading the genome fasta file, you should generate the index files for bowtie2
which is crucial for efficiently reducing computational time and improving alignment quality.

```bash
bowtie2-build --threads 8 ref_genome/GRCh38.fa ref_genome/Bowtie2/genome_GRCh38
```

You will create another index file via samtools to retrieve sequence information based on genomic coordinates or sequence identifiers.

```bash
samtools faidx ref_genome/GRCh38.fa
```

Next, the index file created by samtools will be converted to a ron file (a different format of the index).
This file will be useful when we simulate our samples with boquila.

```bash
python3 scripts/idx2ron.py -i ref_genome/GRCh38.fa.fai -o ref_genome/GRCh38.ron -l genome_idex2ron.log
```

## Quality control

Here, you'll learn how to perform quality control on the NGS data using FastQC, a tool for assessing the quality of sequencing reads.
A nice [github page](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html) for the evaluation of FastQC results.

Activate the `fastqc` environment:

```bash
conda activate fastqc
```

Then you can create a directory for the output of FastQC with `mkdir qc/` command.
Lastly, you can run the command below to execute the tool for both Damage-seq and XR-seq samples:

```bash
fastqc -t 8 --outdir qc/ samples/hela_ds_cpd.fq 

fastqc -t 8 --outdir qc/ samples/hela_xr_cpd.fq 

multiqc qc/ -o qc/
```

## Adaptor handling

This section demonstrates how to handle adaptors in the NGS data using Cutadapt, a tool for trimming adaptor sequences and other contaminants. Activate the `trimming` environment:

```bash
conda activate trimming
```

At this stage, we will perform different tasks for Damage-seq and XR-seq.
In Damage-seq we want to discard all the reads with adaptors
since having the adaptor in Damage-seq reads mean that the read does not contain any damage.
In the case of XR-seq, we will only trim the adaptors from the reads and keep the trimmed ones.

```bash
cutadapt -j 8 -g GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT --discard-trimmed -o results/ds.trim.fq samples/hela_ds_cpd.fq

cutadapt -j 8 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG -o results/xr.trim.fq samples/hela_xr_cpd.fq
```

## Mapping, removing duplicates, quality trimming, and converting to bed

This section covers the steps involved in mapping the preprocessed reads to the reference genome using Bowtie2, converting the mapped reads to BAM format, and extracting BED files.

Initially we will align our reads to the reference genome using the prepared index files.
After that we will convert the output sam files to bam.

```bash
(bowtie2 --threads 8 --seed 1 --reorder -x ref_genome/Bowtie2/genome_GRCh38 -U results/ds.trim.fq -S results/ds.sam)

samtools view -Sbh -o results/ds.bam results/ds.sam

(bowtie2 --threads 8 --seed 1 --reorder -x ref_genome/Bowtie2/genome_GRCh38 -U results/xr.trim.fq -S results/xr.sam)

samtools view -Sbh -o results/xr.bam results/xr.sam
```

In the next part, we will remove the duplicate reads with picard. Activate the `mapping` environment (includes bowtie2, samtools, picard):

```bash
conda activate mapping
```

To run MarkDplicates command of picard (this command will remove the duplicates for us),
you need your files to be ordered by their header.
For that purpose, you should sort the files and then use picard.

```bash
samtools sort -o results/ds.sort.bam results/ds.bam -@ 8 -T results/

samtools sort -o results/xr.sort.bam results/xr.bam -@ 8 -T results/

(picard MarkDuplicates --REMOVE_DUPLICATES true --INPUT results/ds.sort.bam --TMP_DIR results/ --OUTPUT results/ds.dedup.bam --METRICS_FILE results/ds.dedup.metrics.txt)

(picard MarkDuplicates --REMOVE_DUPLICATES true --INPUT results/xr.sort.bam --TMP_DIR results/ --OUTPUT results/xr.dedup.bam --METRICS_FILE results/xr.dedup.metrics.txt)
```

Lastly, we will remove low quality reads (MAPQ score < 20) via samtools and
use bedtools to convert our bam files into bed format.

```bash
# Mapping env for samtools filtering
conda activate mapping
samtools index results/ds.dedup.bam
samtools view -q 20 -b results/ds.dedup.bam > results/ds.q20.bam

samtools index results/xr.dedup.bam
samtools view -q 20 -b results/xr.dedup.bam > results/xr.q20.bam

# Bedops env for BED conversion
conda activate bedops
bedtools bamtobed -i results/ds.q20.bam > results/ds.bed
bedtools bamtobed -i results/xr.q20.bam > results/xr.bed
```

## Sorting BED files

After obtaining the BED files, the reads in the BED files are sorted according to genomic coordinates.

```bash
sort -k1,1 -k2,2n -k3,3n results/ds.bed > results/ds.sorted.bed

sort -k1,1 -k2,2n -k3,3n results/xr.bed > results/xr.sorted.bed
```

## Filtering for chromosomes

Reads that are aligned to regions other than chromosomes 1-22 and X are filtered and the rest is kept for further analysis.

```bash
grep "^chr" results/ds.sorted.bed | grep -v -e "chrY" -e "chrM" > results/ds.sorted.chr.bed

grep "^chr" results/xr.sorted.bed | grep -v -e "chrY" -e "chrM" > results/xr.sorted.chr.bed
```

## Separating strands

The reads in the BED files are separated according to which strand they were mapped on.

```bash
    awk '{if($6=="+"){print}}' results/ds.sorted.chr.bed > results/ds.sorted.chr.plus.bed
    awk '{if($6=="-"){print}}' results/ds.sorted.chr.bed > results/ds.sorted.chr.minus.bed

    awk '{if($6=="+"){print}}' results/xr.sorted.chr.bed > results/xr.sorted.chr.plus.bed
    awk '{if($6=="-"){print}}' results/xr.sorted.chr.bed > results/xr.sorted.chr.minus.bed
```

## Centering the damage site in for Damage-seq reads

Due to the experimental protocol of Damage-seq, the damaged dipyrimidines are located at the two bases upstream of the reads. We use [bedtools flank](https://bedtools.readthedocs.io/en/latest/content/tools/flank.html) and [bedtools slop](https://bedtools.readthedocs.io/en/latest/content/tools/slop.html) to obtain 10-nucleotide long read locations with the damaged nucleotides at 5th and 6th positions.

```bash
    cut -f1,2 ref_genome/GRCh38.fa.fai > ref_genome/sizes.chrom

    bedtools flank -i  results/ds.sorted.chr.plus.bed -g ref_genome/sizes.chrom -l 6 -r 0 > results/ds.plus.flank.bed
    bedtools flank -i results/ds.sorted.chr.minus.bed -g ref_genome/sizes.chrom -l 0 -r 6 > results/ds.minus.flank.bed
```

```bash
    bedtools slop -i results/ds.plus.flank.bed -g ref_genome/sizes.chrom -l 0 -r 4 > results/ds.plus.10.bed
    bedtools slop -i results/ds.minus.flank.bed -g ref_genome/sizes.chrom -l 4 -r 0 > results/ds.minus.10.bed
```

```bash
    cat results/ds.plus.10.bed results/ds.minus.10.bed > results/ds.10.bed
```

### Obtaining FASTA files from BED files

To convert BED files to FASTA format:

```bash
    bedtools getfasta -fi ref_genome/GRCh38.fa -bed results/ds.plus.10.bed -fo results/ds.plus.10.fa -s
    bedtools getfasta -fi ref_genome/GRCh38.fa -bed results/ds.minus.10.bed -fo results/ds.minus.10.fa -s
```

### Filtering Damage-seq data by motif

Damage-seq reads are expected to contain C and T nucleotides at certain position due to the nature of the NER and Damage-seq experimental protocol. Therefore, we use this as an additional filtering step to eliminate the reads that don’t fit this criteria.

```bash
    python3 scripts/fa2bedByChoosingReadMotifs.py -i results/ds.plus.10.fa -o results/ds.dipy.plus.bed -r '.{4}(c|t|C|T){2}.{4}'
    python3 scripts/fa2bedByChoosingReadMotifs.py -i results/ds.minus.10.fa -o results/ds.dipy.minus.bed -r '.{4}(c|t|C|T){2}.{4}'

    cat results/ds.dipy.plus.bed results/ds.dipy.minus.bed > results/ds.dipy.bed
```

### Obtaining read length distribution and read count
In order see if our data contains the reads with lengths expected from the XR-seq or Damage-seq methods, read length distribution of the data should be extracted.

```bash
    awk '{print $3-$2}' results/ds.sorted.chr.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}' > results/ds.len.txt

    awk '{print $3-$2}' results/xr.sorted.chr.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}' > results/xr.len.txt
```

In addition, we count the reads in the BED files and use this count in the following steps. 

```bash
    grep -c "^" results/ds.dipy.bed > results/ds.dipy.count.txt
    grep -c "^" results/xr.sorted.chr.bed > results/xr.count.txt
```

### Assessing nucleotide content
We check the reads for the enrichment of nucleotides and dinucleotides at specific positions.

```bash
    bedtools getfasta -fi ref_genome/GRCh38.fa -bed results/ds.dipy.bed -fo results/ds.dipy.fa -s

    bedtools getfasta -fi ref_genome/GRCh38.fa -bed results/xr.sorted.chr.bed -fo results/xr.fa -s

    python3 scripts/fa2kmerAbundanceTable.py -i results/ds.dipy.fa -k 1 -o results/ds.dipy.nt.txt

    python3 scripts/fa2kmerAbundanceTable.py -i results/xr.fa -k 1 -o results/xr.nt.txt
```

### Generating BigWig files

BigWig file format includes representation of the distribution reads in each genomic window without respect to plus and minus strands. It is a compact file format that is required for many tools in the further analysis steps. 

The first step for generation of BigWig files from BED files is generating BedGraph files. Since the strands of the reads are not taken into account in the BigWig and BedGraph files

```bash
    # BedGraph generation with bedtools (bedops env)
    conda activate bedops
    bedtools genomecov -i results/ds.dipy.plus.bed -g ref_genome/GRCh38.fa.fai -bg -scale $(cat results/ds.dipy.count.txt | awk '{print 1000000/$1}') > results/ds.dipy.plus.bdg
    
    bedtools genomecov -i results/ds.dipy.minus.bed -g ref_genome/GRCh38.fa.fai -bg -scale $(cat results/ds.dipy.count.txt | awk '{print 1000000/$1}') > results/ds.dipy.minus.bdg

    bedtools genomecov -i results/xr.sorted.chr.plus.bed -g ref_genome/GRCh38.fa.fai -bg -scale $(cat results/xr.count.txt | awk '{print 1000000/$1}') > results/xr.plus.bdg        
    
    bedtools genomecov -i results/xr.sorted.chr.minus.bed -g ref_genome/GRCh38.fa.fai -bg -scale $(cat results/xr.count.txt | awk '{print 1000000/$1}') > results/xr.minus.bdg
```

Then, from these BedGraph files, we generate BigWig files. 

```bash
    # Sort BedGraphs (still in bedops env)
    sort -k1,1 -k2,2n results/ds.dipy.plus.bdg > results/ds.dipy.plus.sorted.bdg
    sort -k1,1 -k2,2n results/ds.dipy.minus.bdg > results/ds.dipy.minus.sorted.bdg
    
    sort -k1,1 -k2,2n results/xr.plus.bdg > results/xr.plus.sorted.bdg
    sort -k1,1 -k2,2n results/xr.minus.bdg > results/xr.minus.sorted.bdg
```

```bash
    # Convert to BigWig (ucsc env)
    conda activate ucsc
    bedGraphToBigWig results/ds.dipy.plus.sorted.bdg ref_genome/GRCh38.fa.fai results/ds.dipy.plus.bw
    bedGraphToBigWig results/ds.dipy.minus.sorted.bdg ref_genome/GRCh38.fa.fai results/ds.dipy.minus.bw

    bedGraphToBigWig results/xr.plus.sorted.bdg ref_genome/GRCh38.fa.fai results/xr.plus.bw
    bedGraphToBigWig results/xr.minus.sorted.bdg ref_genome/GRCh38.fa.fai results/xr.minus.bw
```

## Simulating the sample reads
Because CPD and (6-4)PP damage types require certain nucleotides in certain positions, the genomic locations rich in adjacent CC, TC, CT or TT dinucleotides may be prone to receiving more UV damage while other regions that are poor in these dinucleotides receive less damage. Therefore, the sequence contents may bias our analysis results while comparing the damage formation or NER efficiency of two genomic regions. In order to eliminate the effect of  sequence content, we create synthetic sequencing data from the real Damage-seq and XR-seq data, which give us the expected damage counts and NER efficiencies, respectively, from the sequence content of the genomic areas of interest. We use [Boquila](https://github.com/CompGenomeLab/boquila) to generate simulated data.

```bash
    conda activate boquila

    boquila --fasta results/ds.dipy.fa --bed results/ds.sim.bed --ref ref_genome/GRCh38.fa --regions ref_genome/GRCh38.ron --kmer 2 --seed 1 --sens 2 > results/ds.sim.fa

    boquila --fasta results/xr.fa --bed results/xr.sim.bed --ref ref_genome/GRCh38.fa --regions ref_genome/GRCh38.ron --kmer 2 --seed 1 --sens 2 > results/xr.sim.fa
```

```bash
    awk '{if($6=="+"){print}}' results/ds.sim.bed > results/ds.sim.plus.bed
    awk '{if($6=="-"){print}}' results/ds.sim.bed > results/ds.sim.minus.bed

    awk '{if($6=="+"){print}}' results/xr.sim.bed > results/xr.sim.plus.bed
    awk '{if($6=="-"){print}}' results/xr.sim.bed > results/xr.sim.minus.bed
```

The read counts from the simulated Damage-seq and XR-seq data are then used to normalize our real Damage-seq and XR-seq data to eliminate the sequence content bias.

```bash
    grep -c '^' results/xr.sim.bed > results/xr.sim.count.txt

    bedtools genomecov -i results/xr.sim.plus.bed -g ref_genome/GRCh38.fa.fai -bg -scale $(cat results/xr.sim.count.txt | awk '{print 1000000/$1}') > results/xr.sim.plus.bdg
    
    bedtools genomecov -i results/xr.sim.minus.bed -g ref_genome/GRCh38.fa.fai -bg -scale $(cat results/xr.sim.count.txt | awk '{print 1000000/$1}') > results/xr.sim.minus.bdg

    sort -k1,1 -k2,2n results/xr.sim.plus.bdg > results/xr.sim.plus.sorted.bdg
    sort -k1,1 -k2,2n results/xr.sim.minus.bdg > results/xr.sim.minus.sorted.bdg

    bedGraphToBigWig results/xr.sim.plus.sorted.bdg ref_genome/GRCh38.fa.fai results/xr.sim.plus.bw
    bedGraphToBigWig results/xr.sim.minus.sorted.bdg ref_genome/GRCh38.fa.fai results/xr.sim.minus.bw
```

## Plotting length distribution, nucleotide enrichment, and bam correlations

### Plotting read length distribution

We create a histogram of the read length distribution to clearly understand which read length is mostly found in our data. This information is used as a proof about the success of Damage-seq or XR-seq experiment in capturing the right reads. 

```bash
    conda activate rplots
```

``` R
    library(ggplot2)
    read_len <- read.table("results/xr.len.txt")
    colnames(read_len) <- c("length", "counts")
    xrLenDistPlot <- ggplot(data = read_len, aes(x = length, y = counts)) +
        geom_bar(stat='identity') +
        xlab("Read length") +
        ylab("Read count")
    ggsave("results/xrLenDistPlot.png")
```

### Plotting nucleotide enrichment of the reads

Since the reads from Damage-seq and XR-seq data are expected to contain C and T nucleotides in certain positions, we plot the nucleotide enrichment in each position in reads. This also is used as a proof of data quality. 

```R
library(ggplot2)
library(dplyr)
library(tidyr)

nucl_content <- read.table("results/xr.nt.txt", header=TRUE)
nucl_content_gathered <- gather(nucl_content, nucl, count, -kmer)
nucl_content_gathered$nucl <- factor(nucl_content_gathered$nucl,
    levels = c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22","X23","X24","X25","X26","X27"),
    labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27"))
nucl_content_gathered$kmer <- factor(nucl_content_gathered$kmer, levels = c("C","T","G","A"))

xrNuclPlot <- ggplot(data = nucl_content_gathered, aes(x = nucl, y = count, fill = kmer)) +
    geom_bar(position = "fill", stat='identity') +
    xlab("Position") +
    ylab("Count") +
    labs(fill = "") +
    scale_fill_manual(values = c("seagreen3","gray60","steelblue2","steelblue4"))
ggsave("results/xrNuclPlot.png", xrNuclPlot, width = 8, height = 4, dpi = 150)
```

### Bam correlations

Another way to assess the data quality is to compare samples and replicates. Replicates should be alike in terms of genomic distribution of their reads while, for example, CPD samples should be different than (6-4)PP samples.

```bash
    conda activate deeptools

    multiBamSummary bins --bamfiles results/ds.dedup.bam results/xr.dedup.bam --minMappingQuality 20 --outFileName results/bamSummary.npz
    
    plotCorrelation -in results/bamSummary.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o results/bamCorr.png
```

### Plotting read distribution on genes

[deepTools](https://deeptools.readthedocs.io/en/develop/) is an easy way to plot read distribution in certain genomic regions. This tool uses the BigWig files to count the reads in genomic windows and requires the regions of interest in BED format. To plot a line graph of the read  distribution on genes:

```bash
    bigwigCompare --bigwig1 results/xr.plus.bw --bigwig2 results/xr.sim.plus.bw --operation ratio --outFileFormat bigwig --outFileName results/xr.norm_sim.plus.bw

    computeMatrix scale-regions -S results/xr.norm_sim.plus.bw -R ref_genome/GRCh38_genes.bed --outFileName results/xr.norm_sim.plus.on_genes.computeMatrix.out -b 1000 -a 1000 --smartLabels

    plotHeatmap --matrixFile results/xr.norm_sim.plus.on_genes.computeMatrix.out --outFileName results/xr.norm_sim.plus.on_genes.heatmap_k2.png --xAxisLabel "Position with respect to genes" --yAxisLabel "Genes" --kmeans 2 --heatmapWidth 10 --startLabel "Start" --endLabel "End"
```

## Automated workflow

If you prefer an automated version of the steps in this tutorial, see the Snakemake workflow:

- Repository: [CompGenomeLab/xr-ds-seq-snakemake](https://github.com/CompGenomeLab/xr-ds-seq-snakemake)

It orchestrates genome preparation, QC, adapter handling, alignment, BED conversions, Damage-seq motif filtering, coverage track generation, simulations, and downstream summaries. Use it when you need a fully reproducible, end-to-end run with minimal manual intervention.


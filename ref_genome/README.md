# ref_genome/

Reference genome and indexes for the tutorial. Expected files:

- `Celegans.fa` — C. elegans (WBcel235/ce11) reference FASTA
- `Celegans.fa.fai` — FASTA index (samtools)
- `Celegans.ron` — regions index (from `scripts/idx2ron.py`) for boquila
- `sizes.chrom` — two-column chromosome sizes (derived from `.fai`)
- `Bowtie2/` — bowtie2 index files with prefix `genome_Celegans`

You can regenerate these using the commands in the main README.

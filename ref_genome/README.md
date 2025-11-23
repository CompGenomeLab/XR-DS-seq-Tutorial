# ref_genome/

Reference genome and indexes for the tutorial. Expected files (WBcel235/ce11):

- `celeg_genes.bed` — gene annotations in BED format aligned to this FASTA
- `GCF_000002985.6_WBcel235_genomic.fna` — reference FASTA (gunzipped from `.fna.gz`)
- `GCF_000002985.6_WBcel235_genomic.fna.fai` — FASTA index (samtools)
- `WBcel235.ron` — regions index (from `scripts/idx2ron.py`) for boquila
- `sizes.chrom` — two-column chromosome sizes (derived from `.fai`)
- `Bowtie2/` — bowtie2 index files with prefix `genome_Celegans`

You can regenerate these using the commands in the main README.

# scripts/

Helper scripts used in the tutorial:

- `fa2bedByChoosingReadMotifs.py` — filter FASTA reads by regex motif and output BED.
- `fa2kmerAbundanceTable.py` — compute k-mer abundance tables for FASTA sequences.
- `fasta.py`, `sequence.py` — FASTA/sequence utilities.
- `idx2ron.py` — convert `.fai` to `.ron` region file for boquila.
- `plot_kmer_content.R` — plot nucleotide (k=1) and dinucleotide (k=2) content from k-mer tables. Supports optional dinucleotide filter via `-f`.

These are invoked from the main README steps. See script `--help` for usage details.

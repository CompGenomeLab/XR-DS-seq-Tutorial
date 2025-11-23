# samples/

Place the input FASTQ files here. This tutorial assumes:

- `ce_ds.fastq.gz` — Damage-seq reads (C. elegans, wild-type)
- `ce_xr.fastq.gz` — XR-seq reads (C. elegans, wild-type)

Notes:

- Files are gzipped; tools in this tutorial accept `.fastq.gz` directly.
- If you use different filenames, update commands accordingly.

Source:

- Article: Kose, C., Azgari, C., Lindsey-Boltz, L. A., Adebali, O., & Sancar, A. (2025). Chromatin context shapes DNA damage formation and nucleotide excision repair dynamics in Caenorhabditis elegans. Nucleic Acids Research, 53(20), gkaf1080.
- GEO accessions:
  - Damage‑seq example (WT L1 0 h (6‑4)PP Rep1): [GSM8595180](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8595180)
  - XR‑seq example (WT L1 5 min (6‑4)PP Rep1): [GSM8595202](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8595202)

Subsampling:

- The provided tutorial FASTQs are 1,000,000‑read random subsamples of the above datasets, generated with `seqtk` (seed 42).

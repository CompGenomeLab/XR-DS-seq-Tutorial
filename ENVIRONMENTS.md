# Tutorial Conda Environments

This guide provides multiple small conda environments for different stages of the XR-seq and Damage-seq manual analysis. Keeping tools isolated reduces solver conflicts and makes it easy to switch between steps.

Initially, you will create and activate a conda environment to set up the necessary dependencies for running the tutorial.
If you haven't downloaded conda yet, you can find the instructions in the [link](https://conda.io/projects/conda/en/stable/user-guide/install/download.html).

Before creating any environment:

```bash
conda config --set channel_priority strict
conda config --add channels conda-forge
conda config --add channels bioconda
```

Note for macOS Apple Silicon (arm64): some bio-tools may not have native arm64 builds. If you hit resolution issues, use Rosetta, micromamba, or containers for those specific steps. The main tutorial remains the same.

## Environment list

- fastqc: FastQC, MultiQC
- trimming: cutadapt (adapter handling)
- mapping: bowtie2, samtools, picard (alignment, BAM ops, dedup)
- bedops: bedtools (BED conversion and processing)
- ucsc: ucsc_tools (bedGraphToBigWig and utilities)
- deeptools: deepTools (correlations, matrix/plots)
- rplots: R packages for plotting (ggplot2, tidyr, dplyr)
- boquila: boquila (simulation-based controls)

Environments are defined in `envs/*.yaml`. Create them once and activate as needed.

```bash
# Create all environments (run once)
conda env create -f envs/fastqc.yaml
conda env create -f envs/trimming.yaml
conda env create -f envs/mapping.yaml
conda env create -f envs/bedops.yaml
conda env create -f envs/ucsc.yaml
conda env create -f envs/deeptools.yaml
conda env create -f envs/rplots.yaml
conda env create -f envs/boquila.yaml
```

## Switching environments

Activate the environment for the step you are running, then deactivate when moving to the next step.

```bash
conda activate fastqc
# ... run FastQC and MultiQC ...
conda deactivate

conda activate trimming
# ... run cutadapt ...
conda deactivate
```

## Mapping tutorial steps to environments

- QC (FastQC/MultiQC): `fastqc`
  - Commands (from README): `fastqc ...`, `multiqc ...`

- Adapter handling (cutadapt): `trimming`
  - Commands: `cutadapt ...`

- Alignment and BAM processing: `mapping`
  - Commands: `bowtie2 ...`, `samtools ...`, `picard MarkDuplicates ...`

- BED conversion and BED ops: `bedops`
  - Commands: `bedtools bamtobed ...`, `bedtools sort/flank/slop/getfasta ...`

- Coverage conversion utilities: `ucsc`
  - Commands: `bedGraphToBigWig ...`

- QA and genome-wide analyses: `deeptools`
  - Commands: `multiBamSummary`, `plotCorrelation`, `computeMatrix`, `plotHeatmap`, `bigwigCompare`

- Plotting in R: `rplots`
  - Commands: R scripts to generate length distributions and nucleotide enrichment plots

- Simulation-based controls: `boquila`
  - Commands: `boquila ...`

## Verifying tools

After activating an environment, verify the primary tool is available:

```bash
conda activate mapping
bowtie2 --version
samtools --version
picard MarkDuplicates --version 2>/dev/null || echo "Picard is available"
conda deactivate
```

## Updating and removing environments

```bash
# Update packages within an env
conda activate deeptools
conda update --all -y
conda deactivate

# Remove an env if no longer needed
conda env remove -n deeptools
```

## Troubleshooting

- Solver conflicts: try creating the specific environment alone, or use micromamba for faster solves.
- Package not found on arm64: try `conda config --env --set subdir osx-64` inside that env before install, or use a container for that step.
- Path issues after switching envs: make sure to `conda deactivate` before `conda activate` another.



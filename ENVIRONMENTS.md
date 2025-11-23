# Tutorial Conda Environments

This guide provides multiple small conda environments for different stages of the XR-seq and Damage-seq manual analysis. Keeping tools isolated reduces solver conflicts and makes it easy to switch between steps.

Initially, you will create and activate a conda environment to set up the necessary dependencies for running the tutorial.
If you haven't downloaded conda yet, see the [Conda installation guide](https://conda.io/projects/conda/en/stable/user-guide/install/download.html).

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

Alternatively, a single command to create all environments (bash/zsh):

```bash
for i in fastqc trimming mapping bedops ucsc deeptools rplots boquila; do
  conda env create -f "envs/${i}.yaml"
done
```

Verify environments were created:

```bash
conda env list
# You should see at least:
# fastqc, trimming, mapping, bedops, ucsc, deeptools, rplots, boquila
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

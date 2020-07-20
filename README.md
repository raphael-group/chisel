# CHISEL <br/> <sub><u>C</u>opy-number <u>H</u>aplotype <u>I</u>nference in <u>S</u>ingle-cell by <u>E</u>volutionary <u>L</u>inks</sub> #

CHISEL is an algorithm to infer allele- and haplotype-specififc copy numbers in individual cells from low-coverage single-cell DNA sequencing data.
The full description of the algorithm and its application on published cancer datasets are described in

[Simone Zaccaria and Ben Raphael, 2019](https://doi.org/10.1101/837195)

The results of CHISEL on the published low-coverage single-cell DNA sequencing data of cancer patients and all the related analyses are available at

[CHISEL data](https://github.com/raphael-group/chisel-data)

This repository includes detailed instructions for installation and requirements, demos with related tutorial of different CHISEL applitcations, a list of current issues, and contacts.
The CHISEL repository is currently in a preliminary release and improved versions are released frequently.
During this stage, please keep checking for updates.

## Contents ##

1. [Overview](#overview)
    - [Algorithm](#algorithm)
    - [Software](#software)
2. [Setup](#setup)
    - [Standard installation](#standard)
    - [Full automatic installation](#automatic)
    - [Custom installation](#custom)
    - [Additional notes](#additionalnotes)
3. [Usage](#usage)
    - [Required data](#requireddata)
    - [System requirements](#requirements)
    - [Commands](#commands)
    - [Demos](#demos)
    - [Reccomendations and quality control](#reccomendations)
4. [Development](#development)
    - [Recent and important updates](#updates)
    - [Current issues](#currentissues)
5. [Contacts](#contacts)

## Overview
<a name="overview"></a>

### Algorithm
<a name="algorithm"></a>

![](doc/chisel-cartoon.png "CHISEL algorithm")

**The CHISEL algorithm.** **(A)** CHISEL computes RDRs and BAFs in low-coverage (<0.05X per cell) single-cell DNA sequencing data (top left). Read counts from 2000 individual cells (rows) in 5Mb genomic bins (columns) across three chromosomes (grey rectangles in first row) are shown. For each bin in each cell, CHISEL computes the RDR (top) by normalizing the observed read counts. CHISEL computes the BAF in each bin and cell (bottom) by first performing referenced-based phasing of germline SNPs in 50kb haplotype blocks (magenta and green) and then phasing all these blocks jointly across all cells. **(B)** CHISEL clusters RDRs and BAFs globally along the genome and jointly across all cells resulting here in 5 clusters of genomic bins (red, blue, purple, yellow, and grey) with distinct copy-number states. **(C)** CHISEL infers a pair of allele-specific copy numbers for each cluster by determining whether the allele-specific copy numbers of the largest balanced (BAF~0.5) cluster are equal to *{1, 1}* (diploid), *{2, 2}* (tetraploid), or are higher ploidy. **(D)** CHISEL infers haplotype-specific copy numbers *(a, b)* by phasing the allele-specific copy numbers consistently across all cells. **(E)** CHISEL clusters tumor cells into clones according to their haplotype-specific copy numbers. Here, a diploid clone (light gray) and two tumor clones (red and blue) are obtained. A phylogenetic tree describes the evolution of these clones. Somatic single-nucleotide variants (SNVs) are derived from pseudo-bulk samples and placed on the branches of the tree.

### Software
<a name="software"></a>

The current implementation of CHISEL is on `python2.7` and provides different commands (in `bin`) to atuomatically execute the different features/pipelines of CHISEL. See details in [Usage](#usage).

## Setup
<a name="setup"></a>

[CHISEL](https://bioconda.github.io/recipes/chisel/README.html) is distributed as a [bioconda](https://bioconda.github.io/) package.

1. [Standard installation](#standard)
2. [Full automatic installation](#automatic)
3. [Custom installation](#custom)
4. [Additional notes](#additionalnotes)

### Standard installation
<a name="standard"></a>

The recommended installation method for CHISEL is using `conda`.
First, install `conda`, a standard package manager, for example using the compact [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or the complete [Anaconda](https://www.anaconda.com/).
Moreover, `bioconda` [requires](https://bioconda.github.io/user/index.html) to set the following channels in this exact order:
```shell
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Next, `chisel` can be installed with the following one-time command:
```shell
conda install chisel
```

This command installs `chisel` in the base conda environment, however best practices generally recommend, especially in computing servers, to keep distinct packages indepedent.
Therefore, it is preferable to use the following command to install `chisel`:
```shell
conda create -n chisel chisel
```

In this case, `chisel` must be activated before every session with the following command
```shell
source activate chisel
```

Note that all the `conda` and `activate` commands above are located in `${CONDA_HOME}/bin/`.

### Full automatic installation
<a name="automatic"></a>

An automatic BASH script that fully installs `conda` through Miniconda and installs `chisel` on it is provided in [install_full.sh](install_full.sh): simply run `bash install_full.sh` from the CHISEL's home `${CHISEL_HOME}`.
Note that this procedures installs `conda` in `${CHISEL_HOME}/conda` and the following command (fully output at the end of the installation) must be executed before executing `chisel` in every new session:
```shell
source ${CHISEL_HOME}/conda/bin/activate chisel
```

### Custom installation
<a name="custom"></a>

CHISEL is written in `python2.7` and requires few standard python packages and two additional standard softwares.
Once all the requirements have been succesfully installed, `chisel` can be manually installed from `${CHISEL_HOME}` as
```shell
python setup.py install
```

#### > Python packages

CHISEL depends on the following standard python packages, which must be available in the python environment where the user runs CHISEL.

| Package | Tested version | Comments |
|---------|----------------|----------|
| [numpy](https://numpy.org/) | 1.16.1 | Efficient scientific computations |
| [scipy](https://www.scipy.org/) | 1.2.1 | Efficient mathematical functions and methods |
| [pandas](https://pandas.pydata.org/) | 0.20.1 | Dataframe management |
| [matplotlib](https://matplotlib.org/) | 2.0.2 | Basic plotting utilities |
| [seaborn](https://seaborn.pydata.org/) | 0.7.1 | Advanced plotting utilities |

#### > Additional software

CHISEL also requires few standard additional softwares, which must be included in `PATH` (e.g. by executing `export PATH=${SAMTOOLS_HOME}/bin/:${PATH}`)

| Software | Tested version | Comments |
|----------|----------------|----------|
| [AWK](https://en.wikipedia.org/wiki/AWK) | GNU Awk 4.0.2 | Scripting language available by default on most Unix-like OSs and with specific implementation available for any other OS |
| [SAMtools and BCFtools](http://www.htslib.org/doc/)  | 1.9 | Suite of programs for interacting with high-throughput sequencing data |

### Additional notes:
<a name="additionalnotes"></a>

- To be to run the conda commands `conda` and `activate` without prefix, add their path full path to `PATH` by
```shell
export PATH=${CONDA_HOME}/bin/:${PATH}
```

- `bioconda` [requires](https://bioconda.github.io/user/index.html) does not support Windows; however, Windows uses can use either [WSL](https://docs.microsoft.com/en-us/windows/wsl/) or the [custom installation](#custom).

## Usage
<a name="usage"></a>

1. [Required data](#requireddata)
2. [System requirements](#requirements)
3. [Commands](#commands)
4. [Demos](#demos)
5. [Reccomendations and quality control](#reccomendations)

### Required data
<a name="requireddata"></a>

CHISEL requires 4 input data:

1. A **single-cell barcoded BAM** file containing aligned DNA sequencing reads from single cells of a single patient. The BAM file must be indexed and sorted (see [SAMtools](http://www.htslib.org/workflow/#mapping_to_variant) and [GATK Best Practices](https://software.broadinstitute.org/gatk/best-practices/)). Reads must be labelled by a barcode, which is string `PREFIX:[A,C,G,T]+` composed of two parts: (1) a prefix `PREFIX`; and (2) the actual barcode `[A,C,G,T]+` which is a string of letters `A,C,G,T` of arbitrary length which uniquely identifies a single cell. Current CHISEL implementation requires the `PREFIX` to be `CB:Z:` according to the 10X Genomics [format](https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/output/bam). However, the upcoming version of CHISEL will support any PREFIX to accommodate any different technology. Note that the barcode `PREFIX:[A,C,G,T]+` can be either an independent tab-separated field of the SAM alignments or can be incorporated in the name of the sequencing read.

2. The **reference human genome** used to aligned the sequencing reads in the single-cell barcoded BAM. The most-used human reference genome are available at [GRC](https://www.ncbi.nlm.nih.gov/grc/human) or [UCSC](http://hgdownload.cse.ucsc.edu/downloads.html#human). Moreover, the reference genome `${REF}.fa` must be index by SAMtools (i.e. `samtools faidx ${REF}.fa`) and a dictionary must be created (`i.e. samtools dict ${REF}.fa > ${REF}.dict`) in the same folder as `${REF}`.

3. A **matched-normal BAM** file containing sequencing reads from a matched-normal sample from the same patient. When this sample is not available, user can use the corresponding CHISEL command `chisel-pseudonormal.py` to form a **pseudo** matched-normal sample by extracting diploid normal cells. Note that this approach is successfull only when normal diploid cells are present in the barcoded BAM file.

4. A **VCF file with phased germline SNPs** present in the matched-normal sample. CHISEL only requires a VCF with the phase `0|1` or `1|0` within the record of heterozygous phased SNPs and any reference-based phasing method can be used. Here, we recommend a very easy two step procedure: (1) use BCFtools to call germline SNPs from the matched-normal sample (as described [here](http://samtools.github.io/bcftools/howtos/variant-calling.html)); (2) use Eagle2 through the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!pages/home) to phase germline SNPs with respect to the reference panel HRC. Alternatively, Eagle2 can be used [locally](https://data.broadinstitute.org/alkesgroup/Eagle/), and other methods as [SHAPEIT](http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html) or other panels as 1000 Genomes Phase 3 can be applied, too. Note that HRC currently supports only reference `hg19`, therefore user can use [LiftoverVcf from Picard tools](http://broadinstitute.github.io/picard/) or similar methods to convert genomic coordinates to the required `hg19` build, when a different reference genome has been used to align sequencing reads. Note that phased germline SNPs can be also provided in `POS` format, i.e. tab-separated dataframe with field `CHROMOSOME  POSITION  PHASE` where phase is either `0|1` or `1|0`.

### System requirements
<a name="requirements"></a>

CHISEL is highly parallelized in order to make efficient the extensive computations needed for counting sequencing reads from BAM files, often splitting independent computations to parallel processes. We recommend executing CHISEL on multi-processing computing machines. Note that CHISEL automatically uses the maximum number of available cores, however this can be controlled with option `-j` of every command. The minimum system requirements for running the available demos or running every CHISEL command from partial results currently are:
- CPU with >2 virtual cores
- >4GB of RAM
However, when the complete CHISEL pipeline is applied to huge BAM files (>200GB) for millions of germline SNPs across thousands of single cells, CHISEL requirements are limited by the existing requirements of SAMtools and BCFtools, requiring higher amounts of RAM (e.g. >100GB of RAM) and a higher number of parallel cores (e.g. >12 cores). In general, the higher the number of cores is with a sufficient amount of RAM,  the faster the execution will be.

### CHISEL Commands
<a name="commands"></a>

CHISEL offers different commands to run either the entire pipeline with all the steps or only some specific steps. In particular, the latter commands are useful when user wants to re-run some specific steps by varying some of the default parameters or for quality control (see below).
Every command can be run directly when CHISEL has been correctly installed.

| Command | Description | Required input | Output |
|---------|-------------|----------------|--------|
| [`chisel`](man/chisel.md) | Running the complete CHISEL pipeline | The 4 required input data | [Final results and plots](doc/chisel.md) |
| [`chisel-calling`](man/chisel-calling.md) | Running from the inference of copy numbers | Computed RDRs and BAFs | [Final results and plots](doc/chisel-calling.md) |
| [`chisel-cloning`](man/chisel-cloning.md) | Running from the identification of clones | Inferred copy numbers | [Final results and plots](doc/chisel-cloning.md) |
| [`chisel-plotting`](man/chisel-plotting.md) | Running plot generation | Inferred copy numbers and clones | [Final plots](doc/chisel-plotting.md) |
| [`chisel-pseudonormal`](man/chisel-pseudonormal.md) | Extracting diploid cells to form a pseudo matched-normal sample | A barcoded BAM file containing diploid cells and the corresponding reference genome | [A BAM file to be used as a matched-normal sample](doc/chisel-pseudonormal.md) |

Click on the name of each command to obtain a description of all the available parameters.

### Demos
<a name="demos"></a>

Each demo is an exemplary and guided execution of a CHISEL command on available data. The demos are meant to illustrate the complete CHISEL pipeline and the usage of each specific command. In particular, the specific commands allow the user to explore the results when varying the main parameters (without the need to re-execute the entire pipeline), allowing to deal with different features of distinct datasets and to obtain the best-quality results (see recommendations below). Each demo is simultaneously a guided description of the entire example and a BASH script which can be directly executed to run the complete demo. As such, the user can both read the guided description as a web page and run the same script to execute the demo. At this time the following demos are available (more demos will be available soon):

| Demo | Description |
|------|-------------|
| [completeE](demos/completeE/demo-completeE.sh) | Demo of the complete CHISEL pipeline on tumor section E of breast cancer patient S0 |
| [callingE](demos/callingE/demo-callingE.sh) | Demo of `chisel-calling.py` command to re-run pipeline from the inference of copy numbers on tumor section E of breast cancer patient S0 |
| [cloningE](demos/cloningE/demo-cloningE.sh) | Demo of `chisel-cloning.py` command to re-run pipeline from the identification of clones on tumor section E of breast cancer patient S0 |
| [plottingE](demos/plottingE/demo-plottingE.sh) | Demo of `chisel-plotting.py` command to re-run plot generation on tumor section E of breast cancer patient S0 |
| [pseudonormal](demos/pseudonormal/demo-pseudonormal.sh) | Demo of `chisel-pseudonormal.py` command to generate a pseudo matched-normal sample by extracting diploid cells from a barcoded BAM file |

### Reccomendations and quality control
<a name="reccomendations"></a>

The following recommendations guide the user in the process of quality control for the final results. In order to deal with datasets with different features, these recommendations help the user in investigating solutions obtained with different parameters. In particular, some of these guides are especially helpful when analyzing datasets with non-standard features and with noisy or high-variance sequencing data, e.g. due to low number of cells, different sequencing coverages, low number of phased germline SNPs.

| Recommendation | Description |
|----------------|-------------|
| [Ploidy selection](guides/ploidy.md) | Adjusting model-selection parameters to accurately infer ploidy in noisy and high-variance datasets |
| [Identification of clones](guides/clones.md) | Interpreting the identified clones and explore alternative solutions |
| [Global clustering](guides/clustering.md) | Interpreting the global clustering of RDRs and BAFs, and exploring alternative solutions |

## Development
<a name="development"></a>

### Recent and important updates
<a name="updates"></a>

- **[20-Jul-2020]** CHISEL has been deposited and approved in [Bioconda](https://bioconda.github.io/recipes/chisel/README.html); check the new streamlined installation procedure.
- **[22-Jan-2020]** This version introduced two important updates:
    1. CHISEL now use the inferred clones to correct the inferred allele- and haplotype-specific copy numbers, removing noise and outlying errors. As such, CHISEL generates plots with both uncorrected and corrected copy numbers;
    2. CHISEL now has a **sensitivity**, which can be varied through the command `chisel-cloning.py` to enable the accurate inference of cell ploidy even in the case of noisy or high-variance data (e.g. low number of cells with low coverage).
- **[20-Jan-2020]** This version introduced the previously missing command `chisel-pseudonormal.py`. This command can be used to extract sequencing reads from diploid cells to form a pseudo matched-normal sample, which is required by the full CHISEL pipeline. A corresponding demo has also been introduced.

### Current issues
<a name="currentissues"></a>

CHISEL is in active development, please report any issue or question as this could help the development and improvement of CHISEL. Known issues with current version are reported here below.

## Contacts
<a name="contacts"></a>

CHISEL's repository is actively maintained by [Simone Zaccaria](https://simozacca.github.io/), currently a Postodoctoral Research Associate at Princeton University in the research group of prof. Ben Raphael.

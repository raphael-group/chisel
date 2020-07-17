# Demo for WGS data from a cancer patient
: ex: set ft=markdown ;:<<'```shell' #

The following CHISEL demo represents a guided example of the complete CHISEL pipeline starting from an exemplary barcoded [BAM file](https://doi.org/10.5281/zenodo.3950299) publicly available. From this directory, simply run this file through BASH as a standard script to run the complete demo. The demo can also be considered as a guided example of a complete execution and is correspondingly commented.

## Requirements and set up

The demo requires that CHISEL has been succesfully installed with conda. If the custom installation was used, please make sure that you can succesfully run the command `chisel` as well as the required `samtools`, `bcftools`, and `awk`. The demo includes the downloading of all the required files and will terminate in <20 minutes on machine with minimum requirements satisfied.

We gurantee that the running directory in the same directory of the demo and we remove previous results.

```shell
cd $( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )
rm -rf rdr/ baf/ combo/ calls/ clones/ plots/
:<<'```shell' # Ignore this line
```

We also ask the demo to terminate in case of errors and to print a trace of the execution by the following commands
```shell
set -e
set -o xtrace
PS4='[\t]'
:<<'```shell' # Ignore this line
```

## Downloading of data

The demo auomatically downloads the required barcoded single-cell and matched-normal BAM files in `data` folder.

```shell
# Creating data folder
mkdir -p data

# Downloading tumor barcoded BAM file
curl -L 'https://zenodo.org/record/3950299/files/cells.bam?download=1' > data/cells.bam
curl -L 'https://zenodo.org/record/3950299/files/cells.bam.bai?download=1' > data/cells.bam.bai
export TUM="data/cells.bam"

# Downloading matched-normal BAM file
curl -L 'https://zenodo.org/record/3950299/files/normal.bam?download=1' > data/normal.bam
curl -L 'https://zenodo.org/record/3950299/files/normal.bam.bai?download=1' > data/normal.bam.bai
export NOR="data/normal.bam"
:<<'```shell' # Ignore this line
```

Next the corresponding reference genome is downloaded and unpacked

```shell
curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz | gzip -d > data/hg19.fa
samtools faidx data/hg19.fa
samtools dict data/hg19.fa > data/hg19.dict
export REF="data/hg19.fa"
export DIC="data/hg19.dict"
:<<'```shell' # Ignore this line
```

Last, we download the pre-computed list of phased germline SNPs. Note that differently from the [one](https://github.com/raphael-group/chisel-data/raw/master/demos/completeE/phased.HRC.vcf.gz) obtained through the reccommended instructions (i.e. using BCFtools to call germline SNPs and Eagle2 throught the Michigan Imputation Serverve with HRC panel to phase the SNPs) this file only contains the lables `0|1` or `1|0` for every SNP, which is the minimum requirement for CHISEL.

```shell
curl -L 'https://zenodo.org/record/3950299/files/phases.tsv?download=1' > data/phased.tsv
export PHA="data/phased.tsv"
:<<'```shell' # Ignore this line
```

## Run CHISEL

We now run the complete pipeline of CHISEL with the corresponding command `chisel`.

```shell
chisel -t ${TUM} -n ${NOR} -r ${REF} -l ${PHA} --seed 12
exit $?
```

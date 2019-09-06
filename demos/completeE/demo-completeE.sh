# Demo for WGS data from a cancer patient
: ex: set ft=markdown ;:<<'```shell' #

The following CHISEL demo represents a guided example of the complete CHISEL pipeline starting from the barcoded [BAM file](https://support.10xgenomics.com/single-cell-dna/datasets/1.0.0/breast_tissue_E_2k) publicly available from 10X Genomics archive and obtained through 10X Chromium Single Cell CNV Solution for section E of a breast tumor. Simply run this file through BASH as a standard script to run the complete demo. The demo can also be considered as a guided example of a complete execution and is correspondingly commented.

## Requirements and set up

The demo requires that CHISEL has been succesfully installed, such that the python environment called by the command `python2.7` has the required packages, and both `samtools` and `awk` are available in `${PATH}`.

```shell
export CHISEL_HOME="../../" # This is CHISEL home by default, update if needed
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

The demo auomatically downloads the required barcoded single-cell and matched-normal BAM files from 10X Genomics archive through the following commands in `data` folder.

```shell
# Creating data folder
mkdir -p data

# Downloading barcoded single-cell BAM of breast tumor section E
wget -N -c http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-dna/1.0.0/breast_tissue_E_2k/breast_tissue_E_2k_possorted_bam.bam -P data/
wget -N -c http://cf.10xgenomics.com/samples/cell-dna/1.0.0/breast_tissue_E_2k/breast_tissue_E_2k_possorted_bam.bam.bai -P data/
export TUM="data/breast_tissue_E_2k_possorted_bam.bam"

# Downloading matched-normal BAM file as section A
wget -N -c http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-dna/1.0.0/breast_tissue_A_2k/breast_tissue_A_2k_possorted_bam.bam -P data/
wget -N -c http://cf.10xgenomics.com/samples/cell-dna/1.0.0/breast_tissue_A_2k/breast_tissue_A_2k_possorted_bam.bam.bai -P data/
export NOR="data/breast_tissue_A_2k_possorted_bam.bam"
:<<'```shell' # Ignore this line
```

Next the corresponding reference genome is downloaded and unpacked

```shell
export REF="data/refdata-GRCh38-2.1.0/fasta/genome.fa"
export DIC="data/refdata-GRCh38-2.1.0/fasta/genome.dict"
if [[ ! -f "${REF}" || ! -f "${DIC}" ]]; then
    wget -N -c http://cf.10xgenomics.com/supp/genome/refdata-GRCh38-2.1.0.tar.gz -P data/
    tar -xzvf data/refdata-GRCh38-2.1.0.tar.gz -C data/
    rm -f data/refdata-GRCh38-2.1.0.tar.gz
fi
:<<'```shell' # Ignore this line
```

Last, we download the pre-computed VCF with phased SNPs; the VCF has been computed following the reccommended instructions, using BCFtools to call germline SNPs and Eagle2 throught the Michigan Imputation Serverve with HRC panel to phase the SNPs.

```shell
wget -N -c https://github.com/raphael-group/chisel-data/raw/master/demos/completeE/phased.HRC.vcf.gz -P data/
gzip -d data/phased.HRC.vcf.gz
export PHA="data/phased.HRC.vcf"
:<<'```shell' # Ignore this line
```

## Run CHISEL

We now run the complete pipeline of CHISEL with the corresponding command `chisel`.

```shell
python2.7 ${CHISEL_HOME}/bin/chisel.py -t ${TUM} -n ${NOR} -r ${REF} -l ${PHA} --seed 25
exit $?
```

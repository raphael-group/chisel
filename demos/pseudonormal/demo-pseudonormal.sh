# Demo for generating pseudo matched-normal sample
: ex: set ft=markdown ;:<<'```shell' #

The following CHISEL demo represents a guided example of the CHISEL command for generating a pseudo matched-normal sample starting from the barcoded [BAM file](https://support.10xgenomics.com/single-cell-dna/datasets/1.0.0/breast_tissue_A_2k) publicly available from 10X Genomics archive and obtained through 10X Chromium Single Cell CNV Solution for section A of a breast tumor. Simply run this file through BASH as a standard script to run the complete demo. The demo can also be considered as a guided example of a complete execution and is correspondingly commented.

## Requirements and set up

The demo requires that CHISEL has been succesfully installed with conda. The demo includes the downloading of all the required files and will terminate in <20 minutes on machine with minimum requirements satisfied.

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

The demo auomatically downloads the required barcoded single-cell BAM file from 10X Genomics archive through the following commands in `data` folder.

```shell
# Creating data folder
mkdir -p data

# Downloading matched-normal BAM file as section A
curl -L http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-dna/1.0.0/breast_tissue_A_2k/breast_tissue_A_2k_possorted_bam.bam > data/breast_tissue_A_2k_possorted_bam.bam
curl -L http://cf.10xgenomics.com/samples/cell-dna/1.0.0/breast_tissue_A_2k/breast_tissue_A_2k_possorted_bam.bam.bai > data/breast_tissue_A_2k_possorted_bam.bam.bai
export BAM="data/breast_tissue_A_2k_possorted_bam.bam"
:<<'```shell' # Ignore this line
```

Last, the corresponding reference genome is downloaded and unpacked

```shell
export REF="data/refdata-GRCh38-2.1.0/fasta/genome.fa"
export DIC="data/refdata-GRCh38-2.1.0/fasta/genome.dict"
if [[ ! -f "${REF}" || ! -f "${DIC}" ]]; then
    curl -L http://cf.10xgenomics.com/supp/genome/refdata-GRCh38-2.1.0.tar.gz > data/refdata-GRCh38-2.1.0.tar.gz
    tar -xzvf data/refdata-GRCh38-2.1.0.tar.gz -C data/
    rm -f data/refdata-GRCh38-2.1.0.tar.gz
fi
:<<'```shell' # Ignore this line
```

## Run CHISEL

We now run the command `chisel-pseudonormal.py` of CHISEL for generating a pseudo mathched-normal sample by extracting the sequencing reads from diploid cells in the provided barcoded BAM file `${BAM}`.
Specifically, we are required to specify the reference genome `${REF}` and we use the default values of all parameters.
By default, temporary files and the sorted and indexed output BAM `pseudonormal.bam` will be generated in the current directory.

```shell
chisel_pseudonormal ${BAM} -r ${REF}
exit $?
```

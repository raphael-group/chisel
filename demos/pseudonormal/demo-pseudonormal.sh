# Demo for generating pseudo matched-normal sample
: ex: set ft=markdown ;:<<'```shell' #

The following CHISEL demo represents a guided example of the CHISEL command for generating a pseudo matched-normal sample starting from an exemplary barcoded [BAM file](https://doi.org/10.5281/zenodo.3952985) publicly available. Simply run this file through BASH as a standard script to run the complete demo. The demo can also be considered as a guided example of a complete execution and is correspondingly commented.

## Requirements and set up

The demo requires that CHISEL has been succesfully installed with conda. If the custom installation was used, please make sure that you can succesfully run the command `chisel_pseudonormal` as well as the required `samtools`, `bcftools`, and `awk`. The demo includes the downloading of all the required files and will terminate in <20 minutes on machine with minimum requirements satisfied.

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

# Downloading tumor barcoded BAM file
echo "Downloading tumor barcoded BAM file from Zenodo, please be patient as downloading time may vary."
curl -L https://zenodo.org/record/3952985/files/cells.bam?download=1 > data/cells.bam
curl -L https://zenodo.org/record/3952985/files/cells.bam.bai?download=1 > data/cells.bam.bai
export BAM="data/cells.bam"
:<<'```shell' # Ignore this line
```

Last, the corresponding reference genome is downloaded and unpacked

```shell
echo "Downloading human reference genome, please be patient as downloading time may vary."
curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz | gzip -d > data/hg19.fa
samtools faidx data/hg19.fa
samtools dict data/hg19.fa > data/hg19.dict
export REF="data/hg19.fa"
export DIC="data/hg19.dict"
:<<'```shell' # Ignore this line
```

## Run CHISEL

We now run the command `chisel_pseudonormal` of CHISEL for generating a pseudo mathched-normal sample by extracting the sequencing reads from diploid cells in the provided barcoded BAM file `${BAM}`.
Specifically, we are required to specify the reference genome `${REF}` and we use the default values of all parameters.
By default, temporary files and the sorted and indexed output BAM `pseudonormal.bam` will be generated in the current directory.

```shell
chisel_pseudonormal ${BAM} -r ${REF}
exit $?
```

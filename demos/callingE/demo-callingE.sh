# Demo for WGS data from a cancer patient
: ex: set ft=markdown ;:<<'```shell' #

The following CHISEL demo represents a guided example of the CHISEL pipeline starting from the computed RDRs and BAFs (typically the file `combo.tsv` in the folder `combo`) for tumor section E of breast cancer patient S0. Simply run this file through BASH as a standard script to run the complete demo. The demo represent a guided example for the command `chisel-calling` which allows to re-run the inference of copy numbers and can be used to try different parameters, especially related to the inference of tumor ploidy.

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

The demo auomatically downloads the required RDRs and BAFs already computed by the complete CHISEL pipeline in `data` folder.

```shell
# Creating data folder
mkdir -p data

# Downloading RDRs and BAFs computed by CHISEL for tumor section E
curl -L https://github.com/raphael-group/chisel-data/raw/master/demos/callingE/combo.tsv.gz > data/combo.tsv.gz
gzip -df data/combo.tsv.gz
export INPUT="data/combo.tsv"
:<<'```shell' # Ignore this line
```

## Run CHISEL

We now run the command CHISEL command that starts from the inference of copy numbers from RDRs and BAFs.

```shell
chisel_calling ${INPUT} --seed 25
exit $?
```

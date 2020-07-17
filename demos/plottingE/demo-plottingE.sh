# Demo for WGS data from a cancer patient
: ex: set ft=markdown ;:<<'```shell' #

The following CHISEL demo represents a guided example of the CHISEL pipeline starting from the inferred copy numbers (typically the file `calls.tsv` in the folder `calls`) and identified clones (typically the file `mapping.tsv` in the folder `clones`) for tumor section E of breast cancer patient S0, and thus produces the corresponding plots. The demo represent a guided example for the command `chisel-plotting` which allows to re-run the plot generation and can be used to try different parameters to obtain the best format for the results.

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

The demo auomatically downloads the required inferred copy numbers already computed by the complete CHISEL pipeline in `data` folder.

```shell
# Creating data folder
mkdir -p data

# Downloading copy numbers and clones inferred by CHISEL for tumor section E
wget -N https://github.com/raphael-group/chisel-data/raw/master/demos/cloningE/calls.tsv.gz -P data/
gzip -df data/calls.tsv.gz
export INPUT="data/calls.tsv"

wget -N https://github.com/raphael-group/chisel-data/raw/master/demos/plottingE/mapping.tsv.gz -P data/
gzip -df data/mapping.tsv.gz
export MAPP="data/mapping.tsv"
:<<'```shell' # Ignore this line
```

## Run CHISEL

We now run the command CHISEL command that starts from the inference of copy numbers from RDRs and BAFs.

```shell
chisel-plotting ${INPUT} -m ${MAPP}
exit $?
```

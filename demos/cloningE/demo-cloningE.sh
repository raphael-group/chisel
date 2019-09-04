# Demo for WGS data from a cancer patient
: ex: set ft=markdown ;:<<'```shell' #

The following CHISEL demo represents a guided example of the CHISEL pipeline starting from the inferred copy numbers (typically the file `calls.tsv` in the folder `calls`) for tumor section E of breast cancer patient S0, and thus identifies the clones and produces the corresponding plots. Simply run this file through BASH as a standard script to run the complete demo. The demo represent a guided example for the command `chisel-cloning.py` which allows to re-run the inference of clones and can be used to try different parameters to explore different solutions and clustering of cells..

## Requirements and set up

The demo requires that CHISEL has been succesfully installed, such that the python environment called by the command `python2.7` has the required packages.
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

The demo auomatically downloads the required inferred copy numbers already computed by the complete CHISEL pipeline in `data` folder.

```shell
# Creating data folder
mkdir -p data

# Downloading copy numbers inferred by CHISEL for tumor section E

export INPUT="data/calls.tsv"
:<<'```shell' # Ignore this line
```

## Run CHISEL

We now run the command CHISEL command that starts from the inference of copy numbers from RDRs and BAFs.

```shell
python2.7 ${CHISEL_HOME}/bin/chisel-cloning.py ${INPUT} --seed 25
exit $?
```

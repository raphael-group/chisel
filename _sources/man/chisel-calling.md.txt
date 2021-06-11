```shell
usage: chisel-calling.py [-h] [-x RUNDIR] [-A SENSITIVITY] [-P MAXPLOIDY]
                         [-K UPPERK] [--seed SEED] [-j JOBS]
                         [INPUT]

CHISEL command to re-run the inference of allele- and haplotype-specific copy
numbers, cell clustering, and plotting. This steps starts from estimated RDRs
and BAFs.

positional arguments:
  INPUT                 Input file with combined RDR and BAF per bin and per
                        cell (default: combo/combo.tsv)

optional arguments:
  -h, --help            show this help message and exit
  -x RUNDIR, --rundir RUNDIR
                        Running directory (default: current directory)
  -A SENSITIVITY, --sensitivity SENSITIVITY
                        Sensitivity of model selection for ploidy (default: 1,
                        increase this parameter to lower sensitivity to noisy
                        data, adjust this value (e.g. 2, 4, ..., 10, ...) to
                        better deal with high-variance data (e.g. low
                        coverage, small number of cells, low number of phased
                        SNPs, etc...)
  -P MAXPLOIDY, --maxploidy MAXPLOIDY
                        Maximum total copy number to consider for balanced
                        cluster (default: 4, corresponding to a WGD)
  -K UPPERK, --upperk UPPERK
                        Maximum number of bin clusters (default: 100, use 0 to
                        consider maximum number of clusters)
  --seed SEED           Random seed for replication (default: None)
  -j JOBS, --jobs JOBS  Number of parallele jobs to use (default: equal to
                        number of available processors)
```
```shell
usage: chisel_bedding [-h] [-x RUNDIR] [--rawcalls] [--noextending] [-j JOBS]
                      [INPUT]

CHISEL command to generate a BED file for each cell with the corresponding
CHISEL's results.

positional arguments:
  INPUT                 Input file with inferred copy numbers (default:
                        calls/calls.tsv)

optional arguments:
  -h, --help            show this help message and exit
  -x RUNDIR, --rundir RUNDIR
                        Running directory (default: current directory)
  --rawcalls            Use raw copy numbers instead of consensus corrected
                        ones (default: False)
  --noextending         Merge consecutive bins only if they are neighboring
                        (default: False, segments are extended to fill gaps)
  -j JOBS, --jobs JOBS  Number of parallele jobs to use (default: equal to
                        number of available processors)
```

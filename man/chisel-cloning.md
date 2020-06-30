usage: chisel-cloning.py [-h] [-x RUNDIR] [-f MAXDIFF] [-s MINSIZE]
                         [-r REFINEMENT] [--seed SEED]
                         [INPUT]

CHISEL command to run the pipeline starting from inferred copy numbers.

positional arguments:
  INPUT                 Input file with combined RDR and BAF per bin and per
                        cell

optional arguments:
  -h, --help            show this help message and exit
  -x RUNDIR, --rundir RUNDIR
                        Running directory (default: current directory)
  -f MAXDIFF, --maxdiff MAXDIFF
                        Maximum haplotype-specific distance between the genome
                        of cells in the same clone (default: 0.06, when -1 is
                        chosen the maximum cluster method of SciPy is used)
  -s MINSIZE, --minsize MINSIZE
                        Minimum number of cells in a subpopulation to define a
                        clone (default: 14)
  -r REFINEMENT, --refinement REFINEMENT
                        Maximum difference to assign noisy cells to the
                        closest clone (default: 0.0, note that 1.0 can be used
                        to force the assigment of all cells)
  --seed SEED           Random seed for replication (default: None)

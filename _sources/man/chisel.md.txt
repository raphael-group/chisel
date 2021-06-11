```shell
usage: chisel.py [-h] [-x RUNDIR] -t TUMOR -n NORMAL -r REFERENCE -l
                 LISTPHASED [-b SIZE] [-k BLOCKSIZE] [-c CHROMOSOMES]
                 [-m MINREADS] [-p MAXPLOIDY] [-K UPPERK]
                 [--bcftools BCFTOOLS] [--samtools SAMTOOLS]
                 [--cellprefix CELLPREFIX] [--cellsuffix CELLSUFFIX]
                 [--seed SEED] [-j JOBS]

CHISEL command to run the complete pipeline starting from the 4 required data:
(1) Barcoded single-cell BAM; (2) Matched-normal BAM; (3) Reference genome;
(4) Phased VCF.

optional arguments:
  -h, --help            show this help message and exit
  -x RUNDIR, --rundir RUNDIR
                        Running directory (default: current directory)
  -t TUMOR, --tumor TUMOR
                        Barcoded single-cell BAM file
  -n NORMAL, --normal NORMAL
                        Matched-normal BAM file
  -r REFERENCE, --reference REFERENCE
                        Reference genome
  -l LISTPHASED, --listphased LISTPHASED
                        Phased SNPs file (lines of heterozygous germline SNPs
                        must contain either 0|1 or 1|0)
  -b SIZE, --size SIZE  Bin size, with or without "kb" or "Mb"
  -k BLOCKSIZE, --blocksize BLOCKSIZE
                        Size of the haplotype blocks (default: 50kb, use 0 to
                        disable)
  -c CHROMOSOMES, --chromosomes CHROMOSOMES
                        Space-separeted list of chromosomes between apices
                        (default: "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8
                        chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17
                        chr18 chr19 chr20 chr21 chr22")
  -m MINREADS, --minreads MINREADS
                        Minimum number total reads to select cells (default:
                        100000)
  -p MAXPLOIDY, --maxploidy MAXPLOIDY
                        Maximum total copy number to consider for balanced
                        cluster (default: 4, corresponding to a WGD)
  -K UPPERK, --upperk UPPERK
                        Maximum number of bin clusters (default: 100, use 0 to
                        consider maximum number of clusters)
  --bcftools BCFTOOLS   Path to the directory to "bcftools" executable,
                        required in default mode (default: bcftools is
                        directly called as it is in user $PATH)
  --samtools SAMTOOLS   Path to the directory to "samtools" executable,
                        required in default mode (default: samtools is
                        directly called as it is in user $PATH)
  --cellprefix CELLPREFIX
                        Prefix of cell barcode field in SAM format (default:
                        CB:Z:)
  --cellsuffix CELLSUFFIX
                        Suffix of cell barcode field in SAM format (default:
                        none)
  --seed SEED           Random seed for replication (default: None)
  -j JOBS, --jobs JOBS  Number of parallele jobs to use (default: equal to
                        number of available processors)
```

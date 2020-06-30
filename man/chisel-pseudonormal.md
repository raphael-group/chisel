usage: chisel-pseudonormal.py [-h] -r REFERENCE [-x RUNDIR] [-e THRESHOLD]
                              [-b SIZE] [-c CHROMOSOMES] [-m MINREADS]
                              [--samtools SAMTOOLS] [-j JOBS]
                              [--tmpdir TMPDIR] [-n NORMAL]
                              INPUT

CHISEL command to generate a pseudo-matched normal sample by extracting
diploid cells from a barcoded single-cell BAM file.

positional arguments:
  INPUT                 Barcoded single-cell BAM file

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        Reference genome
  -x RUNDIR, --rundir RUNDIR
                        Running directory (default: current directory)
  -e THRESHOLD, --threshold THRESHOLD
                        Minimum fraction of diploid genome to select diploid
                        cells (default: 0.9)
  -b SIZE, --size SIZE  Bin size, with or without "kb" or "Mb"
  -c CHROMOSOMES, --chromosomes CHROMOSOMES
                        Space-separeted list of chromosomes between apices
                        (default: "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8
                        chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17
                        chr18 chr19 chr20 chr21 chr22")
  -m MINREADS, --minreads MINREADS
                        Minimum number total reads to select cells (default:
                        100000)
  --samtools SAMTOOLS   Path to the directory to "samtools" executable,
                        required in default mode (default: samtools is
                        directly called as it is in user $PATH)
  -j JOBS, --jobs JOBS  Number of parallele jobs to use (default: equal to
                        number of available processors)
  --tmpdir TMPDIR       Temporary directory in running directory (default:
                        _TMP)
  -n NORMAL, --normal NORMAL
                        Name of the generated pseudo matched-normal BAM file
                        (default: pseudonormal.bam)

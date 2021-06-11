```shell
usage: chisel_prep [-h] [-r REFERENCE] [-x RUNDIR] [-o OUTPUT]
                   [--rexpname REXPNAME] [--rexpread REXPREAD]
                   [--noduplicates] [--keeptmpdir]
                   [--barcodelength BARCODELENGTH] [--bcftools BCFTOOLS]
                   [--samtools SAMTOOLS] [--bwa BWA] [-j JOBS] [--seed SEED]
                   INPUT [INPUT ...]

CHISEL command to create a barcoded BAM file from single-cell FASTQs (or gz-
compressed FASTQs), single-cell BAMs, or a `RG:Z:`-barcoded BAM files without
`CB:Z:` tags. When single-cell FASTQs or BAMs are provided a CELL name is
assigned to each file (through either filename or table) and the same cell
barcode will be assigned to all corresponding reads, but a different RG tag as
they are considered as different repetitions of sequencing of the same cell.
Specifically, when a table of inputs is not provied, for FASTQs each CELL name
is extracted from the filename through the provided regular expression
(default matches Illumina standard format), for BAMs basename is used as CELL
name. When single-cell FASTQs are provided a READ value is also assigned to
each file (through either filename or table) and files with the same filename
when removing READ values are considered as pairs of sequencing read mates.
Input files, CELL names, and possible READ values can be provided through a
table of inputs.

positional arguments:
  INPUT                 Input FASTQs, BAMs, or TSV file with different
                        behaviors: .........................................
                        (1) FASTQs -- specified in a directory DIR as
                        `DIR/*.fastq` or `DIR/*.fastq.gz` -- will be barcoded
                        and aligned with (optionally) marked duplicates into a
                        barcoded BAM file; .................................
                        (2) BAMs -- specified in a directory DIR as
                        `DIR/*.bam` -- will be barcoded and aligned with
                        (optionally) marked duplicates into a barcoded BAM
                        file; ..............................................
                        (3) a single BAM file with unique cells names in the
                        field `RG:Z:` will be converted into a barcoded BAM
                        file with the additional `CB:Z:` tag; ..............
                        (4) a tab-separated table of inputs (TSV with optional
                        header starting with `#`) with two columns: the first
                        column is an input file (FASTQ or BAM) and the second
                        column is the corresponding cell name. When FASTQs are
                        provided, a third column can be optionally specified
                        to indicate the read name in paired-end sequencing,
                        e.g., indicating either R1 or R2 for the first or
                        second mate of paired-end reads, respectively. If a
                        third column is not present, FASTQs are assumed to be
                        from single-end sequencing.

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        Reference genome, which is mandatory in FASTQ mode
                        (default: None)
  -x RUNDIR, --rundir RUNDIR
                        Running directory (default: current directory)
  -o OUTPUT, --output OUTPUT
                        Output name in running directory (default:
                        barcodedcells.bam)
  --rexpname REXPNAME   Regulare expression to extract cell name from input
                        FASTQ filenames (default:
                        `(.*)_S.*_L.*_R[1|2]_001.fastq.*`)
  --rexpread REXPREAD   Regulare expression to extract cell name from input
                        FASTQ filenames (default:
                        `.*_S.*_L.*_(R[1|2])_001.fastq.*`)
  --barcodeonly         Only compute barcodes but do not run aligning pipeline
                        (default: False)
  --noduplicates        Do not perform marking duplicates and recalibration
                        with Picard tools (default: False)
  --keeptmpdir          Do not erase temporary directory (default: False)
  --barcodelength BARCODELENGTH
                        Length of barcodes (default: 12)
  --bcftools BCFTOOLS   Path to the directory to "bcftools" executable
                        (default: in $PATH)
  --samtools SAMTOOLS   Path to the directory to "samtools" executable
                        (default: in $PATH)
  --bwa BWA             Path to the directory to "bwa" executable (default: in
                        $PATH)
  -j JOBS, --jobs JOBS  Number of parallele jobs to use (default: equal to
                        number of available processors)
  --seed SEED           Random seed for replication (default: None)
```

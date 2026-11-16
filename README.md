# chisel_hmm

A fork of [CHISEL](https://github.com/raphael-group/chisel) ([Zaccaria and Raphael, 2021](https://doi.org/10.1038/s41587-020-0661-6)), which infers allele- and haplotype-specific copy numbers of single cells from low-coverage scDNA-seq. The fork adds variable-length bins, HMM clustering of bins, a new cloner and a per-cell HMM.

## Install

Python 2.7 with CHISEL's dependencies (e.g. from CHISEL's bioconda package), plus samtools and bcftools:

```shell
git clone -b hmm https://github.com/kevinyu2/chisel_hmm.git && cd chisel_hmm
pip install -e .
```

## Required data

1. **Barcoded BAM** of the single cells, sorted and indexed; every read carries its cell barcode (`CB:Z:` by default). `chisel_prep` makes one from per-cell FASTQs or BAMs.
2. **Reference genome** the reads were aligned to, with its `.fai` and `.dict`.
3. **Phased germline SNPs**: a VCF of heterozygous SNPs phased `0|1` / `1|0` (e.g. called from the normal and phased with Eagle2 - use UGP (https://github.com/raphael-group/Universal-Genotyping-Pipeline)).
4. **Matched-normal BAM** (`chisel_nonormal` runs without one).

## Algorithm

- **CHISEL**: computes the RDR and BAF of every cell in every bin, clusters bins across all cells, calls allele-specific copy numbers per cluster and cell (choosing the ploidy), phases them into haplotype-specific copy numbers, and (with `--cloner`) groups cells into clones.
- **Variable-length bins**: 1 kb windows are merged until every bin has enough SNP-covering reads over all cells, so bins carry similar BAF information.
- **HMM bin clustering** (`--hmm`): a HATCHet2-style Gaussian HMM along the genome replaces k-medians; its variances come from each bin's read counts, so short noisy bins do not form their own clusters.
- **Cloner** (`--cloner`): cells are split top-down wherever a principal component beats a permutation noise floor; the clones are then pooled for the HMM and the calls.
- **EXPERIMENTAL Per-cell HMM** (`chisel_cellhmm`): every cell gets its own state in every bin, the cells sharing K (RDR, BAF) states started from copy-number grid points.

## Important options

| option | what it does |
|---|---|
| `--bins FILE` | use these bins (CHR START END), e.g. UGP's |
| `--varbins --minsnpreadspercell N` | compute variable-length bins with N SNP-covering reads per cell per bin (we recommend just running UGP and submitting the --bins file) |
| `--hmm --hmmk K` | HMM bin clustering with K clusters |
| `--hmm --hmmmink A --hmmmaxk B` | HMM, K chosen in A..B |
| `--cloner` | group cells into clones (default: off, every cell is called on its own reads) |
| `--minclone N` | smallest clone (default 30) |
| `--clusters FILE` | (`chisel_calling`) take bin clusters from a file instead of computing them |
| `-K K` | (`chisel_cellhmm`) number of states |

## Examples

```shell
# full pipeline on given variable-length bins, HMM with K chosen in 5-12
chisel -t cells.bam -n normal.bam -r hg38.fa -l phased.vcf -x run/ --bins bins.tsv --hmm --hmmmink 5 --hmmmaxk 12

# re-call an existing run with 40 HMM clusters, and with the cloner
chisel_calling run/combo/combo.tsv -x run_hmm40/ --hmm --hmmk 40
chisel_calling run/combo/combo.tsv -x run_cloner/ --hmm --hmmmink 24 --hmmmaxk 48 --cloner

# per-cell HMM with 6 states
chisel_cellhmm run/combo/combo.tsv -K 6 -o calls.tsv
```

Outputs are in the run folder: `combo/combo.tsv` (RDR and BAF per cell and bin), `calls/calls.tsv` (copy numbers), `clones/mapping.tsv` (with `--cloner`) and `plots/`.

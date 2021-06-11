# Command `chisel.py`

The CHISEL command `chisel.py` runs the entire CHISEL pipeline starting from the required inputs (e.g. BAM files).
During the executiong, the command creates six folders which contain the temporary and final results produced by the 5 distinct steps of CHISEL.

## Estimating RDRs

This step aims to estimate the RDR for every genomic bin in each cell.
Moreover, it selects the barcodes that correspond to cells using a specified threshold on the minimum number of reads.
This step creates a folder `rdr` with three files:

1. `total.tsv`: a TSV dataframe containing the number of sequencing reads observed for every selected cell. More specifically, the fields are:
   1. `CELL`: the name of a cell or the name `normal` indicating the matched-normal sample;
   2. `TOTAL`: the total number of sequencing reads observed for the cell.
2. `rdr.tsv`: a TSV dataframe containg the estimated RDRs with the following fields:
   1. `CHROMOSOME`: the name of a chromosome;
   2. `START`: the starting coordinate of a genomic bin;
   3. `END`: the ending coordinate of the genomic bin;
   4. `CELL`: the name of a cell;
   5. `NORMAL`: the number of sequencing reads from the matched-normal sample for the bin;
   5. `COUNT`: the number of sequencing reads from the cell `CELL` in the bin;
   6. `RDR`: the estimated RDR.
3. `log`: a logging file of the execution of this step (optional).

## Estimating BAF

This step aims to estimate the BAF for phased germline heterozygous SNPs in the selected cells.
This step creates a folder `baf` with two files:

1. `baf.tsv`: a TSV dataframe with the following fields:
   1. `CHROMOSOME`: the name of a chromosome;
   2. `POS`: a genomic position in the chromosome `CHROMOSOME` for a germline heterozygous SNP;
   3. `CELL`: the name of a cell;
   4. `A-COUNT`: the number of observed sequencing reads from the haplotype A of the SNP;
   4. `B-COUNT`: the number of observed sequencing reads from the haplotype B of the SNP.
2. `log`: a logging file of the execution of this step (optional).

## Combining RDRs and BAFs

This step aims to combine the RDRs and BAFs for the selected bins in the selected cells.
This step creates a folder `combo` with two files:

1. `combo.tsv`: a TSV dataframe with the following fields:
   1. `CHROMOSOME`: the name of a chromosome;
   2. `START`: the starting coordinate of a genomic bin;
   3. `END`: the ending coordinate of the genomic bin;
   4. `CELL`: the name of a cell;
   5. `NORMAL`: the number of sequencing reads from the matched-normal sample for the bin;
   6. `COUNT`: the number of sequencing reads from the cell `CELL` in the bin;
   7. `RDR`: the estimated RDR for the bin in the cell `CELL`;
   8. `A-COUNT`: the number of observed sequencing reads from the haplotype A of the SNP;
   9. `B-COUNT`: the number of observed sequencing reads from the haplotype B of the SNP;
   10. `BAF`: the B-allele frequency estimated for the bin in the cell `CELL`.
2. `log`: a logging file of the execution of this step (optional).

## Calling

This step aims to infer the ploidy of each cell and, after global clustering of RDRs and BAFs, to infer the allele- and haplotype-specific copy numbers for every bin in every cell.
This step creates a folder `calls` with two files

1. `calls.tsv`: a TSV dataframe with the following fields:
   1. `CHROMOSOME`: the name of a chromosome;
   2. `START`: the starting coordinate of a genomic bin;
   3. `END`: the ending coordinate of the genomic bin;
   4. `CELL`: the name of a cell;
   5. `NORMAL`: the number of sequencing reads from the matched-normal sample for the bin;
   6. `COUNT`: the number of sequencing reads from the cell `CELL` in the bin;
   7. `RDR`: the estimated RDR for the bin in the cell `CELL`;
   8. `A-COUNT`: the number of observed sequencing reads from the haplotype A of the SNP;
   9. `B-COUNT`: the number of observed sequencing reads from the haplotype B of the SNP;
   10. `BAF`: the B-allele frequency estimated for the bin in the cell `CELL`;
   11. `ALLELECN`: dash-separated ordered pair of the inferred haplotype-specific copy numbers for the bin in the cell `CELL`.
2. `log`: a logging file of the execution of this step (optional).

## Cloning

This steps aims to infer the clones by clustering cells based on the inferred haplotype-specific copy numbers and selecting the clusters that correspond to actual clones.
This step creates a folder `clones` with two files:

1. `mapping.tsv`: a TSV dataframe with the following fields:
   1. `CELL`: the name of a selected cell;
   2. `CLUSTER`: the cluster where the cell `CELL` has been assigned;
   3. `CLONE`: the clone of the cell `CELL`, however it corresponds to `None` if the cells is classified as noisy.
2. `log`: a logging file of the execution of this step (optional).

Moreover, this step introduces a new field (right-most field) in the file `calls.tsv` which is `CORRECTED_CNS` and corresponds to the final haplotype-specific copy numbers estimated after consensus of cells in the same clone.

## Plotting

This step generate several useful plots about the results, which are fully described [here](chisel-plotting.md).

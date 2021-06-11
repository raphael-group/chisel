# Command `chisel-cloning.py`

The CHISEL command `chisel-cloning.py` runs the CHISEL pipeline starting from the already estimated allele- and haplotype-specific copy numbers.
To do this, this command requires to have the folder `calls` with the files and formats described [here](chisel.md).
This command is particularly useful if the user would like to re-run the CHISEL's inference of tumor clones to adapt to datasets with particularly high levels of noise and variance.
Examples of usage of this command for QC is described [here](../guides/clones.md).

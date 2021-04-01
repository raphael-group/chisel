# Identification of clones

CHISEL infers clones by clustering cells based on the inferred haplotype-specific copy numbers and selecting the clusters that correspond to actual clones.
This selection is indeed required because the data are noisy and minor differences between cells may indicate errors as well as small clusters may indicate noisy cells with bad sequencing.
This identification is controlled by two parameters:
- `f`: the maximum fraction of the genome with different haplotype-specific copy numbers for cells in the same clone (default: 0.06);
- `s`: the minimum number of cells in a clone (default: 14).

The values of these two parameters have been calibrated for the expected number of cells and sequencing coverage of 10X Genomics datasets.
However, when analyzing datasets with different number of cells, different sequencing coverage, or particularly noisy datasets, the default values of these parameters may not be appropriate.
Therefore, when the user observes a outlying high number of noisy cells or too few inferred clones (even 0), it is important to vary these values to explore different solutions.

Given the inferred clones with the previous parameters, there is one additional parameter that can be used to adjust the classification of noisy cells: `-r`, which controls the refinement of the identified clones and allows the user to include more "noisy" cells into the identified clones. Specifically, every cell that has a fraction of the genome with different haplotype-specific copy numbers lower than then value of `r` will be included into the clones. Therefore, the user can user increasingly higher values to force the inclusion of more noisy cells into the inferred clones, for example `-r 0.2`, `-r 0.3`, `-r 0.4`, etc. Note that `-r 1` will force every cell to be assigned to a clone.

These tasks can be performed very efficiently and easily by using the [CHISEL command `chisel-cloning.py`](../doc/chisel-cloning.md), which allows the user to only re-execute the inference of clones and the generation of plots very efficiently from the already inferred haploty-specific copy numbers.
As such, using this command, the user can attempt to use different combinations of the parameters, varying the maximum difference `f` (e.g. `-f 0.1`, `-f 0.12`, `-f 0.15`, ...) and the minimum number `-s` of cells to select the clones (either increasing like `-s 20`, `-s 30`, ... or decreasing like `-s 3`, `-s 2`, according to the total number of cells).
More details on adjusting and selecting reasonable values of these parameters are available in the [CHISEL's manuscript](https://doi.org/10.1101/837195).

# Command `chisel-calling.py`

The CHISEL command `chisel-calling.py` runs the CHISEL pipeline starting from the already estimated RDRs and BAFs.
To do this, this command requires to have the folder `combo` with the files and formats described [here](chisel.md).
This command is particularly useful if the user would like to re-run the CHISEL pipeline without the extensive re-computation of RDRs and BAFs but using different values of some of the main parameters, including:

1. `-A`: varying sensitiviy of the model selection criterion for cell ploidy: in case of particularly noisy datasets or with particularly high variance, the estimation of cell ploidy may be more challenging and it may needed to increase the sensitivity of the selection (e.g. 2, 3, 4, ...);
2. `-K`: varying the maximum number of clusters allowed in the global clustering of RDRs and BAFs: choosing values lower than the default (i.e. 100) generally allows to reduce presence of noisy CNAs at the cost of lower resolution;
3. `-P`: varying the maximum value allowed for cell ploidy, since the default is 4 which generally corresponds to at most one WGD.


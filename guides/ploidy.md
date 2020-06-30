# Ploidy selection

CHISEL infers the ploidy of each cell from the estimated RDRs and BAFs by using a model-selection criterion which has been calibrated for observing a sufficient number of reads on average across all bins and cells (further details are reported in the [CHISEL's manuscript](https://doi.org/10.1101/837195)).
In the case of particularly noisy cells or datasets with a particularly low sequencing coverage or high variance, the inference can thus be more challenging.
A sign which generally indicates potential issues in the inference of cell ploidies is the observation of a substatial number of cells with different ploidies (and thus with completely different copy numbers).
However, CHISEL provides a parameter to adjust the sensitivity of the model-selection criterion for dealing with these cases.

For QC purposes, the reccommendation is to analyze the allele-specific copy numbers inferred by CHISEL, for example using the corresponding [plots](../doc/chisel-plotting.md).
If a substantial number of cells with different ploidies has been inferred, the reccommendation is to analyze how the results change by re-running the inference of copy numbers varying the sensitivity of the model-selection criterion.
Specifically, the CHISEL [command `chisel-calling.py`](../doc/chisel-calling.py) can be used to do this very efficiently by varying the sensitivity with the argument `-A` to re-run the CHISEL inference without the need of re-estimating RDRs and BAFs (which generally is the most time-consuming step).
The use can thus analyse the inferred tumor ploidies by increase the sensitivity with values of `-A 2`, `-A 3`, `-A 4`...
The inference of different ploidies is well supported by the data if the results do not substantially change when increasing the sensitivity, otherwise the results obtained with higher sensitivity are more likely.

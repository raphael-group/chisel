# Clustering

CHISEL globally clusters the estimated RDRs and BAFs by using a k-means algorithm and model-selection criterion based on the elbow method to select the best number of clusters (further details are reported in the [CHISEL's manuscript](https://doi.org/10.1101/837195)).
In order to do this, CHISEL fixes the maximum number of clusters (default value is 100).
However, the value can be too high when analyzing very noisy datasets, since the highe levels of noise in the data can be misinterpreted as actual signal and can lead to overfitting.
Therefore, the user can assess the levels of variance and noise by using the [BAF and RDR plots](chisel-plotting.md) and, when analyzing datasets with very high levels of variance, it may lower the maximum number of clusters to avoide overfitting.
A possible signal that may indicate overiffint is for example the observation of many outlying and noisy CNAs from the inferred allele-specific copy numbers, i.e. observing many cells with isolated and small CNAs.
For such QC purposes, the user can use the CHISEL [command `chisel-calling.py`](chisel-calling.py) to vary the maximum number of clusters with the argument `-K` to re-run the CHISEL inference without the need of re-estimating RDRs and BAFs (which generally is the most time-consuming step).
The user can thus quantify the presence of noisy CNAs by varying the value of this parameter, for example `-K 80`, `-K 60`, `-K 40`...

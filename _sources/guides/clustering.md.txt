# Clustering

CHISEL globally clusters the estimated RDRs and BAFs by using a k-means algorithm and model-selection criterion based on the elbow method to select the best number of clusters (further details are reported in the [CHISEL's manuscript](https://doi.org/10.1101/837195)).
In order to do this, CHISEL fixes the maximum number of clusters (default value is 100).
However, the value can be too high when analyzing very noisy datasets since the high levels of noise in the data can be misinterpreted and may lead to overfitting.
Therefore, the user can assess the levels of variance and noise by using the [BAF and RDR plots](../chisel-plotting.html): in particular, high levels of noise can be immediately noted when a clear clustering structure is missing from such plots. 
As such, when analyzing datasets with very high levels of variance, the user can lower the maximum number of clusters to avoid overfitting.
Another possible signal that may indicate overiffing is for example the inference of many outlying and noisy CNAs, i.e. observing many cells with isolated and small CNAs.
For such QC purposes, the user can use the CHISEL [command `chisel-calling.py`](../chisel-calling.html) to vary the maximum number of clusters with the argument `-K` to re-run the CHISEL inference without the need of re-estimating RDRs and BAFs (which generally is the most time-consuming step).
The user can thus quantify the presence of noisy CNAs when varying the value of this parameter, for example `-K 80`, `-K 60`, `-K 40`...

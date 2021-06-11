# Command `chisel-plotting.py`

The CHISEL command `chisel-plotting.py` generates several useful plots, which can be used to inspect the inferred results or for quality control.
More specifically, this command generates 15 plots.

## Main plots

### Allele-specific copy numbers

This plot (`allelecn.png`) depicts the allele-specific copy numbers inferred by CHISEL for every cell (rows) along the entire genome (y-axis), with cells/rows colored according to the inferred clone and columns/genomic bins colored according to the chromosome.
The colors of the heatmap represent the difference pairs of allele-specific copy numbers and a full description of the color map used is available in the [CHISEL manuscript](https://doi.org/10.1101/837195).

### Corrected allele-specific copy numbers

This plot (`allelecn-corrected.png`) depicts the allele-specific copy numbers inferred by CHISEL and corrected using the inferred clones for every cell (rows) along the entire genome (y-axis), with cells/rows colored according to the inferred clone and columns/genomic bins colored according to the chromosome.
The colors of the heatmap represent the difference pairs of allele-specific copy numbers and a full description of the color map used is available in the [CHISEL manuscript](https://doi.org/10.1101/837195).

### Haplotype-specific copy numbers

This plot (`haplotypecn.png`) depicts the haplotype-specific copy numbers inferred by CHISEL for every cell (rows) along the entire genome (y-axis), with cells/rows colored according to the inferred clone and columns/genomic bins colored according to the chromosome.
The colors of the heatmap represent the haplotype of the allele with fewer copies, such that green and magenta represent haplotype A and B, respectively.
Note that balanced regions (allele with the same number of copies) are represented in white.
Further descriptions of the color map used are available in the [CHISEL manuscript](https://doi.org/10.1101/837195).

### Corrected haplotype-specific copy numbers

This plot (`haplotypecn-corrected.png`) depicts the haplotype-specific copy numbers inferred by CHISEL and corrected using the inferred clones for every cell (rows) along the entire genome (y-axis), with cells/rows colored according to the inferred clone and columns/genomic bins colored according to the chromosome.
The colors of the heatmap represent the haplotype of the allele with fewer copies, such that green and magenta represent haplotype A and B, respectively.
Note that balanced regions (allele with the same number of copies) are represented in white.
Further descriptions of the color map used are available in the [CHISEL manuscript](https://doi.org/10.1101/837195).

## Useful additional plots

### BAF and RDR plots

This plot (`rbplot_mirrored.png`) shows the global clusters of RDRs and BAFs inferred for a random sample of a certain number of cells (by default 20 cells).
Each plot corresponds to a different cell, with each plot depicting the bins (each point) which are represented by the corresponding values of `|0.5 - mirrored BAF|` (x-axis) and RDR (y-axis), and are colored according to the corresponding cluster; note that colors are consistent across cells.

### Clustered RDR

This plot (`crdr.png`) shows the estimated RDR and their cluster for a random sample of a certain number of cells (by default 20 cells).
Each plot corresponds to a different cell, with each plot depicting the bins (each point) which are represented by the corresponding values of RDR (y-axis) along the entire genome (x-axis) and are colored according to the corresponding cluster; note that colors are consistent across cells.

### Clustered mirrored BAF

This plot (`cbaf.png`) shows the estimated BAF and their cluster for a random sample of a certain number of cells (by default 20 cells).
Each plot corresponds to a different cell, with each plot depicting the bins (each point) which are represented by the corresponding values of `|0.5 - mirrored BAF|` (y-axis) along the entire genome (x-axis) and are colored according to the corresponding cluster; note that colors are consistent across cells.

### Total copy numbers

This plot (`totalcn.png`) is an heatmap that shows the total copy numbers inferred by CHISEL for every cell (rows) along the entire genome (y-axis), with cells/rows colored according to the inferred clone and columns/genomic bins colored according to the chromosome.
Each point of the heatmap this represents the total copy numbers inferred by CHISEL, such that grey represents 2 copies, blue colors represent <2 copies with darker colors corresponding to smaller values, and red colors represent >2 copies with darker colors corresponding to higher values.

### Corrected total copy numbers

This plot (`totalcn-corrected.png`) is an heatmap that shows the total copy numbers inferred by CHISEL and corrected using the inferred clones for every cell (rows) along the entire genome (y-axis), with cells/rows colored according to the inferred clone and columns/genomic bins colored according to the chromosome.
Each point of the heatmap this represents the total copy numbers inferred by CHISEL, such that grey represents 2 copies, blue colors represent <2 copies with darker colors corresponding to smaller values, and red colors represent >2 copies with darker colors corresponding to higher values.

### LOH

This plot (`loh.png`) is an heatmap that shows the LOH inferred by CHISEL for every cell (rows) along the entire genome (y-axis), with cells/rows colored according to the inferred clone and columns/genomic bins colored according to the chromosome.
Each point of the heatmap is colored according to the absence (white) or presence (black) of a LOH in the corresponding cell and bin.

### Corrected LOH

This plot (`loh-corrected.png`) is an heatmap that shows the LOH inferred by CHISEL and corrected using the inferred clones for every cell (rows) along the entire genome (y-axis), with cells/rows colored according to the inferred clone and columns/genomic bins colored according to the chromosome.
Each point of the heatmap is colored according to the absence (white) or presence (black) of a LOH in the corresponding cell and bin.

### A-specific copy numbers

This plot (`Aspecificcn.png`) is an heatmap that shows the copy numbers inferred by CHISEL for haplotype A for every cell (rows) along the entire genome (y-axis), with cells/rows colored according to the inferred clone and columns/genomic bins colored according to the chromosome.
The colors of the heatmap are the same for the total copy-numbers (see above).

### Corrected A-specific copy numbers

This plot (`Aspecificcn-corrected.png`) is an heatmap that shows the copy numbers inferred by CHISEL for haplotype A and corrected using the inferred clones for every cell (rows) along the entire genome (y-axis), with cells/rows colored according to the inferred clone and columns/genomic bins colored according to the chromosome.
The colors of the heatmap are the same for the total copy-numbers (see above).

### B-specific copy numbers

This plot (`Bspecificcn.png`) is an heatmap that shows the copy numbers inferred by CHISEL for haplotype B for every cell (rows) along the entire genome (y-axis), with cells/rows colored according to the inferred clone and columns/genomic bins colored according to the chromosome.
The colors of the heatmap are the same for the total copy-numbers (see above).

### Corrected B-specific copy numbers

This plot (`Bspecificcn-corrected.png`) is an heatmap that shows the copy numbers inferred by CHISEL for haplotype B and corrected using the inferred clones for every cell (rows) along the entire genome (y-axis), with cells/rows colored according to the inferred clone and columns/genomic bins colored according to the chromosome.
The colors of the heatmap are the same for the total copy-numbers (see above).



# Command `chisel-pseudonormal.py`

The CHISEL command `chisel-pseudonormal.py` implements the method integrated in CHISEL for generating a pseudo matched-normal sample by extracting diploid cells from a barcoded BAM file.
The command simply required as input a barcoded BAM file and the corresponding reference genome; detailed descriptions of the required input are available [here](../man/chisel-pseudonormal.md).
After the execution, the command generates a new BAM file by only merging sequencing reads from diploid cells; thus the resulting BAM file can be used as a paseudo matched-normal sample to run the entire [CHISEL's pipeline](chisel.md).
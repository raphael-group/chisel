usage: chisel-plotting.py [-h] [-m CLONEMAP] [-f FIGFORMAT] [-s SAMPLE]
                          [--excludenoisy] [--gridsize GRIDSIZE]
                          [--plotsize PLOTSIZE] [--clussize CLUSSIZE]
                          [--xmax XMAX] [--xmin XMIN] [--ymax YMAX]
                          [--ymin YMIN]
                          [INPUT]

CHISEL command to re-create the plots.

positional arguments:
  INPUT                 Input file with inferred copy numbers (default:
                        calls/calls.tsv)

optional arguments:
  -h, --help            show this help message and exit
  -m CLONEMAP, --clonemap CLONEMAP
                        Clone map (default: not used, the cells will be
                        clustered for plotting purposes)
  -f FIGFORMAT, --figformat FIGFORMAT
                        Format of output figures (default: png, the only other
                        option is pdf)
  -s SAMPLE, --sample SAMPLE
                        Number of cells to sample (default: 20)
  --excludenoisy        Exclude noisy cells from plots (default: False)
  --gridsize GRIDSIZE   Grid dimenstions specified as comma-separated numbers
                        (default: 12,6)
  --plotsize PLOTSIZE   Plot dimenstions for RDR-BAF plots, specified as
                        comma-separated numbers (default: 5,1.5)
  --clussize CLUSSIZE   Grid dimenstions for clustered plots, specified as
                        comma-separated numbers (default: 5,3)
  --xmax XMAX           Maximum x-axis value (default: None)
  --xmin XMIN           Minimum x-axis value (default: None)
  --ymax YMAX           Maximum x-axis value (default: None)
  --ymin YMIN           Minimum x-axis value (default: None)

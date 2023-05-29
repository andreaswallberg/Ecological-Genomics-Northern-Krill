# Scripts for analysis of population structure

## ld_prune.py

Script to prune a vcf file from variants that are in linkage. The
script was used as a pre-processing step for principal components
analysis and run with the following parameters:

    python ld_prune.py --subsample-fraction 0.05 --threshold 0.1 --as-bed input.vcf.gz

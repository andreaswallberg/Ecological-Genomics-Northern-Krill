# Scripts for processing DNA-methylation output from f5c

## methylation_frequencies_subdivide_repeats_bed.pl
This script parses DNA methylation frequencies provided in this tabular format:

    # chromosome      start   end     num_motifs_in_group     called_sites    called_sites_methylated methylated_frequency    group_sequence
    # ctg1    77      93      3       15      0       0.000   GCTTCCGCTTGTGTCGCTCCACGTGCC

In cases when the CpG locus / group_sequence contains more than one CpG site, the scripts splits them into individual sites with the same frequencies.

The scripts requires a BED file with repeat coordinates. It partitions the CpG sites and methylation frequencies according to whether they occur in repeats or not. However, the BED file can empty and all methylation data will then be written out as non-repeated.

The methylation call file can be provided in as regular text or be gzipped.

### Usage:

    ./methylation_frequencies_subdivide_repeats_bed.pl repeats.bed methylation_calls.tsv.gz

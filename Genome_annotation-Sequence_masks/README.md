### Scripts for generating genome masks and associated files

## gtf2gene_regions.comprehensive.pl

This script generates two per-base gene region masks for the genome using a table of sequence lengths and a GTF/GFF3 file with gene and CDS coordinates.

The length information can either be provided using a tabular "lengths" file (using "-length"), which assumes a header line is present on the first line and that sequence name and length are provided in the second and third columns, respectively. Alternatively, a BED file can be provided (using "-bed" instead), which is assumed to have sequence name and length in the first and second columns and to have no header line. More than one such file can be provided and both types can be used at the same time.

It prints one mask with the parsed gene regions and a second mask where in accessible sites have been set to "0". The character states are:

    0 => inaccessible
    1 => intergenic
    2 => intron
    3 => three_prime_utr
    4 => exon
    5 => five_prime_utr
    6 => cds

To generate the second file it needs a per-base accessibility mask that specifys whether the site is accessible or not, assuming the following character states:

    0 => inaccessible
    1 => accessible

### Usage example:

    gtf2gene_regions.comprehensive.pl \
        --length genome_sequence_lengths.csv \
        --gtf gene_annotation.gff3 \
        --coverage genome_mask_accessible.fasta \
        --output genome_mask_region

The script generates intergenic strings of "1"s for every sequence and then populate them with gene regions according to the annotations.

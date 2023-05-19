# Scripts for processing DNA-methylation output from f5c

## methylation_frequencies_subdivide_repeats_bed.pl

This script parses DNA methylation frequencies provided in this tabular format:

    # chromosome      start   end     num_motifs_in_group     called_sites    called_sites_methylated methylated_frequency    group_sequence
    # ctg1    77      93      3       15      0       0.000   GCTTCCGCTTGTGTCGCTCCACGTGCC

In cases when the CpG locus / group_sequence contains more than one CpG site, the scripts splits them into individual sites with the same frequencies.

The scripts requires a BED file with repeat coordinates. It partitions the CpG sites and methylation frequencies according to whether they occur in repeats or not. However, the BED file can empty and all methylation data will then be written out as non-repeated.

The methylation call file can be provided in as regular text or be gzipped.

### Usage example:

    ./methylation_frequencies_subdivide_repeats_bed.pl repeats.bed methylation_calls.tsv.gz


## methylation_frequencies2gene_region_repeats_coverage.genome_mask.pl

This scripts also parses DNA methylation provided in this tabular format and splits CpG loci into individual sites (as above).

It reads up to three genome masks in FASTA format:

1. A per-base gene region mask (required), assuming the following character states:
        
        1 => intergenic
        2 => intron
        3 => 3_UTR
        4 => exon
        5 => 5_UTR
        6 => cds

2. A per-base repeat mask (optional), assuming the following character states:

        0 => non-repeated
        1 => repeated

3. A per-base coverage file specifying whether the site is accessible or not (optional), assuming the following character states:
    
        0 => inaccessible
        1 => accessible
    
It can also read a tabular text file ("-hetero") that specifies particular sites to exclude from the analysis (optional), for example heterozygous CpG locations. It assumes the format (1-based positions):
    
    name_of_sequence    position

### Usage example:

    ./methylation_frequencies2gene_region_repeats_coverage.genome_mask.pl \
	-dna 10 \
	-regions genome_mask_region.fasta \
        -repeats genome_mask_repeats.fasta \
        -coverage genome_mask_accessible.fasta \
        -hetero heterozygous_sites.tsv \
	-output methylation_calls.tsv.gz.out \
	-repeat_location 1 \
	-methylation methylation_calls.tsv.gz
        
* The "-methylation" option can specify more than one file.
* The "-repeat_location" option means:
    0 = skip CpGs in repeats
    1 = only consider CpGs in repeats
    2 = do not distinguish between CpGs in repeats or non-repeat sequences (i.e. look at all CpGs)


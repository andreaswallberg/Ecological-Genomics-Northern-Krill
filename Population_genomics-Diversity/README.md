# Scripts for computing levels and patterns of genetic diversity

These scripts process VCFs to estimate levels and patterns of genetica variation (e.g. Watterson's theta, Pi and Tajima's D).

## vcf2watterson_pi_full_seqs_fast_coverage_gene_regions.pl

This script estimates Watterson's theta (population mutation rate), Pi (nucleotide diversity) and Tajima's D.

The code for Tajima's D originates from BioPerl and was originally written by Jason Stajich

https://metacpan.org/release/CJFIELDS/BioPerl-1.6.924/source/Bio/PopGen/Statistics.pm

The implementarion here differs by reusing some already calculated variables.

### Usage example:

    time vcf2watterson_pi_full_seqs_fast_coverage_gene_regions.pl \
        --vcf phased_imputed_annotated_snps.vcf.gz \ # Input VCF
        --output phased_imputed_annotated_snps.vcf.gz.diversity \ # Output basename
        --group \ # Populations/groups to estimate variation statistics for
            at=samples.at.csv \
            me=samples.me.csv \
        --window 1000 10000 100000 \ # Size(s) of non-overlapping windows to estimate variation across
        --coverage genome_mask_accessible_sites.fasta \ # A genome mask (this can be a binary-state accessibility mask or one encoding different genomic regions such as CDS, intron using multiple states)
        --seqs subset.bed \ # A bed file specifying which sequences to estimate FST values from (OPTIONAL, if not specified it will scan all sequences in the allele count file)
        --region \ # Which genomic regions to look at (refers back to the states encoded in the provided genome mask)
            1=intergenic \
            2=intron \
            3=three_prime_utr \
            4=exon \
            5=five_prime_utr \
            6=cds \
        --verbose

Additional options and other aspects of usage:
- The "--vcf" argument can take more than one VCF file.
- Samples should not be included in more than one group
- The "--window" argument can take any number of window sizes (at the cost of performance)
- If a "-genes" argument is provided together with a GFF/GTF file, the script will compute the average levels of diversity at different distances upstream/downstream of genes.
- The "groups" argument specifies the label of the group/population and points to a simple csv file with the names of the samples for that group (one sample name per line).

## xp2windows.pl

This script compiles window-based cross-population EHH summary statistics for normalized XP-nSL or XP-EHH estimates produced by selscan.
Selscan is here: https://github.com/szpiech/selscan

It produces tabular output with the following format:

- CHROM = name of sequence
- START = start position of the window
- STOP = stop position of the window
- N = number of SNPs in the window
- N_CRIT = number of critical SNPs (values higher than 2 or lower than 2)
- PROP_CRIT = proportion of critical SNPs
- MIN = SNP with lowest value
- MAX = SNP with highest value
- MEAN = Mean value of all SNPs in the window
 
### Usage example:

    xp2windows.pl \
        1000 \ # Window size
        xp.out \ # Basename of output
        `ls -v *.xpnsl.out.norm` # list of normalized files to look at

## consolidate_window_stats.pl

This is a helper script to consolidate multiple sources of window-based estimates, such as FST, diversity and cross-population XP-EHH/nSL statistics computed with other scripts in this collection. The script aggregates the data for every window and produces a new tabular output file. It assumes that the window-resolution is the same across all datasets. It essentially builds a flat-file tabular database of estimates.

### Usage example:


    ./consolidate_window_stats.pl \
        -output consolidated_windows.1000bp.tsv \ # The tabular output file
        -chromosome_list subset.bed \ # A bed file specifying which sequences to print info for
        -fst_windows fst_windows.1000bp.csv \ # Window-based FST values
        -xp_windows xp.out.1000bp.csv \ # cross-population EHH estimates
        -diversity_windows phased_imputed_annotated_snps.vcf.gz.diversity.wattersons_theta_pi_tajimas_D.window_1000.any.csv \ # Diversity estimates (overall for whole dataset)
        -diversity_comparative_windows phased_imputed_annotated_snps.vcf.gz.diversity.wattersons_theta_pi_tajimas_D.window_1000.any.csv \ # Diversity estimates (subdivided by populations)
        -diversity_region_windows ALL.regions.consolidated_windows.1000bp.csv.filtrered.csv \ # Diversity estimates (subdivided by genomics regions, such as CDS, intron etc.)
        
The downstream **get_fields.pl** script can be used to extract particular fields:

    ./get_fields.pl \
        -table consolidated_windows.1000bp.tsv \ # The consolidated tabular file
        -fields 0 1 4 10 \ # The fields/columns to extract
        -inverse_fields 10 \ # Multiplies values in this column by -1
        -rename_headers \ # Renames output headers for convenience
            4=FST \
        -tag "fst_diff.at_vs_me" Adds a tag to the output file name (which is based on the input file name)
        
 ## stat_over_fst2mean_ci.pl
 
This script computes basic stats from window-based data to show associations between two parameters, such as FST and diversity or gene content. It partitions response-variable data according to a parameter scaled between 0 to 1, such as FST, using a user-specified number of categories. For every partition of data it computes the mean and 95% confidence intervals using non-parametric bootstrapping.

### Usage example:

    ./stat_over_fst2mean_ci.pl \
        "abs" \
        5 \
        2 3 \
        fst_and_diversity_windows_1000bp.tsv \
        > fst_and_diversity_windows_1000bp.tsv.out

Here the arguments are, in turn:
- How to treat or filter the data for the response-variable (i.e. "Y-axis" data): any=include any data point; abs=use the absolute value of any data point; neg=use only negative values; pos=use only positive values
- How many categories of the "X-axis" or "explanatory" parameter that observations should be partitioned into.
- Which column contains X-axis data and which contains Y-axis data
- The tabular input file with the data

### Usage example:

    ./stat2distance_from_genes.95.pl \
        -window 1000 \ # Window size
        -genes genes.gff3 \ # Gene coordinates in GFF format
        -stats diversity_windows_1000bp.tsv \ # Window-based estimates
        -fields 3 4 5 \ # Which fields to extract data from (can be more than one as in this example)
        -output diversity_windows_1000bp.tsv.distances_away_from_genes.out

## stat2distance_from_genes.pl
 
This script computes the average and 95% non-parametric bootstrap interval of precalculated window-based features such as levels of genetic variation at increasing distaces away from genes. It reads a set of gene coordinates from a GFF file and window-based estimates from a tabular text file. It assumes the first row of the table contains column headers.

It compiles all data points each given distance away from all genes into genome-wide distance-classes (e.g. 1000, 2000, 3000, ...) and bootstraps the observations for each class to generate confidence intervals.

### Usage example:

    stat2distance_from_genes.95.pl \
        -window 1000 \ # Window size
        -genes genes.gff3 \ # Gene coordinates in GFF format
        -stats diversity_windows_1000bp.tsv \ # Window-based estimates
        -fields 3 4 5 \ # Which fields to extract data from (can be more than one as in this example)
        -output diversity_windows_1000bp.tsv.distances_away_from_genes.out


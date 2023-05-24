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

Aspects of usage:
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

    time xp2windows.pl \
    1000 \ # Window size
    xp.out \ # Basename of output
    `ls -v *.xpnsl.out.norm` # list of normalized files to look at


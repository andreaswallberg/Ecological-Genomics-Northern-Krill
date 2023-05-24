# Scripts for computing levels and patterns of genetic diversity

These scripts process VCFs to estimate levels and patterns of genetica variation (e.g. Watterson's theta, Pi and Tajima's D).

## vcf2watterson_pi_full_seqs_fast_coverage_gene_regions.pl

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
        
  The "groups" argument specifies the label of the group/population and points to a simple csv file with the names of the samples for that group (one sample name per line).

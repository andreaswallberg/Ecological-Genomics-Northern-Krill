# Scripts for inferring demographic history and haplotype ages

These scripts convert datasets into formats suitable to infer historical fluctuations of effective population size and allele/haplotype ages.

## vcf2msmc_coverage.pl

This script reads a VCF file and exports data in a format compatible with MSMC. By providing a list of samples in a simple text file, it can produce datasets with any number of haplotypes. It uses a accessibility mask to keep track on the number of accessible sites observed between SNPs.

The format is described in more detail on page of the original tool:
https://github.com/stschiff/msmc/blob/master/guide.md

# Usage example:

    vcf2msmc_coverage.pl \
            --vcf snps_annotated.vcf.gz \ # A VCF file with SNPs
            --output snps_annotated.vcf.gz.out_for_msmc \ # Basename of the output file
            --group msmc_run_label=samples.csv \ # A "label" pointing to a corresponding list of samples (one sample per line)
            --coverage genome_mask_accessible.fasta \ # A binary-state genome mask (0=inaccessible; 1=accessible) to keep track on the number of accessible sites
            --seqs sequences.gt_500kbp.bed \ # A BED file specifying which sequences to consider
            --verbose \ # Print some extra progress statements
          
 

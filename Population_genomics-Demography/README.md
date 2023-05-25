# Scripts for inferring demographic history and haplotype ages

These scripts convert datasets into formats suitable to infer historical fluctuations of effective population size and allele/haplotype ages.

## vcf2msmc_coverage.pl

# Usage example:

  vcf2msmc_coverage.pl \
          --vcf snps_annotated.vcf.gz \ # A VCF file with SNPs
          --output snps_annotated.vcf.gz.out_for_msmc \ # Basename of the output file
          --group msmc_for_sample=samples.csv \ # A "run id tag" pointing to a list of samples (one sample per line)
          --coverage genome_mask_accessible.fasta \ # A binary-state genome mask (0=inaccessible; 1=accessible) to keep track on the number of observed sites
          --seqs sequences.gt_500kbp.bed \ # A BED file specifying which sequences to consider
          --verbose \ # Print some extra progress statements
          
 

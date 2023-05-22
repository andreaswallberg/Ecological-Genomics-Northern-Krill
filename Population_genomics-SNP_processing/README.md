# Scripts for processing SNPs

These scripts process variants in VCF files produced by FreeBayes (but are likely compatible with other callers such as GATK).

## vcf2filtered_vcf_by_coverage.pl

This script takes a per-base genome mask with binary states (0=inaccessible; 1=accessible) in FASTA format and filters SNPs in VCF format such that only those that fall inside accessible regions are kept.

The VCF input can be compressed with gzip.

It generates three VCF files:

- ".keep.vcf"       -> VCF file with SNPs to keep (detected at accessible sites)
- ".removed.vcf"    -> VCF file with SNPs to that were removed (detected at inaccessble sites)
- ".discarded.vcf"  -> VCF file with "invalid" variants (i.e. involving multiple alleles or multiple nucleotides) 

In addition, a BED file needs to be specified with the set of genome sequences that should be processed. SNPs on other sequences are ignored.

Usage example:

    vcf2filtered_vcf_by_coverage.pl \
      --vcf snps.vcf \
      --coverage genome_mask_accessible_sites.fasta \
      --seqs genome_sequences.bed \
      --verbose

## vcf_biallelic2fasta.pl

This script filters SNPs in VCF format with the aim to keep only those that:

- are biallelic
- meet a minimum and maximum treshhold for total sequencing depth ("--min_depth" and "--max_depth")
- meet a minimum treshold for the proportion of genotyped samples at a SNP ("--min_fill_position")

It outputs a new VCF with variants, as well as variants written in FASTA and GENO (e.g. for EIGENSTRAT) formats.

Usage example:
    
    vcf_biallelic2fasta.pl \
    --input snps.keep.vcf.gz \ # The VCF input
    --output snps.keep.vcf.gz.biallelic.FILTERED \ # Basename of the output files
    --min_fill_sample 0.5 \ # Minimum genotyping rate for a sample to be included in the FASTA output
    --min_fill_position 0.5 \
    --min_depth 94
    
The FASTA output contains one string per sample. An example:

    >bar_1
    TCGAAT
    >bar_10
    TTGGAA

Here, the sample "bar_1" is heterozygous at three SNPs (T/C, G/A, A/T), while bar_10 is homozygous at all three SNPs. Heterozygous genotypes are always sorted with the reference allele first.

The GENO format contains data in SNP per line and one column per sample. The symbol encoding means:

      0 => zero copies of reference allele.
      1 => one copy of reference allele.
      2 => two copies of reference allele.
      9 => missing data.

More information about the GENO format is available here: https://reich.hms.harvard.edu/software/InputFileFormats


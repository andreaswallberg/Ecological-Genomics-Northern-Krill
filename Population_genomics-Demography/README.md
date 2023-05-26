# Scripts for inferring demographic history and haplotype ages

These scripts convert datasets into formats suitable to infer historical fluctuations of effective population size and allele/haplotype ages.

## vcf2psmc_coverage.pl

This script reads a VCF file and exports data in the FASTA-like "psmcfa" format compatible with PSMC. The script overlays the pattern of heterozygous genotypes from a single specimen on a pre-calculated genome accessibility mask (per-window resolution).

The genome-mask should have one symbol per window:
- N = missing data
- T = data present

For every window in which there is at least on heterozygous SNP, the script overlays the symbol "K" on top of the genome-mask. 

The format is described in more detail on the page of the original tool:
https://github.com/lh3/psmc

    ./vcf2psmc_coverage.pl \
        --vcf snps.vcf.gz \ # The VCF file (can be gzipped)
        --windows 100 \ # The window-resolution
        --output snps.vcf.gz.psmcfa \ # The output file
        --sample samplename \ # The name of the sample in the VCF file to use
        --coverage genome_mask_accessible_100bp_windows.fasta \
        --seqs sequences.gt_500kbp.bed \ # A BED file specifying which sequences to consider
        --verbose \ # Print some extra progress statements

The small helper script **coverage_to_windows.pl** can be used to convert a per-base binary-state genome mask (0=inaccessible; 1=accessible) in FASTA format into the necessary window-based format:

    ./coverage_to_windows.pl 100 genome_mask_accessible.fasta > genome_mask_accessible_100bp_windows.fasta

## vcf2msmc_coverage.pl

This script reads a VCF file and exports data in a format compatible with MSMC. By providing a list of samples in a simple text file, it can produce datasets with any number of haplotypes. It uses a accessibility mask to keep track on the number of accessible sites observed between SNPs.

The format is described in more detail on the page of the original tool:
https://github.com/stschiff/msmc/blob/master/guide.md

# Usage example:

    ./vcf2msmc_coverage.pl \
            --vcf snps_annotated.vcf.gz \ # A VCF file with SNPs (can be gzipped)
            --output snps_annotated.vcf.gz.out_for_msmc \ # Basename of the output file
            --group msmc_run_label=samples.csv \ # A "label" pointing to a corresponding list of samples (one sample per line)
            --coverage genome_mask_accessible.fasta \ # A binary-state genome mask (0=inaccessible; 1=accessible) to keep track on the number of accessible sites (per-base resolution)
            --seqs sequences.gt_500kbp.bed \ # A BED file specifying which sequences to consider
            --verbose \ # Print some extra progress statements
          
 

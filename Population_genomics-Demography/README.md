# Scripts for inferring demographic history and haplotype ages

These scripts convert datasets into formats suitable to infer historical fluctuations of effective population size, recombination rates and allele/haplotype ages.

## vcf2psmc_coverage.pl

This script reads a VCF file and exports data in the FASTA-like "psmcfa" format compatible with PSMC. The script overlays the pattern of heterozygous genotypes from a single specimen on a pre-calculated genome accessibility mask (per-window resolution).

The genome-mask should have one symbol per window:
- N = missing data
- T = data present

For every window in which there is at least on heterozygous SNP, the script overlays the symbol "K" on top of the genome-mask. 

The format is described in more detail on the page of the original tool:
https://github.com/lh3/psmc

### Usage example:

    ./vcf2psmc_coverage.pl \
        --vcf snps.vcf.gz \ # The VCF file (can be gzipped)
        --windows 100 \ # The window-resolution
        --output snps.vcf.gz.psmcfa \ # The output file
        --sample mysample \ # The name of the sample in the VCF file to use
        --coverage genome_mask_accessible_100bp_windows.fasta \
        --seqs sequences.gt_500kbp.bed \ # A BED file specifying which sequences to consider
        --verbose \ # Print some extra progress statements

The small helper script **coverage_to_windows.pl** can be used to convert a per-base binary-state genome mask (0=inaccessible; 1=accessible) in FASTA format into the necessary window-based format:

### Usage example:

    ./coverage_to_windows.pl 100 genome_mask_accessible.fasta > genome_mask_accessible_100bp_windows.fasta

Here the first argument specifies the size of the window and the second argument is the accessibility mask.

## vcf2msmc_coverage.pl

This script reads a VCF file and exports data in a format compatible with MSMC. By providing a list of samples in a simple text file, it can produce datasets with any number of haplotypes. It uses a accessibility mask to keep track on the number of accessible sites observed between SNPs.

The format is described in more detail on the page of the original tool:
https://github.com/stschiff/msmc/blob/master/guide.md

### Usage example:

    ./vcf2msmc_coverage.pl \
            --vcf snps_annotated.vcf.gz \ # A VCF file with SNPs (can be gzipped)
            --output snps_annotated.vcf.gz.out_for_msmc \ # Basename of the output file
            --group msmc_run_label=samples.csv \ # A "label" pointing to a corresponding list of samples (one sample per line)
            --coverage genome_mask_accessible.fasta \ # A binary-state genome mask (0=inaccessible; 1=accessible) to keep track on the number of accessible sites (per-base resolution)
            --seqs sequences.gt_500kbp.bed \ # A BED file specifying which sequences to consider
            --verbose \ # Print some extra progress statements

## vcf2ismc_coverage_merged.pl

This script reads multiple VCF files (using one or more lists) and exports data for use with ISMC to infer recombination rates. To this end, it produces:

1. A merged VCF file containing only the genotypes of the desired sample.
2. A corresponding set of genome sequences in FASTA format in which inaccessible sites have been masked by "N"s.

### Usage example:

 		./vcf2ismc_coverage_merged.pl \
			--files \
				vcfs_in_list1 \ # The first list of sequences and VCF files
				vcfs_in_list2 \ # A second list of sequences and VCF files
			--output merged_data \ # Basename of the output
			--sample mysample \ # The name of the sample
			--coverage genome_mask_accessible_mysample.fasta \ # Genome-mask, should be specific to this sample
			--seqs sequences.gt_500kbp.bed \ # A BED file specifying which sequences to consider
			--fasta genome.fasta \ # The main genome assembly file to get the genome sequences from
			--verbose \ # Print some extra progress statements

The lists should be tabular and the first column should contain the name of the sequence and the second column should contain the location of the VCF file that contains the SNPs for that sequence:

	seq_1	vcfs/my_data.vcf

## get_minor_major_allele_all.pl

This script and **vcf2recode_minor_vcf.pl** are used to prepare data for GEVA, a tool to infer the ages of alleles:
https://github.com/pkalbers/geva

For each population, the major allele is assumed to be ancestral and the minor allele to be derived.

This script reads a tabular pairwise FST-matrix file that also contains allele frequencies for two populations and produces tables with the major and minor allele at every SNP position. Three files are produced, specifying the major/minor allele for:

- The first population
- The second population
- Across both all data

The code to generate such a file is here: https://github.com/andreaswallberg/Ecological-Genomics-Northern-Krill/tree/main/Population_genomics-Divergence

### Usage example:

	./get_minor_major_allele_all.pl \
		sequences.tsv \ # A list of sequences to consider (one sequence per line)
		pop1 \ # The name of the first population in the file
		pop2 \ # The name of the second population in the file
		populations.allele_counts.GT.csv.divergence.matrix.tsv \ # The FST-matrix file (can be gzipped)

## vcf2recode_minor_vcf.pl

This script reads the major/minor allele table and recodes VCF files such that REF will always be the major allele and ALT the minor allele at every SNP.

	./vcf2recode_minor_vcf.pl \
		-tag pop1 \ # A tag to assign to the output file name, which is otherwise based on thename of the VCF input file
		-minor populations.allele_counts.GT.csv.divergence.matrix.tsv.major_minor.pop1.tsv \ # The file with the major/minor alleles
		-files vcfs.list \ # A file-of-files that specify the VCF files to recode


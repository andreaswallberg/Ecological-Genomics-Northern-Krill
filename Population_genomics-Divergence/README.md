# Scripts for computing population divergence

These scripts process allele count tables to estimate divergence between population using F-statistics such as FST per-SNP or per-window.

See https://github.com/andreaswallberg/Ecological-Genomics-Northern-Krill/tree/main/Population_genomics-SNP_processing for how to generate allele frequency tables.

## allele_counts2fst_matrix.reynolds.pl

This script computes FST per SNP (Weir-Cockerham) and/or per window (Reynolds)

    allele_counts2fst_matrix.reynolds.pl \
      --input populations.allele_counts.GT.csv \ # Allele counts
      --output populations.allele_counts.GT.csv.divergence \ # Basename of output
      --coverage $COVFILE \ # A genome mask (this can be a binary-state accessibility mask or one encoding different genomic regions such as CDS, intron using multiple states)
      --seqs $SUBSET \ # A bed file specifying which sequences to estimate FST values from (OPTIONAL, if not specified it will scan all sequences in the allele count file)
      --window 100 \ # Compute window-based FST (OPTIONAL)
      --groups \ # Groups together populations
        at=can,mai,bar,ice,nor,sva,swe \ # First population
        me=brc \ # Second population (there can be pairwise comparisons among more than two populations)
      --region \ # Which genomic regions to look (refers back to the states encoded in the provided genome mask)
        1=any # In a binary-state accessibility mask, only consider SNPs that fall on a position in the genome encoded by a "1"
      
 Some additional and optional arguments are:
 
      --no_print_snps 1 # This skips the printing of individual SNPs, such that only window-based FSTs are printed. This is faster and used when per-SNP estimates are not needed.
      --keep snps.tsv # A tabular file specifying which SNPs to inspect (for example, SNPs remaining after pruning for LD or those in a specific region)
      
The per-SNP output looks as follows:

    seq_s_1	1982	T	C	at/me:0.0246:148:1.0000:1.0000:1.0000|at,134,119.0000,15.0000,0.8881,0.1119|me,14,14.0000,0.0000,1.0000,0.0000
    seq_s_1	1994	G	A	at/me:0.0503:148:1.0000:1.0000:1.0000|at,134,124.0000,10.0000,0.9254,0.0746|me,14,11.0000,3.0000,0.7857,0.2143
    seq_s_1	2005	A	T	at/me:0.0158:148:1.0000:1.0000:1.0000|at,134,112.0000,22.0000,0.8358,0.1642|me,14,10.0000,4.0000,0.7143,0.2857
    seq_s_1	2012	T	C	at/me:0.0121:148:1.0000:1.0000:1.0000|at,134,111.0000,23.0000,0.8284,0.1716|me,14,10.0000,4.0000,0.7143,0.2857
    seq_s_1	2013	C	A	at/me:0.0121:148:1.0000:1.0000:1.0000|at,134,111.0000,23.0000,0.8284,0.1716|me,14,10.0000,4.0000,0.7143,0.2857
    seq_s_1	2023	A	T	at/me:0.0000:148:1.0000:1.0000:1.0000|at,134,131.0000,3.0000,0.9776,0.0224|me,14,14.0000,0.0000,1.0000,0.0000
    seq_s_1	2024	C	A	at/me:0.0000:148:1.0000:1.0000:1.0000|at,134,129.0000,5.0000,0.9627,0.0373|me,14,13.0000,1.0000,0.9286,0.0714
    seq_s_1	2027	T	A	at/me:0.0000:148:1.0000:1.0000:1.0000|at,134,129.0000,5.0000,0.9627,0.0373|me,14,13.0000,1.0000,0.9286,0.0714
    seq_s_1	2028	C	A	at/me:0.0000:148:1.0000:1.0000:1.0000|at,134,129.0000,5.0000,0.9627,0.0373|me,14,13.0000,1.0000,0.9286,0.0714
    seq_s_1	2033	T	G	at/me:0.0000:148:1.0000:1.0000:1.0000|at,134,131.0000,3.0000,0.9776,0.0224|me,14,14.0000,0.0000,1.0000,0.0000
    seq_s_1	2039	T	C	at/me:0.0000:148:1.0000:1.0000:1.0000|at,134,131.0000,3.0000,0.9776,0.0224|me,14,14.0000,0.0000,1.0000,0.0000

- The first four columns specify SNP position and ref/alt alleles
- This fifth column splits into three major subfields on "|": one about the pairwise comparison and two with metadata about each population.

      name of contrast -> at/me
      FST -> 0.0246
      Total chromosome count in dataset -> 148 (74 diploids)
      Proportion of genotyped samples across the whole dataset -> 1.0000
      Proportion of genotyped samples in population 1 -> 1.0000
      Proportion of genotyped samples in population 2 -> 1.0000

      1st subfield (at/me:0.0246:148:1.0000:1.0000:1.0000)

      name of contrast (at/me)
      FST of SNP (0.0246)
      Sample size (148)
      Proportion of observed data given overall sample size (1.0000), <1 if there are missing genotypes.
      Proportion of observed data given sample size of population 1 (1.0000)
      As above but for population 2 (1.0000)
      
      2nd and 3rd subfields (at,134,119.0000,15.0000,0.8881,0.1119 and me,14,14.0000,0.0000,1.0000,0.0000)

      name of population
      sample size
      number of observed reference alleles
      number of observed alternate alleles
      frequency of reference allele
      frequency of alternate allele

Window-based files contain FST estimates across non-overlapping windows. They include a header on the first line and two columns per population:

CHROM = name of sequence
POS = window start position
N_(contrast) = number of SNPs in the window
FST_(contrast) = average Reynold's FST of the window.

    CHROM	POS	N_at_vs_me	FST_at_vs_me
    seq_s_1	0		
    seq_s_1	1000	2	0.0838
    seq_s_1	2000	85	0.0931
    seq_s_1	3000	255	0.0929
    seq_s_1	4000	26	0.0587
    seq_s_1	5000		
    seq_s_1	6000	50	0.1053
    seq_s_1	7000	185	0.076
    seq_s_1	8000	169	0.1067
    seq_s_1	9000	30	0.0954


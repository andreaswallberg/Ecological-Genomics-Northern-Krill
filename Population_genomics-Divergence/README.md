# Scripts for computing population divergence

These scripts process allele count tables to estimate divergence between population using F-statistics such as FST per-SNP or per-window.

See https://github.com/andreaswallberg/Ecological-Genomics-Northern-Krill/tree/main/Population_genomics-SNP_processing for how to generate allele frequency tables.

## allele_counts2fst_matrix.reynolds.pl

This script computes FST per SNP (Weir-Cockerham) and/or per window (Reynolds). It sets FST values for SNPs <0 to 0.

### Usage example:

    allele_counts2fst_matrix.reynolds.pl \
      --input populations.allele_counts.GT.csv \ # Allele counts
      --output populations.allele_counts.GT.csv.divergence \ # Basename of output
      --coverage genome_mask_accessible_sites.fasta \ # A genome mask (this can be a binary-state accessibility mask or one encoding different genomic regions such as CDS, intron using multiple states)
      --seqs subset.bed \ # A bed file specifying which sequences to estimate FST values from (OPTIONAL, if not specified it will scan all sequences in the allele count file)
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

Example from the first line:

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

## allele_counts2fst_matrix.reynolds.genes.pl

This script re-uses much of the code from **allele_counts2fst_matrix.reynolds.pl** but computes the average FST (Reynolds estimator) across genes and their internal regions (e.g. CDS, exon, intron) from the allele count input. It requires a GFF file with non-redundant gene coordinates.

### Usage example:

	time allele_counts2fst_matrix.reynolds.genes.pl \
		--input populations.allele_counts.GT.csv \ # Allele counts
		--output populations.allele_counts.GT.csv.divergence.genes \ # Basename of output
		--gff genes.gff3 \ # A GFF file with gene coodinates
		--coverage genome_mask_gene_regions.fasta \ # A genome mask encoding different genomic regions such as CDS, intron using multiple states.
		--seqs $GENOMEDIR/subsets/genic_sequences/MULTIPLE.short_reads.hisat2.sorted.bam.stringtie2.merged.1x_coverage.gtf/1.m_norvegica.main_w_mito.fasta.genic_sequences.bed.80_subsets \
        --groups \ # Groups together populations
            at=can,mai,bar,ice,nor,sva,swe \ # First population
            me=brc \ # Second population (there can be pairwise comparisons among more than two populations)
		--region \ # Which genomic regions to look (refers back to the states encoded in the provided genome mask)
			1=intergenic \
			2=intron \
			3=three_prime_utr \
			4=exon \
			5=five_prime_utr \
			6=cds \
		--no_print_snps 1 
        
In addition to printing per-SNP and per-window outputs, this script also prints per-gene output for every contrast.
 
The per-gene output looks as follows:

    NR	CHROM	GENE	START	STOP	LENGTH	FLANKING_FST	FLANKING_N	GENE_FST	GENE_N	EXON_FST	EXON_N	CDS_FST	CDS_N	INTRON_FST	INTRON_N	UTR_FST	UTR_N
    1	seq_s_1	REF_TRIN_14_02212_XLOC_012882	568696	606697	38002	0.121876054563044	501	0.0695167578182325	10	0.0695167578182325	10	0.0695167578182325	10				
    2	seq_s_1	REF_STRG_1_4_XLOC_012878	445986	625695	179710	0.0691069254600813	974	0.166466056484738	5125	0.475972915450141	24	0.475972915450141	24	0.164544484297914	5101		
    3	seq_s_1	REF_STRG_1_6_XLOC_012891	626372	627898	1527	0.081422662300597	124										
    4	seq_s_1	REF_TRIN_9_00309_XLOC_012885	866836	1053628	186793	0.0409807815995226	4176	0.0466079341456129	7398	0.0403659153090983	46	0.0403659153090983	46	0.0466378843826125	7352		
    5	seq_s_1	REF_STRG_1_13_XLOC_012886	1566795	1673962	107168	0.0402054060614511	2516	0.0584927303984401	2060	0.0288981554016114	39	0.0442113954736102	5	0.0588392026475325	2021	0.0231608221574779	34
    6	seq_s_1	REF_TRIN_7_01050_XLOC_012887	2066169	2066579	411	0.0449776682579952	4230	0.019891981638239	36	0.019891981638239	36	0.019891981638239	36				
    7	seq_s_1	REF_STRG_1_2_XLOC_012894	2065701	2068097	2397			0.034228264801881	262	0.034228264801881	262	0.0272115684790773	117			0.0412927310438469	145
    8	seq_s_1	REF_STRG_1_3_XLOC_012888	2071257	2084957	13701	0.0634939161703186	1492	0.0510506102953343	928	0.019699743455261	10	0.0199177489541814	8	0.0513751334604216	918	0.00909014250779997	2
    9	seq_s_1	REF_STRG_1_119_XLOC_012889	2340340	2529430	189091	0.0509611700876408	4579	0.0546942248756524	5636					0.0546942248756524	5636		
    10	seq_s_93	REF_STRG_1_1166_XLOC_112702	14781	711797	697017			0.0702911953986204	16063	0.0265099562719092	16	0.0419223853472152	4	0.0703155251942406	16047	0.0207186996030281	12

Flanking FST is computed from in the region 50-100 kbp upstream and downstream of a gene (hardcoded into the script).



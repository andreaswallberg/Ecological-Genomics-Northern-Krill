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

## vcf2allele_counts.pl

This script estimates SNP allele frequencies for populations by grouping together samples present in a VCF file.

The allele frequencies can be based on directly on the called genotypes ("GT") or on the genotype likelihoods ("GL") or on the genotype probablities ("GP").

The output is a tabular file that specifies the allele counts for every population at every SNP site. It generates two output files per "method" (GT, GL or GP), one containing the allele counts and the other containing the allele frequencies. The GL and GP methods incorporate some level of uncertainty in the genotyping.

The main advantage of this format is that it represents a condensed version of the variation data that is faster to parse downstream of this point compared to repeatedly parsing the VCF.

Usage example:

        time vcf2allele_counts.pl \
        --input phased_imputed_snps.vcf.gz \ # Input VCFs
        --output phased_imputed_snps.vcf.gz.alleles \ # Basename of output
        --gt \ # Output based on genotypes (0/0, 0/1, 1/1)
        --gp \ # Output based on genotype probabilities
        --groups \
            bar=samples.bar.csv \
            brc=samples.brc.csv \
            can=samples.can.csv \
            ice=samples.ice.csv \
            mai=samples.mai.csv \
            nor=samples.nor.csv \
            sva=samples.sva.csv \
            swe=samples.swe.csv
                
 The "groups" argument specifies the label of the group/population and points to a csv file with the names of the samples for that group (one sample per line).
 
 The output for the "GT" "allele_count" file looks as follows for 10 SNPs:
 
    CHROM	POS	REF	ALT	N_WITH_GT	TOT	TOT	bar	bar	brc	brc	can	can	ice	ice	mai	mai	nor	nor	sva	sva	swe	swe
    seq_s_1	1982	T	C	74	133	15	18	2	14	0	19	1	15	5	20	0	18	2	15	3	14	2
    seq_s_1	1994	G	A	74	135	13	18	2	11	3	19	1	19	1	18	2	20	0	16	2	14	2
    seq_s_1	2005	A	T	74	122	26	15	5	10	4	20	0	15	5	16	4	15	5	18	0	13	3
    seq_s_1	2012	T	C	74	121	27	15	5	10	4	19	1	15	5	16	4	15	5	18	0	13	3
    seq_s_1	2013	C	A	74	121	27	15	5	10	4	19	1	15	5	16	4	15	5	18	0	13	3
    seq_s_1	2023	A	T	74	145	3	20	0	14	0	17	3	20	0	20	0	20	0	18	0	16	0
    seq_s_1	2024	C	A	74	142	6	18	2	13	1	20	0	19	1	19	1	20	0	17	1	16	0
    seq_s_1	2027	T	A	74	142	6	18	2	13	1	20	0	19	1	19	1	20	0	17	1	16	0
    seq_s_1	2028	C	A	74	142	6	18	2	13	1	20	0	19	1	19	1	20	0	17	1	16	0
    seq_s_1	2033	T	G	74	145	3	20	0	14	0	17	3	20	0	20	0	20	0	18	0	16	0

Columns 1-4 specify the site and the ref/alt alleles. The 5th column specifies the number of genotyped individuals across the whole dataset. The 6th and 7th columns specify the counts for the reference and alternate alleles across the total dataset ("TOT"), followed by the counts for each population/group.

The corresponding "GT" "allele_freqs" file looks as follows:

    CHROM	POS	REF	ALT	N_WITH_GT	TOT	TOT	bar	bar	brc	brc	can	can	ice	ice	mai	mai	nor	nor	sva	sva	swe	swe
    seq_s_1	1982	T	C	74	0.89865	0.10135	0.9	0.1	1	0	0.95	0.05	0.75	0.25	1	0	0.9	0.1	0.83333	0.16667	0.875	0.125
    seq_s_1	1994	G	A	74	0.91216	0.08784	0.9	0.1	0.78571	0.21429	0.95	0.05	0.95	0.05	0.9	0.1	1	0	0.88889	0.11111	0.875	0.125
    seq_s_1	2005	A	T	74	0.82432	0.17568	0.75	0.25	0.71429	0.28571	1	0	0.75	0.25	0.8	0.2	0.75	0.25	1	0	0.8125	0.1875
    seq_s_1	2012	T	C	74	0.81757	0.18243	0.75	0.25	0.71429	0.28571	0.95	0.05	0.75	0.25	0.8	0.2	0.75	0.25	1	0	0.8125	0.1875
    seq_s_1	2013	C	A	74	0.81757	0.18243	0.75	0.25	0.71429	0.28571	0.95	0.05	0.75	0.25	0.8	0.2	0.75	0.25	1	0	0.8125	0.1875
    seq_s_1	2023	A	T	74	0.97973	0.02027	1	0	1	0	0.85	0.15	1	0	1	0	1	0	1	0	1	0
    seq_s_1	2024	C	A	74	0.95946	0.04054	0.9	0.1	0.92857	0.07143	1	0	0.95	0.05	0.95	0.05	1	0	0.94444	0.05556	1	0
    seq_s_1	2027	T	A	74	0.95946	0.04054	0.9	0.1	0.92857	0.07143	1	0	0.95	0.05	0.95	0.05	1	0	0.94444	0.05556	1	0
    seq_s_1	2028	C	A	74	0.95946	0.04054	0.9	0.1	0.92857	0.07143	1	0	0.95	0.05	0.95	0.05	1	0	0.94444	0.05556	1	0
    seq_s_1	2033	T	G	74	0.97973	0.02027	1	0	1	0	0.85	0.15	1	0	1	0	1	0	1	0	1	0

The corresponding "GP" "allele_count" file looks like follows:

    CHROM	POS	REF	ALT	N_WITH_GT	TOT	TOT	bar	bar	brc	brc	can	can	ice	ice	mai	mai	nor	nor	sva	sva	swe	swe
    seq_s_1	1982	T	C	74	121.006	26.994	15.922	4.078	12.955	1.047	16.905	3.095	13.599	6.401	16.828	3.172	16.493	3.507	15.314	2.686	12.99	3.008
    seq_s_1	1994	G	A	74	118.225	29.723	15.577	4.421	10.329	3.669	16.693	3.299	16.677	3.313	15.373	4.615	17.863	2.131	13.468	4.526	12.245	3.749
    seq_s_1	2005	A	T	74	114.551	33.441	14.151	5.849	10.034	3.966	18.921	1.079	13.526	6.474	14.715	5.281	13.922	6.074	16.753	1.249	12.529	3.469
    seq_s_1	2012	T	C	74	113.088	34.908	13.959	6.041	9.987	4.013	18.708	1.29	13.32	6.68	14.507	5.491	13.594	6.408	16.614	1.384	12.399	3.601
    seq_s_1	2013	C	A	74	113.339	34.661	13.981	6.019	9.99	4.01	18.76	1.238	13.354	6.646	14.543	5.457	13.642	6.358	16.652	1.35	12.417	3.583
    seq_s_1	2023	A	T	74	143.455	4.549	19.817	0.183	14	0	16.684	3.318	19.752	0.25	19.71	0.29	19.855	0.145	17.776	0.224	15.861	0.139
    seq_s_1	2024	C	A	74	137.149	10.847	17.108	2.892	12.999	1.001	19.442	0.556	18.27	1.73	18.252	1.748	19.595	0.403	15.919	2.081	15.564	0.436
    seq_s_1	2027	T	A	74	137.356	10.64	17.135	2.861	13	1	19.479	0.521	18.298	1.702	18.283	1.717	19.621	0.379	15.956	2.044	15.584	0.416
    seq_s_1	2028	C	A	74	137.208	10.79	17.114	2.888	13	1	19.455	0.545	18.274	1.724	18.259	1.739	19.606	0.394	15.93	2.07	15.57	0.43
    seq_s_1	2033	T	G	74	143.578	4.43	19.83	0.17	14	0	16.702	3.298	19.771	0.229	19.73	0.272	19.87	0.132	17.801	0.203	15.874	0.126

The corresponding "GP" "allele_freqs" file looks like follows:

    CHROM	POS	REF	ALT	N_WITH_GT	TOT	TOT	bar	bar	brc	brc	can	can	ice	ice	mai	mai	nor	nor	sva	sva	swe	swe
    seq_s_1	1982	T	C	74	0.81761	0.18239	0.7961	0.2039	0.92522	0.07478	0.84525	0.15475	0.67995	0.32005	0.8414	0.1586	0.82465	0.17535	0.85078	0.14922	0.81198	0.18802
    seq_s_1	1994	G	A	74	0.7991	0.2009	0.77893	0.22107	0.73789	0.26211	0.83498	0.16502	0.83427	0.16573	0.76911	0.23089	0.89342	0.10658	0.74847	0.25153	0.7656	0.2344
    seq_s_1	2005	A	T	74	0.77404	0.22596	0.70755	0.29245	0.71671	0.28329	0.94605	0.05395	0.6763	0.3237	0.7359	0.2641	0.69624	0.30376	0.93062	0.06938	0.78316	0.21684
    seq_s_1	2012	T	C	74	0.76413	0.23587	0.69795	0.30205	0.71336	0.28664	0.93549	0.06451	0.666	0.334	0.72542	0.27458	0.67963	0.32037	0.9231	0.0769	0.77494	0.22506
    seq_s_1	2013	C	A	74	0.7658	0.2342	0.69905	0.30095	0.71357	0.28643	0.93809	0.06191	0.6677	0.3323	0.72715	0.27285	0.6821	0.3179	0.92501	0.07499	0.77606	0.22394
    seq_s_1	2023	A	T	74	0.96926	0.03074	0.99085	0.00915	1	0	0.83412	0.16588	0.9875	0.0125	0.9855	0.0145	0.99275	0.00725	0.98756	0.01244	0.99131	0.00869
    seq_s_1	2024	C	A	74	0.92671	0.07329	0.8554	0.1446	0.9285	0.0715	0.9722	0.0278	0.9135	0.0865	0.9126	0.0874	0.97985	0.02015	0.88439	0.11561	0.97275	0.02725
    seq_s_1	2027	T	A	74	0.92811	0.07189	0.85692	0.14308	0.92857	0.07143	0.97395	0.02605	0.9149	0.0851	0.91415	0.08585	0.98105	0.01895	0.88644	0.11356	0.974	0.026
    seq_s_1	2028	C	A	74	0.92709	0.07291	0.85561	0.14439	0.92857	0.07143	0.97275	0.02725	0.91379	0.08621	0.91304	0.08696	0.9803	0.0197	0.885	0.115	0.97313	0.02687
    seq_s_1	2033	T	G	74	0.97007	0.02993	0.9915	0.0085	1	0	0.8351	0.1649	0.98855	0.01145	0.9864	0.0136	0.9934	0.0066	0.98872	0.01128	0.99213	0.00788

 

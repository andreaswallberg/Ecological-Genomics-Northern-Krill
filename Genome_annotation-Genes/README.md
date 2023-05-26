# Scripts for processing gene annotations

Note: These scripts processes GTF/GFF files derived from GFFCOMPARE as the upstream tool, and may not be compatible with other GTFs/GFFs.

## get_best_blast_hit_for_locus.pl

This script parses tabular BLAST / DIAMOND output ("-outfmt 6") and generates a list with the best BLAST hits per locus. It assumes the common locus tag shared among accessions, e.g. among isoforms, are part of the accession name and preceeds the first dot:

        XLOC_000022.1.p1	XP_027237621.1	77.9	810	172	4	11	820	3	805	0.0e+00	1239.9
        XLOC_000022.1.p1	XP_037775203.1	76.7	810	182	4	11	820	3	805	0.0e+00	1220.7

Here "XLOC_000022" is the common locus.

### Usage example:

        ./get_best_blast_hit_for_locus.pl BLAST.out > BLAST.out.best.csv

## keep_best_isoform.pl

This scripts works downstream of **get_best_blast_hit_for_locus.pl** and filters a GTF file based on the table of "best" isoforms, i.e. it produces a non-redundant set of annotations.

### Usage example:

    .././keep_best_isoform.pl \
        BLAST.out.best.csv \
        genes.gtf \
        1> genes_best.gtf \
        2> genes_best.tsv

## parse_tracking.pl

This script has a narrow and study-specific scope. It parses all loci detected with GFFCOMPARE by reading the ".tracking" file and keeps only loci that fulfil the criteria below. It was designed to provide a set of well-supported comparative gene evidence derived from mapping multiple species to the genome.

To keep a locus, it must be supported by:
- At least two genes with two or more exons that together are >= 500 bp long
- A signicant BLAST/DIAMOND hit with a score >= 100
- Detected in two or species

### Usage example:

    ./parse_tracking.pl \
        gffcompare.tracking \
        BLAST.outfmt6

The species names are hard-coded into the script, making this script specific to this study.

Downstream of this script, **get_and_rename_transcripts.pl** (also study-specific) processes and renames transcripts in GTFs/GFFs:

### Usage example:

    ./get_and_rename_transcripts.pl \
        COMPARATIVE \
        gffcompare.tracking.keep_transcripts.csv \
        genes.gff3
        
       ./get_and_rename_transcripts.pl \
        COMPARATIVE \
        gffcompare.tracking.keep_transcripts.csv \
        genes.gtf

## rename_multiple_transcripts_on_chosen.pl

This script has a narrow and study-specific scope. It relabels all isoforms belonging to the same locus based on the "best" isoform identified above. It preserves the original evidence for each isoform in the "oId" tag in the GTF file, using the following labels:

    "REF_STRG_*" => a locus with the "best" gene model derived from expressed Illumina RNA-seq or Nanopore cDNA (RNA evidence)
    "REF_TRIN_*" => a locus with the "best" gene model derived assembled and mapped Trinity RNA-seq transcripts (RNA evidence)
    "COM_SPAL_*" => a locus with the "best" gene model derived mapped gene models from other species (comparative evidence)

### Usage example:

    ./rename_multiple_transcripts_on_chosen.pl \
        genes_best.tsv \
        genes_best.gtf \
        genes_best.gtf.relabeled.gtf



## gff_update_stop_codons.pl

This script adds missing stop codons to coding sequences in GFF files if they are supported by the genome sequence, i.e. if they occur just after the end of the open reading frame.

It requires providing genome sequences and CDS sequences in separate FASTA files, as well as a GFF file:

### Usage example:

    ./gff_update_stop_codons.pl \
        genome.fasta \
        genes.gff3 \
        cds.fa

## fix_coordinates.pl

A small script used to fix coordinates such that start coordinates for features in a GFF are always smaller than stop coordinates, regardless of strand orientation of the feature.

        ./fix_coordinates.pl my_file.gff > my_file_fixed.gff

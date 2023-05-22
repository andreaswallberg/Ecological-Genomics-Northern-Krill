# Scripts for processing gene annotations

## gff_update_stop_codons.pl

This script adds missing stop codons to coding sequences in GFF files if they are supported by the genome sequence, i.e. if they occur just after the end of the open reading frame.

It requires providing genome sequences and CDS sequences in separate FASTA files, as well as a GFF file:

### Usage example:

    gff_update_stop_codons.pl \
        genome.fasta \
        genes.gff3 \
        cds.fa


# Scripts for basic processing tasks of genome sequences

## split_fasta_chunks.pl

This script splits a FASTA file into multiple files of a given size. It reads the information on STDIN. The script does not split sequences so the actual output size of chunks may vary depending on the lengths of input sequences.

### Usage example:

    cat genome_file.fasta | ./split_fasta_chunks.pl genome_file_split 100000000

The first argument is the basename of the output file(s) and the second argument is the chunk size of the outputs.

## rewrite_fastq_chromium.pl

This script reads 10X Chromium barcoded reads in FASTQ format and appends a string to the BX:Z barcode tag

### Usage example:

    ./rewrite_fastq_chromium.pl <string> reads.fastq

The first argument is the string and the second argument is the FASTQ file. The file needs to be decompressed. Alterantively, the script can read a decompressed stream of sequence data on STDIN. It prints to STDOUT.

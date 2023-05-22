# Scripts for basic processing tasks of genome sequences

## split_fasta_chunks.pl

This script splits a FASTA file into multiple files of a given size (bp). It reads the information on STDIN. The script does not split sequences so the actual output size of chunks may vary depending on the lengths of input sequences.

### Usage example:

    cat genome_file.fasta | ./split_fasta_chunks.pl genome_file_split 100000000

The first argument is the basename of the output file(s) and the second argument is the chunk size of the outputs.

## split_fasta_N_files.pl

This script splits a FASTA file into multiple files of a given number (n). The script does not split sequences so the actual output size of chunks may vary depending on the lengths of input sequences.

### Usage example:

    split_fasta_N_files.pl 10 genome_file.fasta

The first argument is the number of output file(s) and the second argument is the input genome file. The output file names are derived from the input.

## rewrite_fastq_chromium.pl

This script reads 10X Chromium barcoded reads in FASTQ format and appends a string to the BX:Z barcode tag

### Usage example:

    ./rewrite_fastq_chromium.pl <string> reads.fastq

The first argument is the string and the second argument is the FASTQ file. The file needs to be decompressed. Alterantively, the script can read a decompressed stream of sequence data on STDIN. It prints to STDOUT.

## chromium_barcode_analysis_simple.pl

This script filters a BAM file with the purpose to keep only edge reads that map within some distance to the edge of sequences.

### Usage example:

        time ./chromium_barcode_analysis_simple.pl \
            -mq 20 \
            -md 20000 \
            -rl 142 \
            -bam short_reads.bam
Here, "-mq" is mapping quality, "-md" is mapping distance and "-rl" is the average read length.

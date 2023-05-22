#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;

# Number of chunks to make

my $out_n = shift @ARGV;

# Base chunk name on first FASTA file

my $out_name = $ARGV[ 0 ];

# Total sequence length

my $seq_length = 0;

# Scan all FASTA files for the total sequence size

foreach my $fasta_file ( @ARGV ) {

    open ( my $in , "<" , $fasta_file ) or die "Unable to open $fasta_file: $!";
    
    while (<$in>) {
    
        chomp;
        
        unless ( $_ =~ m/^\>/ ) {
        
            $seq_length += length $_;
        
        }
    
    }
    
    close $in;

}

# Get the chunk size

my $out_size = $seq_length / $out_n;

# Output file number

my $i = 0;

# Main FASTA chunk output

my $out;

# Companion list files

my $list_out;
my $list_csv_out;

# Temporary size counter that is reset when we reach the size limit for a chunk

my $tmp_size = 0;

# Scan all FASTA files again

foreach my $fasta_file ( @ARGV ) {

    open ( my $in , "<" , $fasta_file ) or die "Unable to open $fasta_file: $!";

    while (<$in>) {

        chomp;
        
        if ( $_ =~ m/^\>(.+)/ ) {
        
            # If we have passed the size limit, or we are at the beginning
            # of the process
        
            if ( $tmp_size >= $out_size or not defined $out ) {
            
                $i++;
            
                # DEBUG: print STDERR "Starting output file $i: ${out_name}.${i}.fasta\n";
            
                open ( $out , ">" , "${out_name}.${i}.fasta" ) or die "$!";
                open ( $list_out , ">" , "${out_name}.${i}.fasta.list" ) or die "$!";
                open ( $list_csv_out , ">" , "${out_name}.${i}.fasta.list.csv" ) or die "$!";
            
                $tmp_size = 0;
            
            }
            
            print $out ">$1\n";
            
            # Insert a space if this is not the first sequence seen for this chunk
            
            if ( $tmp_size > 0 ) {
            
                print $list_out " ";
            
            }
            
            print $list_out "$1";
            
            print $list_csv_out "$1\n";
            
        }
        
        else {
        
            print $out "$_\n";
            
            $tmp_size += length $_;
        
        }

    }
    
    close $in;

}

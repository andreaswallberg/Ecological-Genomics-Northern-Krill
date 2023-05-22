#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $best_file = shift @ARGV;

my $best_table = {};

open ( my $best_in , "<" , $best_file ) or die "$!";

while (<$best_in>) {

    chomp;
    
    my ( $gene , $transcript , $score ) = split ( /\t/ , $_ );
    
    $transcript =~ s/\.p\d+$//;
    
    $best_table->{ $gene } = $transcript;

}

my $gff_file = shift @ARGV;

my $exon_table = {};

open ( my $gff_in , "<" , $gff_file ) or die "$!";

while (<$gff_in>) {

	chomp;
	
    # seq_s_68404	transdecoder	transcript	16777	53953	.	-	.	transcript_id "REFERENCE_00070056"; gene_id "XLOC_036425"; oId "STRG.ref_mixed.95802.1.p1"; tss_id "TSS58692"; num_samples "1";
	
	if ( $_ =~ m/\sexon\s(\d+)\s(\d+).+transcript_id\s\"([^\"]+)\"\;\sgene_id\s\"([^\"]+)\"\;/ ) {
	
		my ( $start , $stop , $isoform , $gene ) = ( $1 , $2 , $3 , $4 );
		
		if ( $start > $stop ) {
		
            ( $start , $stop ) = ( $stop , $start );
		
		}
		
		my $length = $stop - $start + 1;

		$exon_table->{ $gene }->{ $isoform } += $length;
	
	}

}

close $gff_in;

my $longest_table = {};

foreach my $gene ( sort keys %{ $exon_table } ) {

    my @isoforms = keys %{ $exon_table->{ $gene } };
    
    my @sorted_isoforms = sort {
        $exon_table->{ $gene }->{ $b } <=>
        $exon_table->{ $gene }->{ $a }
    } @isoforms;
    
    $longest_table->{ $sorted_isoforms[ 0 ] } = 
    $exon_table->{ $gene }->{ $sorted_isoforms[ 0 ] };

}

my $to_print = 0;

open ( $gff_in , "<" , $gff_file ) or die "$!";

while (<$gff_in>) {

	chomp;
	
	if ( $_ =~ m/\stranscript\s(\d+)\s(\d+).+transcript_id\s\"([^\"]+)\"\;\sgene_id\s\"([^\"]+)\"\;/ ) {
	
		my ( $start , $stop , $isoform , $gene ) = ( $1 , $2 , $3 , $4 );
		
		$to_print = 0;
		
		if ( defined $best_table->{ $gene } ) {
		
            if ( $isoform eq $best_table->{ $gene } ) {
            
                say STDERR "$gene\t$isoform\tbest";
            
                $to_print = 1;
            
            }
		
		}
		
		elsif ( defined $longest_table->{ $isoform } ) {
		
            say STDERR "$gene\t$isoform\tlongest";
		
			$to_print = 1;
		
		}
	
	}
	
	if ( $to_print ) {
	
		say $_;
	
	}

}

#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $gff_file = shift @ARGV;
my $repeat_file = shift @ARGV;

my @ids = sort { $b <=> $a } @ARGV;

open ( my $gff_in , "<" , $gff_file ) or die "$!";

my $repeat_table = {};

while (<$gff_in>) {

	chomp;
	
	# seq_c_19089	LTR_retriever	repeat_region	259605	262296	.	?	.	ID=repeat_region_1;Name=seq_c_19089:259608..262293;Classification=LTR/unknown;Sequence_ontology=SO:0000657;ltr_identity=1.0000;Method=structural;motif=TGCA;tsd=TTT

	if ( $_ =~ m/^(\S+)\s\S+\s\S+\s(\d+)\s+(\d+)\s+.+ltr_identity\=([^\;]+)/ ) {
	
		my $repeat = $1 . "_" . $2 . "_" . $3;
	
		$repeat_table->{ $repeat } = $4;
	
		# say "Adding $repeat $4";
	
	}

}

open ( my $repeat_in , "<" , $repeat_file ) or die "$!";

my $out_table = {};

my $repeat;
my $out;

while (<$repeat_in>) {

	chomp;
	
	if ( $_ =~ m/^\>/ ) {
	
		$repeat = undef;
	
	}
	
	if ( $_ =~ m/^\>(\S+)\_(c|s)\_(\d+)\_(\d+)\_(\d+)/ ) {
	
		$repeat = $1 . "_" . $2 . "_" . $3 . "_" . $4 . "_" . $5;
	
		say $repeat;
	
		$out = undef;
	
		if ( defined $repeat_table->{ $repeat } ) {
	
			my $id = $repeat_table->{ $repeat };
			
			say $id;
			
			foreach my $tmp_id ( @ids ) {
			
				if ( $id >= $tmp_id ) {
				
					if ( defined $out_table->{ $tmp_id }->{ 'out' } ) {
					
						$out = $out_table->{ $tmp_id }->{ 'out' };
					
					}
					
					else {
					
						open ( $out , ">" , $repeat_file . ".id_${tmp_id}.fasta" ) or die "$!";
					
						$out_table->{ $tmp_id }->{ 'out' } = $out;
					
					}
					
					say "$repeat with id $id goes to $tmp_id ...";
				
					say $out_table->{ $tmp_id }->{ 'n' }++;
				
					last;
				
				}
			
			}
	
		}
	
	}
	
	if ( defined $repeat and defined $out ) {
	
		say $out $_;
	
	}
	
	else {
	
		die $_;
	
	}

}

foreach my $tmp_id ( @ids ) {

	if ( defined $out_table->{ $tmp_id }->{ 'n' } ) {
	
		say "Repeat identity class $tmp_id n=" , $out_table->{ $tmp_id }->{ 'n' };
	
	}

}

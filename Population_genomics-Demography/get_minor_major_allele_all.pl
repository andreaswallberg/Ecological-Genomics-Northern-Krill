#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $chrom_table = {};

my $chrom_file = shift @ARGV;

open ( my $chrom_in , "<" , $chrom_file ) or die "$!";

while (<$chrom_in>) {

	chomp;
	
	$chrom_table->{ $_ } = 1;

	# say $_;
	
}

my $first_pop = shift @ARGV;
my $second_pop = shift @ARGV;

my $snp_file = shift @ARGV;

my $in;

if ( $snp_file =~ m/\.gz$/ ) {

	open ( $in , "zcat $snp_file |" ) or die "$!";

}

else {

	open ( $in , "<" , $snp_file ) or die "$!";

}

open ( my $first_out , ">" , "${snp_file}.major_minor.${first_pop}.tsv" ) or die "$!";

say $first_out "CHROM\tPOS\tMAJOR\tMINOR";

open ( my $second_out , ">" , "${snp_file}.major_minor.${second_pop}.tsv" ) or die "$!";

say $second_out "CHROM\tPOS\tMAJOR\tMINOR";

open ( my $all_out , ">" , "${snp_file}.major_minor.all.tsv" ) or die "$!";

say $all_out "CHROM\tPOS\tMAJOR\tMINOR";

while (<$in>) {

	chomp;

	if ( $_ =~ m/^(\S+)/ ) {
	
		next unless defined $chrom_table->{ $1 };
	
	}
	
	if ( $_ =~ m/^(\S+)\s(\d+)\s(\S)\s(\S)\s\S+\|\S+\,\d+\,([^\,]+)\,([^\,]+)\S+\|\S+\,\d+\,([^\,]+)\,([^\,]+)/ ) {
	
		my ( $chrom , $pos , $ref , $alt , $first_ref , $first_alt , $second_ref , $second_alt ) = ( $1 , $2 , $3 , $4 , $5 , $6 , $7 , $8 );
		
		# say "$chrom , $pos , $ref , $alt , $first_ref , $first_alt , $second_ref , $second_alt";
	
		if ( $first_ref >= $first_alt ) {
		
			say $first_out "$chrom\t$pos\t$ref\t$alt";
		
		}
		
		elsif ( $first_alt > $first_ref ) {
		
			say $first_out "$chrom\t$pos\t$alt\t$ref";
		
		}
		
		if ( $second_ref >= $second_alt ) {
		
			say $second_out "$chrom\t$pos\t$ref\t$alt";
		
		}
		
		elsif ( $second_alt > $second_ref ) {
		
			say $second_out "$chrom\t$pos\t$alt\t$ref";
		
		}
		
		my $all_ref = $first_ref + $second_ref;
		my $all_alt = $first_alt + $second_alt;
		
		if ( $all_ref >= $all_alt ) {
		
			say $all_out "$chrom\t$pos\t$ref\t$alt";
		
		}
		
		elsif ( $all_alt > $all_ref ) {
		
			say $all_out "$chrom\t$pos\t$alt\t$ref";
		
		}
	
	}
	
}

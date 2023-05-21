#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $out_name = shift @ARGV;
my $out_size = shift @ARGV;

my $i = 0;

my $out;
my $list_out;
my $list_csv_out;

my $tmp_size = 0;

while (<>) {

	chomp;
	
	if ( $_ =~ m/^\>(.+)/ ) {
	
		if ( $tmp_size >= $out_size or not defined $out ) {
# 		
# 			if ( defined $list_out ) {
# 			
# 				say $list_out "";
# 			
# 			}
		
			$i++;
		
			say "Starting output file $i: ${out_name}.${i}.fasta";
		
			open ( $out , ">" , "${out_name}.${i}.fasta" ) or die "$!";
			open ( $list_out , ">" , "${out_name}.${i}.fasta.list" ) or die "$!";
			open ( $list_csv_out , ">" , "${out_name}.${i}.fasta.list.csv" ) or die "$!";
		
			$tmp_size = 0;
		
		}
		
		say $out ">$1";
		
		if ( $tmp_size > 0 ) {
		
			print $list_out " ";
		
		}
		
		print $list_out "$1";
		
		say $list_csv_out $1;
		
	}
	
	else {
	
		say $out $_;
		
		$tmp_size += length $_;
	
	}

}

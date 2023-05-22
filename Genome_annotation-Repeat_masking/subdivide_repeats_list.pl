#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $tag = shift @ARGV;
my $other_tag = shift @ARGV;
my $list_file = shift @ARGV;
my $repeat_file = shift @ARGV;

open ( my $list_in , "<" , $list_file ) or die "$!";

my $repeat_table = {};

while (<$list_in>) {

	chomp;
	
	$repeat_table->{ $_ } = 1;

}

open ( my $repeat_in , "<" , $repeat_file ) or die "$!";

open ( my $t_out , ">" , $repeat_file . ".${tag}.fasta" ) or die "$!";
open ( my $o_out , ">" , $repeat_file . ".${other_tag}.fasta" ) or die "$!";

my $out;

while (<$repeat_in>) {

	chomp;
	
	if ( $_ =~ m/^\>(.+)/ ) {
	
		if ( defined $repeat_table->{ $1 } ) {
		
			$out = $t_out;
		
		}
		
		else {
		
			$out = $o_out;
		
		}
	
	}
	
	say $out $_;

}


#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $win_size = shift @ARGV;

while (<>) {

	chomp;

	if ( $_ =~ m/^\>/ ) {
	
		print STDERR "$_";
	
		say $_;
	
	}
	
	else {
	
		my $seq = $_;
		
		my $length = length $seq;
		
		print STDERR "\t" , $length;
		say STDERR "\t" , int ( $length / $win_size );
		
		for ( my $i = 0; $i < $length; $i += $win_size ) {
		
			my $subseq = substr ( $seq , $i , $win_size );
			
			if ( $subseq =~ m/^0+$/ ) {
			
				print "N";
			
			}
			
			else {
			
				print "T";
			
			}
		
		}
		
		say "";
	
	}
	
}

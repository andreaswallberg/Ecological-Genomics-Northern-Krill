#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $tag = shift @ARGV;

while (<>) {

	chomp;
	
	if ( $_ =~ m/^\>(.+)/ ) {
	
		my $header = $1;
		$header =~ s/\#/${tag}\#/;
	
		$header = substr ( $header , 0 , 49 );
	
		say ">$header";
	
	}
	
	else {
	
		say $_;
	
	}

}

#!/usr/bin/perl -w

use strict;
use warnings;
use 5.010;

my $tag_num = shift @ARGV;

while (<>) {

	chomp;
	
	if ( $_ =~ m/(^\@.+\:\S+\:\S+\sBX\:Z\:\S+)\-\d+.*$/ ) {
	
		my $header = $1;
		
		say $header , "-${tag_num}";
	
	}
	
	if ( $_ =~ m/(^\@.+\:\S+\:\S+\sBX\:Z\:\S+)$/ ) {
	
		my $header = $1;
		
		say $header , "-${tag_num}";
	
	}
	
	else {
	
		say $_;
	
	}
	
}

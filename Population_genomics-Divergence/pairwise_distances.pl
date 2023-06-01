#!/usr/bin/perl

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $dist_table = {};
my $seq_table = {};

my @headers;

my $out_base = shift @ARGV;

open ( my $out , ">" , $out_base . ".distance.phy" ) or die "$!";
open ( my $i_out , ">" , $out_base . ".distance.phy.info" ) or die "$!";

foreach my $fasta_file ( @ARGV ) {

	my $in;
	
	if ( $fasta_file =~ m/\.gz$/ ) {

		open ( $in , "zcat $fasta_file |" ) or die "$!";
	
	}
	
	else {
	
		open ( $in , "<" , $fasta_file ) or die "$!";

	}
	
	say STDERR "Reading data from $fasta_file ...";
	
	my $tmp_header;
	
	while (<$in>) {

		chomp;

		if ( $_ =~ m/^\>(.+)/ ) {
		
			$tmp_header = $1;
		
			say STDERR "\tReading sequene for $tmp_header ...";
		
			unless ( defined $seq_table->{ $tmp_header } ) {
			
				push @headers , $tmp_header;
			
			}
		
		}
		
		else {
		
			my $seq = uc $_;
			$seq =~ s/[NRYSWJNBDHV\.\?\-]/Z/g;
		
			$seq_table->{ $tmp_header } .= $seq;
		
		}
		
	}

}

compute_pairwise_distances();
print_distances();

sub compute_pairwise_distances {

	my $sum;
	my $n;
	my $length_seq;
	my $sum_aligned;
	
	for ( my $i = 0; $i < @headers - 1; $i++ ) {

		my $header = $headers[ $i ];
		
		say STDERR "Comparing $header (" , $i + 1 , ") to all others ...";
		
		my $seq = $seq_table->{ $header };
		
		my $n_unknown = $seq =~ tr/Z//;
		
		$length_seq = length $seq;
		
		for ( my $j = $i + 1; $j < @headers; $j++ ) {
		
			my $other_header = $headers[ $j ];
			
			say STDERR "\t...$other_header";
			
			my $tmp_seq = $seq;
			
			my $other_seq = $seq_table->{ $other_header };
			
			my $tmp_other_seq = $other_seq;
			
			my $other_n_unknown = $tmp_other_seq =~ tr/Z//;
			
			my $length_other_seq = length $tmp_other_seq;

			# Bitwise operators
			
			my $diff = $tmp_seq ^ $tmp_other_seq;

			( my $mask = $diff ) =~ tr{\x00}{\xff}c;

			$tmp_seq &= $mask;
			$tmp_other_seq &= $mask;
					
			my $n_differences = $tmp_seq =~ tr/[ACGT]//;
			my $n_unknown_diff_other = $tmp_other_seq =~ tr/Z//;

			my $tot = $length_seq - $n_unknown - $n_unknown_diff_other;
			
			my $distance = ( $n_differences - $n_unknown_diff_other ) / $tot;

			say "$header\t$other_header\t$length_seq\t$n_unknown\t$tot\t$n_differences\t$n_unknown_diff_other\t$distance";
					
			$dist_table->{ $header }->{ $other_header } = $distance;
			$dist_table->{ $other_header }->{ $header } = $distance;
			
			$sum += $distance * 2;
			$n += 2;
			
			$sum_aligned += $tot * 2;
			
		}

	}
	
	my $avg = $sum / $n;
	
	say $i_out "AVG_DIFF_PROP\tAVG_DIFF_N\tAVG_ALIGNED_POSITIONS\tALL_POSITIONS\tALIGNED_PROP\tCOMPARISONS\tCOMPARISONS_PER_SAMPLE";

	say $i_out $avg , "\t" , ( $avg * $sum_aligned ) / $n , "\t" , $sum_aligned / $n , "\t" , $length_seq , "\t" , ( $sum_aligned / $n ) / $length_seq , "\t" , $n , "\t" , $n / @headers;
	
}

sub print_distances {
	
	# Print phylip distance matrix
	
	say $out scalar @headers;
	
	for ( my $i = 0; $i < @headers; $i++ ) {

		my $header = $headers[ $i ];

		print $out $header;

		for ( my $j = 0; $j < @headers; $j++ ) {
	
			my $other_header = $headers[ $j ];
			
			my $dist;
			
			if ( $header eq $other_header ) {
			
				$dist = 0;
			
			}
			
			else {
			
				$dist = $dist_table->{ $header }->{ $other_header };
			
			}
			
			printf $out "\t%.4f" , $dist;
	
		}
		
		say $out "";
		
	}
	
}

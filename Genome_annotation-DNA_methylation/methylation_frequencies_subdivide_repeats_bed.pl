#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;
use Data::Dumper;

# Get the repeats file (this can be an empty file)

my $repeats_file = shift @ARGV;

# Open the GTF file for reading

say "Reading repeat coordinates from $repeats_file ...";

open ( REPEATS_IN , "<" , $repeats_file ) or die "$!";

my $repeat_table = {};

my $sum_repeats = 0;

# Read the file

while (<REPEATS_IN>) {

	my $line = $_;

	chomp $line;
	
	next if $line =~ m/^\#/;

	# Grab info from the repeats line
	
	if ( $line =~ m/^\>(\S+)\:(\d+)\-(\d+)/ or $line =~ m/^(\S+)\s(\d+)\s(\d+)/ ) {
	
		my ( $contig , $start , $stop ) = ( $1 , $2 , $3 );
	
		unless ( defined $repeat_table->{ $contig } ) {
		
			say "Repeats in contig $contig ($start , $stop)...";
		
			# last if $contig eq "ctg1001";
		
		}
	
		push @{ $repeat_table->{ $contig } } , [ $start , $stop ];
	
		$sum_repeats += $stop - $start + 1;
	
	}
	
}

say "Read repeat intervals spanning $sum_repeats bases";

# Open the methylation frequency files for reading

foreach my $methylation_file ( @ARGV ) {

	my $methylation_in;
	
	if ( $methylation_file =~ m/\.gz$/ ) {
	
        open ( $methylation_in , "zcat $methylation_file |" ) or die "$!";
	
	}
	
	else {
	
        open ( $methylation_in , "<" , $methylation_file ) or die "$!";
	
	}
	
	say "Reading methylation frequencies from $methylation_file ...";
	
	open ( MET_REPEATS_OUT , "| gzip -f -c > ${methylation_file}.repeats.csv.gz -" ) or die "$!";
	open ( MET_NOREPEATS_OUT , "| gzip -f -c > ${methylation_file}.not_repeats.csv.gz -" ) or die "$!";
	
	# Remove the header line
	
	my $header = <$methylation_in>;
	chomp $header;
	
	say MET_REPEATS_OUT $header;
	say MET_NOREPEATS_OUT $header;
	
	# Loop over the file
	
	while (<$methylation_in>) {

		my $line = $_;
	
		chomp $line;
		
		my @data = split ( /\s+/ , $line );
		
		# chromosome      start   end     num_motifs_in_group     called_sites    called_sites_methylated methylated_frequency    group_sequence
		# ctg1    77      93      3       15      0       0.000   GCTTCCGCTTGTGTCGCTCCACGTGCC
		
		# Get information about this CpG location
		
		my (
		
			$contig , $start , $stop , $num_motifs ,
			$called_sites , $freq , $seq
			
		) = @data[ 0 , 1 , 2 , 3 , 4 , 6 , 7 ];
		
		# last if $contig eq "ctg1001";
		
		$freq = $freq * 100;
		
		# Swap coordinates around if they are in opposite order
	
		if ( $start > $stop ) {
		
			( $start , $stop ) = ( $stop , $start );
		
		}
		
		my $avg_dna_depth = $called_sites / $num_motifs;

		# say "Checking $contig position $start to $stop: $seq $freq ...";
		
		# Fix the contig name :-)
		
		$contig =~ s/\_pilon$//;
		
		my @positions;
		
		# Check if it is only a single CpG site
		
		if ( $num_motifs == 1 ) {
		
			push @positions , $start;
		
		}
		
		# If there is more than one, split them into individual sites
		
		elsif ( $num_motifs > 1 ) {
		
			# Make a temporary copy of the sequence motif
		
			my $tmp_seq = $seq;
		
			# say "$num_motifs\t$tmp_seq";
		
			# Temporarily recode the sequence
			
			$tmp_seq =~ s/CG/M/g;
		
			# Cut away sequence before the first CG (or M)
			
			$tmp_seq =~ s/^[^M]+M/M/;
		
			# Recode it back
		
			$tmp_seq =~ s/M/CG/g;
		
			# say "$num_motifs\t$tmp_seq";
		
			my $seq_length = length $tmp_seq;
			
			# Get the exact positions for every methylated CG site
			
			for ( my $i = 0; $i < $seq_length - 1; $i++ ) {
			
				my $tmp_pos = $start + $i;
				
				my $base = substr ( $tmp_seq , $i , 1 );
				my $next_base = substr ( $tmp_seq , $i + 1, 1 );
			
				if ( $base eq "C" and $next_base eq "G" ) {
				
					push @positions , $tmp_pos;
					
					# say "Adding $tmp_pos";
				
				}
			
			}
		
		}
		
		# say "CpG sites at: @positions";
		
		# Check every position (may be more than one)
		
		my $in_repeat = 0;
		my $not_in_repeat = 0;
		
		foreach my $pos ( @positions ) {
		
			# Check if there are repeats known for this contig
		
			if ( defined $repeat_table->{ $contig } ) {
			
				foreach my $interval_ref ( @{ $repeat_table->{ $contig } } ) {
				
					if ( $pos >= $interval_ref->[ 0 ] and $pos <= $interval_ref->[ 1 ] ) {
					
						say "Position $contig $pos is within $interval_ref->[ 0 ] - $interval_ref->[ 1 ]";
					
						$in_repeat = 1;
						last;
					
					}
				
				}
				
				unless ( $in_repeat ) {
				
					say "Position $contig $pos is not within any repeat";
				
					$not_in_repeat = 1;
				
				}
			
			}
			
			else {
			
				say "No repeat on contig $contig";
			
				$not_in_repeat = 1;
			
			}
		
		}
		
		if ( $in_repeat ) {
		
			say MET_REPEATS_OUT $line;
		
		}
		
		if ( $not_in_repeat ) {
		
			say MET_NOREPEATS_OUT $line;
		
		}
	
	}

}

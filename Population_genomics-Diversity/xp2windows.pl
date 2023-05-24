#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $win_size = shift @ARGV;

my $out_tag = shift @ARGV;

say "Checking files:\t@ARGV";

open ( my $norm_out , ">" , $out_tag . ".windows.${win_size}bp.csv" ) or die "$!";

say $norm_out "CHROM\tSTART\tSTOP\tN\tN_CRIT\tPROP_CRIT\tMIN\tMAX\tMEAN";

foreach my $norm_file ( @ARGV ) {

	if ( $norm_file =~ m/^[^\.]+\.([^\.]+)\./ ) {

		my $chrom = $1;
	
		open ( my $norm_in , "<" , $norm_file ) or die "$!";
		
		my $header = <$norm_in>;
		
		my $stat_list = [];
		
		while (<$norm_in>) {
		
			chomp;
			
			my @data = split ( /\t/ , $_ );
			
			my ( $pos , $norm , $crit ) = @data[ 1 , 8 , 9 ];
			
			my $pos_i = $pos / $win_size;
			
			$stat_list->[ $pos_i ]->{ 'sum' } += $norm;
			$stat_list->[ $pos_i ]->{ 'n' }++;
			
			if ( $crit != 0 ) {
			
				$stat_list->[ $pos_i ]->{ 'crit' }++;
			
			}
			
			if ( defined $stat_list->[ $pos_i ]->{ 'min' } ) {
			
				if ( $norm < $stat_list->[ $pos_i ]->{ 'min' } ) {
				
					$stat_list->[ $pos_i ]->{ 'min' } = $norm;
				
				}
			
			}
			
			else {
			
				$stat_list->[ $pos_i ]->{ 'min' } = $norm;
			
			}
			
			if ( defined $stat_list->[ $pos_i ]->{ 'max' } ) {
			
				if ( $norm > $stat_list->[ $pos_i ]->{ 'max' } ) {
				
					$stat_list->[ $pos_i ]->{ 'max' } = $norm;
				
				}
			
			}
			
			else {
			
				$stat_list->[ $pos_i ]->{ 'max' } = $norm;
			
			}
		
		}
		
		my $j = 1;
		
		for ( my $i = 0; $i < @{ $stat_list }; $i++ ) {
		
			print $norm_out "$chrom\t$j\t" , $j + $win_size;
			
			$stat_list->[ $i ]->{ 'n' } = 0 unless defined
			$stat_list->[ $i ]->{ 'n' };
			
			$stat_list->[ $i ]->{ 'crit' } = 0 unless defined
			$stat_list->[ $i ]->{ 'crit' };
			
			my $n = $stat_list->[ $i ]->{ 'n' };
			my $crit = $stat_list->[ $i ]->{ 'crit' };
			
			print $norm_out "\t$n\t$crit";
			
			if ( $stat_list->[ $i ]->{ 'n' } ) {
			
				print $norm_out "\t" , $crit / $n;
			
				print $norm_out "\t" , $stat_list->[ $i ]->{ 'min' };
				print $norm_out "\t" , $stat_list->[ $i ]->{ 'max' };
				print $norm_out "\t" , $stat_list->[ $i ]->{ 'sum' } /
				$n;
				
			}
			
			else {
			
				print $norm_out "\t\t\t\t";
			
			}
			
			say $norm_out "";
		
			$j += $win_size;
		
		}
		
	}
}

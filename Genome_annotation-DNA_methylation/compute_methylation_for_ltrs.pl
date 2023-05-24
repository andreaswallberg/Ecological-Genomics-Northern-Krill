#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $repeat_table = {};

open ( my $list_in , "<" , $ARGV[ 0 ] ) or die "$!";

while (<$list_in>) {

	chomp;
	
	my ( $repeat , $id ) = split ( /\t/ , $_ );
	
	if ( $repeat =~ m/^(.+)\_(\d+)\_(\d+)$/ ) {
	
		push @{ $repeat_table->{ $1 } } , [ $2 , $3 , $repeat , undef , undef , $id ];
	
	}

}

open ( my $methylation_in , "<" , $ARGV[ 1 ] ) or die "$!";

say "Reading methylation coordinates from " , $ARGV[ 1 ];

my $m_header = <$methylation_in>;

my $i = 1;

while (<$methylation_in>) {

	chomp;
	
	if ( $_ =~ m/^(\S+)/ ) {
	
		my $chrom = $1;
		
		if ( defined $repeat_table->{ $chrom } ) {
	
			my ( $chrom , $pos , $type , $dna , $methylation ) = split ( /\t/ , $_ );
			
			foreach my $repeat_ref ( @{ $repeat_table->{ $chrom } } ) {
			
				if ( $pos >= $repeat_ref->[ 0 ] and $pos <= $repeat_ref->[ 1 ] ) {
				
					my $name = $repeat_ref->[ 2 ];
				
					$repeat_ref->[ 3 ] += $methylation;
					$repeat_ref->[ 4 ] ++;
					
					push @{ $repeat_ref->[ 6 ] } , $methylation;
				
					say "Methylation point $i ...";
				
					last;
				
				}
			
			}
	
		}

	}
	
	$i++;
	
}

my $dist_table = [];

open ( my $out , ">" , $ARGV[ 0 ] . ".95.tsv" ) or die "$!";

say $out "CHROM\tREPEAT\tREPEAT_START\tREPEAT_STOP\tREPEAT_LENGTH\tIDENTITY\tMETHYLATION_N\tMETHYLATION";

foreach my $chrom ( sort keys %{ $repeat_table } ) {

	foreach my $repeat_ref ( @{ $repeat_table->{ $chrom } } ) {

		my $methylation_n = "";
		my $methylation = "";

		if ( defined $repeat_ref->[ 4 ] ) {
			
			$methylation_n = $repeat_ref->[ 4 ];
			$methylation = $repeat_ref->[ 3 ] / $methylation_n;
		
			my $id = $repeat_ref->[ 5 ];
		
			say $out $chrom , "\t" , $repeat_ref->[ 2 ] , "\t" ,
			$repeat_ref->[ 0 ] , "\t" , $repeat_ref->[ 1 ] , "\t" ,
			$repeat_ref->[ 1 ] - $repeat_ref->[ 0 ] + 1 , "\t" ,
			$id , "\t" ,
			$methylation_n , "\t" , $methylation;
			
			my $id_cat = $id / 5;
			
			# Per Cpg site
			
			$dist_table->[ $id_cat ]->{ 'n_cpgs' } += $methylation_n;
			$dist_table->[ $id_cat ]->{ 'n_cpgs_methylation' } += $repeat_ref->[ 3 ];
			push @{ $dist_table->[ $id_cat ]->{ 'n_cpgs_methylations' } } , @{ $repeat_ref->[ 6 ] };
			
			# Per repeat
			
			$dist_table->[ $id_cat ]->{ 'n' }++;
			$dist_table->[ $id_cat ]->{ 'n_methylation' } += $methylation;
			push @{ $dist_table->[ $id_cat ]->{ 'n_methylations' } } , $methylation;
			
		}
		
	}
	
}

open ( my $d_out , ">" , $ARGV[ 0 ] . ".dist.95.tsv" ) or die "$!";

say $d_out "IDENTITY\tN_REPEATS\tAVG_METHYLATION_N_REPEATS\tLOW\tHIGH\tN_CPGS\tAVG_METHYLATION_N_CPGS\tLOW\tHIGH";

for ( my $i = 0; $i <= 20; $i++ ) {

	my $div = $i * 5;
	
	my $n_repeats = 0;
	my $n_repeats_methylation = "";
	
	my $n_cpgs = 0;
	my $n_cpgs_methylation = "";
	
	my $lower = "";
	my $upper = "";
	
	my $cpg_lower = "";
	my $cpg_upper = "";
	
	if ( defined $dist_table->[ $i ]->{ 'n' } ) {
	
		$n_repeats = $dist_table->[ $i ]->{ 'n' };
		$n_repeats_methylation = $dist_table->[ $i ]->{ 'n_methylation' } / $n_repeats;
		
		say "Bootstrapping for divergence interval $div ...";
		
		if ( defined $dist_table->[ $i ]->{ 'n_methylations' } ) {
		
			( $lower , $upper ) = bootstrap( $dist_table->[ $i ]->{ 'n_methylations' } );
		
			$lower = $n_repeats_methylation - $lower;
			$upper = $upper - $n_repeats_methylation;
		
		}
		
		$n_cpgs = $dist_table->[ $i ]->{ 'n_cpgs' };
		$n_cpgs_methylation = $dist_table->[ $i ]->{ 'n_cpgs_methylation' } / $n_cpgs;
		
		if ( defined $dist_table->[ $i ]->{ 'n_cpgs_methylations' } ) {
		
			( $cpg_lower , $cpg_upper ) = bootstrap( $dist_table->[ $i ]->{ 'n_cpgs_methylations' } );
		
			$cpg_lower = $n_cpgs_methylation - $cpg_lower;
			$cpg_upper = $cpg_upper - $n_cpgs_methylation;
		
		}
	
	}
	
	say $d_out "$div\t$n_repeats\t$n_repeats_methylation\t$lower\t$upper\t$n_cpgs\t$n_cpgs_methylation\t$cpg_lower\t$cpg_upper";

}

sub bootstrap {

	my ( $vals_ref ) = @_;
	
	my @vals = @{ $vals_ref };
		
	if ( @vals > 1 ) {
		
		my $n_vals = @vals;
		
		my @reps;
		
		my $n_reps = 200;
		
		for ( my $i = 0; $i <= $n_reps; $i++ ) {
		
			my $sum;
			
			for ( my $j = 0; $j <= $n_vals; $j++ ) {
			
				$sum += $vals[ int ( rand ( $n_vals ) ) ];
			
			}
			
			push @reps , $sum / $n_vals;
		
		}
		
		my @sorted_reps = sort { $a <=> $b } @reps;
		
		my $lower = $sorted_reps[ $n_reps * 0.025 ];
		my $upper = $sorted_reps[ $n_reps * 0.975 ];
		
		return ( $lower , $upper );
	
	}
	
	else {
	
		return ( $vals[ 0 ] , $vals[ 0 ] );
	
	}

}

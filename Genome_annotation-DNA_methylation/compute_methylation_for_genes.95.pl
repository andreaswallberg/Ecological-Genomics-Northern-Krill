#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $keep_table = {};

open ( my $keep_in , "<" , $ARGV[ 0 ] ) or die "$!";

my $header = <$keep_in>;

while (<$keep_in>) {

	chomp;
	
	if ( $_ =~ m/^(\S+)\s(\S+)\s(\S+)/ ) {
	
		my ( $gene , $source , $type ) = ( $1 , $2 , $3 );
		
		if ( $type eq "GENE" ) {
		
			$keep_table->{ $gene } = 0;
		
		}
	
	}

}

open ( my $isoforms_in , "<" , $ARGV[ 1 ] ) or die "$!";

while (<$isoforms_in>) {

	chomp;
	
	my ( $gene , $isoforms ) = split ( /\t/ , $_ );
	
	if ( defined $keep_table->{ $gene } ) {
	
		$keep_table->{ $gene } = $isoforms;
	
	}

}

open ( my $genes_in , "<" , $ARGV[ 2 ] ) or die "$!";

my $genes_header = <$genes_in>;
chomp $genes_header;

my @genes_headers = split ( /\t/ , $genes_header );

my $methylation_table = [];
my $isoform_table = [];

my $region_table = {};

while (<$genes_in>) {

	chomp;
	
	my @data = split ( /\t/ , $_ );
	
	my $gene = $data[ 1 ];
	
	if ( defined $keep_table->{ $gene } and $keep_table->{ $gene } ) {
	
		# CHROM	GENE	3_UTR_N	3_UTR_FREQ	5_UTR_N	5_UTR_FREQ	ALL_N	ALL_FREQ	cds_N	cds_FREQ	exon_N	exon_FREQ	genic_N	genic_FREQ	intergenic_N	intergenic_FREQ	intron_N	intron_FREQ

		for ( my $i = 2; $i < @data; $i += 2 ) {
		
			my $region = $genes_headers[ $i ];
			
			$region =~ s/\_N//;
		
			$region_table->{ $region } = 1;
		
			my $n = $data[ $i ];
			my $freq = $data[ $i + 1 ];
			
			next unless defined $n;
			next unless defined $freq;
			
			if ( $n =~ m/\d/ and $freq =~ m/\d/ ) {
			
				$freq = 99 if $freq >= 100;
			
				$freq /= 10;
			
				my $isoforms = $keep_table->{ $gene };
				
				push @{ $methylation_table->[ $freq ]->{ $region }->{ 'isoforms' } } , $isoforms;
				
				$methylation_table->[ $freq ]->{ $region }->{ 'isoform' } += $isoforms;
				
				$methylation_table->[ $freq ]->{ $region }->{ 'n_genes' }++;
			
			}
			
		}
	
	}

}

my @sorted_regions = sort keys %{ $region_table };

open ( my $out , ">" , $ARGV[ 2 ] . ".isoforms_over_methylation.95.tsv" ) or die "$!";

say $out "METHYLATION" , map { "\t" . "${_}_N" . "\t" . "${_}_ISOFORMS" . "\t" . "${_}_LOW" . "\t" . "${_}_HIGH" } @sorted_regions;

for ( my $i = 0; $i <= 10; $i++ ) {

	print $out $i * 10;

	foreach my $region ( @sorted_regions ) {
	
		my $avg_isoforms = "";
		
		if ( defined $methylation_table->[ $i ]->{ $region }->{ 'isoform' } ) {
		
			my $n = $methylation_table->[ $i ]->{ $region }->{ 'n_genes' };
		
			my $mean = $methylation_table->[ $i ]->{ $region }->{ 'isoform' } / $n;
		
			my ( $low , $high ) = bootstrap( $methylation_table->[ $i ]->{ $region }->{ 'isoforms' } );
		
			print $out "\t$n\t" , $mean , "\t" , $mean - $low , "\t" , $high - $mean;
		
			if ( lc $region eq "exon" ) {
			
				open ( my $e_out , ">" , $ARGV[ 2 ] . ".isoforms_over_methylation.${region}.${i}.95.tsv" ) or die "$!";
			
				say $e_out $_ foreach @{ $methylation_table->[ $i ]->{ $region }->{ 'isoforms' } };
			
			}
		
		}
	
	}
	
	say $out "";

}

sub bootstrap {

	my ( $vals_ref ) = @_;
	
	my @vals = @{ $vals_ref };
		
	if ( @vals > 1 ) {
		
		my $n_vals = @vals;
		
		my @reps;
		
		my $n_reps = 1000;
		
		for ( my $i = 0; $i <= $n_reps; $i++ ) {
		
			my $sum;
			
			my $n = 0;
			
			for ( my $j = 0; $j < $n_vals; $j++ ) {
			
                my $val = $vals[ int ( rand ( $n_vals ) ) ];
			
				$sum += $val;
				
				$n++;
			
			}
			
			push @reps , $sum / $n_vals;
		
		}
		
		my @sorted_reps = sort { $a <=> $b } @reps;
		
		# say STDERR "@sorted_reps";
		
		my $lower = $sorted_reps[ $n_reps * 0.025 ];
		my $upper = $sorted_reps[ $n_reps * 0.975 ];
		
		return ( $lower , $upper );
	
	}
	
	else {
	
		return ( $vals[ 0 ] , $vals[ 0 ] );
	
	}
	
}

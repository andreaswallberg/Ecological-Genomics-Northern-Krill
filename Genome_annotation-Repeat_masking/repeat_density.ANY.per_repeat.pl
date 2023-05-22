#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;
use Data::Dumper;

my $region_table = {};
my $count_table = {};
my $repeat_table = {};

my $cluster = shift @ARGV;
my $genome_size = shift @ARGV;

my $region_length_tot = $genome_size;
my $repeat_length_tot = 0;

my $tag = shift @ARGV;
my $in_file = shift @ARGV;

$count_table->{ 'ANY' }->{ 'length' } = $region_length_tot;

foreach my $repeat_file ( @ARGV ) {

	open ( my $in , "<" , $repeat_file ) or die "$!";

	say "Reading repeat file $repeat_file ...";
	
	while (<$in>) {
	
		chomp;
		
		$_ =~ s/^\s+//g;
		$_ =~ s/\s+$//g;
		
		next if $_ =~ m/^\s*\#/;
		next unless $_ =~ m/^\d/;
		
		my @data;
		
		if ( $_ =~ m/\t/ ) {
		
			@data = split ( /\t/ , $_ );
		
		}
		
		else {
		
			@data = split ( /\s+/ , $_ );
		
		}
		
		my ( $repeat , $chrom , $start , $stop );
		
		# BLAST input = 12
		# AluI	GroupUN_108	95.48	177	0	3	14	189	2982	2813	3e-73	276
		
		if ( @data == 12 ) {
		
			( $repeat , $chrom , $start , $stop ) =
			@data[ 0 , 1 , 8 , 9 ];
		
		}
		
		# Repeatmasker input = 14
		# 258 19.41 0.00 0.75            Group11_2     34564     34698 (8450413) +              (TAAA)n   Simple_repeat      44     177     (3)

		if ( @data >= 14 ) {
		
			( $repeat , $chrom , $start , $stop ) =
			@data[ 10 , 4 , 5 , 6 ];
		
			$repeat .= "-" . $data[ 9 ]; 
		
		}
	
		else {
		
			die "Ooops: $_";
		
		}
		
		if ( defined $repeat ) {

			if ( $cluster eq "narrow" ) {

				# Cluster repeats narrowly by DNA/TcMar, DNA/Kolobok and orders etc
			
				# $repeat =~ s/\-.*//g;
			
			}
			
			if ( $cluster eq "superfamily" ) {

				# Cluster repeats narrowly by DNA/TcMar, DNA/Kolobok etc
			
				$repeat =~ s/\-.*//g;
			
			}
			
			else {
		
				# Cluster repeats broadly by DNA, LINE, LTR, SINE etc
		
				$repeat =~ s/\/.*//g;
			
			}
			
			$repeat =~ s/\?//g;
		
# 			if ( $repeat eq "Low_complexity" ) {
# 			
# 				$repeat = "Simple_repeat";
# 			
# 			}
# 			
# 			elsif ( $repeat =~ m/SINE/ ) {
# 			
# 				$repeat = "SINE";
# 			
# 			}
# 			
# 			elsif ( $repeat =~ m/RNA/ ) {
# 			
# 				$repeat = "RNA";
# 			
# 			}
		
			( $start , $stop ) = sort { $a <=> $b } ( $start , $stop );
		
			my $repeat_length = $stop - $start + 1;
		
			$repeat_table->{ $repeat }->{ 'length' } += $repeat_length;
			$repeat_table->{ $repeat }->{ 'n' }++;
			
			$repeat_length_tot += $repeat_length;
		
			foreach my $type ( 'ANY' ) {

				my $overlap_length =
				$stop - $start + 1;
						
				$count_table->{ $type }->{ 'repeat' }->
				{ $repeat }->{ 'n' }++;
				
				$count_table->{ $type }->{ 'repeat' }->
				{ $repeat }->{ 'length' } +=
				$overlap_length;
			
				$repeat_table->{ $repeat }->{ 'type' }->
				{ $type }->{ 'n' }++;
				
				$repeat_table->{ $repeat }->{ 'type' }->
				{ $type }->{ 'length' } +=
				$overlap_length;
			
			}
			
			say $repeat_length_tot;
			
			# last if $repeat_length_tot >= 1_000_000_000;
		
		}
	
	}
	
}

if ( $region_length_tot and $repeat_length_tot ) {

	open ( my $out , ">" , "${in_file}.repeats_in_ANY_region.${tag}.per_repeat.csv" ) or die "$!";

	say "${in_file}.repeats_in_ANY_region.${tag}.csv";
	
	say $out
	"TYPE\tTYPE_REPEAT\tREPEAT\tN\tREPEAT_LENGTH\tREPEAT_LENGTH_TOT\tPROP\tALL_REPEAT_LENGTH\t" ,
	"REGION_LENGTH\tGENOME_LENGTH\tLENGTH_PROP\tREPEAT_PROP\tBIAS";

	my $sum_repeat_table = {};
	
	foreach my $type ( keys %{ $count_table } ) {
		
		foreach my $repeat (
		
			sort keys %{ $count_table->{ $type }->{ 'repeat' } }
		
		) {
		
			my $repeat_n =
			$count_table->{ $type }->{ 'repeat' }->{ $repeat }->{ 'n' };
		
			$repeat_n = 0 unless defined $repeat_n;
		
			my $repeat_length_in_type =
			$count_table->{ $type }->{ 'repeat' }->{ $repeat }->{ 'length' };
		
			$repeat_length_in_type = 0 unless defined $repeat_length_in_type;
		
			my $repeat_length = $repeat_table->{ $repeat }->{ 'length' };
		
			my $region_length = $count_table->{ $type }->{ 'length' };
			
			my $repeat_prop_in_type = $repeat_length_in_type / $repeat_length;
			
			my $length_prop = $region_length / $region_length_tot;
			
			my $bias = $repeat_prop_in_type / $length_prop;
			
			say $out "$type\t${type}:${repeat}\t$repeat\t$repeat_n\t$repeat_length_in_type\t$repeat_length\t" ,
			$repeat_prop_in_type , "\t" ,
			$repeat_length_tot , "\t" ,
			$region_length , "\t" , $region_length_tot , "\t" ,
			$length_prop , "\t" ,
			$repeat_prop_in_type , "\t" ,
			$bias;
			
			$sum_repeat_table->{ $repeat }->{ 'repeat_length' } = $repeat_length;
			$sum_repeat_table->{ $repeat }->{ 'repeat_n' } += $repeat_n;

			$sum_repeat_table->{ $repeat }->{ 'type' }->{ $type } = {
			
				'prop' => $repeat_prop_in_type ,
				'bias' => $bias ,
				
			};

		
		}
		
	}
	
	open ( my $sum_out , ">" , "${in_file}.sum_repeats_in_ANY_regions.${tag}.per_repeat.csv" ) or die "$!";

	# Placed vs unplaced
	
	my @types = ( 'ANY' );
	
	say $sum_out "REPEAT\tSIZE\tN\tINFO\tLABEL\t" ,
	join ( "\t" , @types ) , "\t" ,
	join ( "\t" , map { $_ . "_BIAS" } @types );
	
	foreach my $repeat ( sort keys %{ $sum_repeat_table } ) {

		my $tot_length = $repeat_table->{ $repeat }->{ 'length' };
		my $tot_n = $repeat_table->{ $repeat }->{ 'n' };
	
		print $sum_out "${repeat}\t${tot_length}\t${tot_n}";
	
		my $print_length;
		
		if ( $tot_length >= 1_000_000 ) {
		
			$print_length = int ( ( $tot_length + 500_000 ) / 1_000_000 );
		
			$print_length .= "Mbp";
		
		}
		
		elsif ( $tot_length >= 1_000 ) {
		
			$print_length = int ( ( $tot_length + 500 ) / 1_000 );
		
			$print_length .= "kbp";
		
		}
		
		else {
		
			$print_length = $tot_length . "bp";
		
		}
		
		$print_length .= "; n=$tot_n";
		
		print $sum_out "\t$print_length\t$repeat ($print_length)";
		
		foreach my $type ( @types ) {
		
			my $prop = $sum_repeat_table->{ $repeat }->
			{ 'type' }->{ $type }->{ 'prop' };
			
			$prop = 0 unless defined $prop;
		
			print $sum_out "\t" , $prop * 100;
		
		}
		
		foreach my $type ( @types ) {

			my $bias = $sum_repeat_table->{ $repeat }->
			{ 'type' }->{ $type }->{ 'bias' };
			
			$bias = 0 unless defined $bias;
		
			print $sum_out "\t" , $bias;
		
		}
		
		say $sum_out "";
	
	}
	
}

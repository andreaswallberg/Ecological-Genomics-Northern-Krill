#!/usr/bin/perl -w

# Enforce some programming rules/syntax

use strict;
use warnings;
use 5.010;
use Data::Dumper;
use Getopt::Long;

# Option table and defaults

my $opt = {

	'group_table' => {} ,

	'coverage' => [] ,
	'hetero' => [] ,
	
	'methylation' => [] ,
	
};

GetOptions(
    $opt ,

    # Input
    # -----

    # Minimum coverage
    'dna=i' ,
    
    # Region sequence and repeat coordinates
    'regions=s' ,
    'repeats=s' ,
    
    # Determines how to treat repeats

    # 0 = skip CpGs in repeats
    # 1 = only consider CpGs in repeats
    # 2 = do not distinguish between CpGs in repeats or non-repeat sequences (i.e. look at all CpGs)
    
    'repeat_location=i' ,
    
    # Coverage files
    
    'coverage=s{,}' ,
    
    # Heterozygous CpG locations to ignore
    
    'hetero=s{,}' ,
    
    # Methylation file(s)
    
    'methylation=s{,}' ,
    
    # Output
    # ------
    
    # Output name (string)
    'output|out|o=s' ,
    
    'verbose|v' ,
    
    'no_print_snps|i' ,

);

my $min_dna_coverage = $opt->{ 'dna' };

# Get the stem base name/paths of the output files

my $out_base = $opt->{ 'output' };

my $gene_table = {};
my $repeat_table = {};
my $region_table = {};
my $coverage_table = {};
my $hetero_table = {};
my $methylation_table = {};

my @contigs;
my $max_region_i = 0;

# Per-CpG site output

open ( my $s_out , ">" , "${out_base}.cpg_sites.csv" ) or die "$!";

say $s_out "CONTIG\tPOS\tREGION\tDNA\tMETHYLATION_FREQ";

read_regions();
read_repeats();
read_coverage();
read_hetero();
read_methylation();

sub mean {

	my ( $data_ref ) = @_;
	
	my $mean = 0;
	
	map { $mean += $_ } @{ $data_ref };
	
	$mean /= @{ $data_ref };

	return $mean;
	
}

sub median {

	my ( $data_ref ) = @_;
	
	my $n = @{ $data_ref };
	
	my @sorted_data = sort { $a <=> $b } @{ $data_ref };
	
	my $n_freqs = @sorted_data;
	
	my $median;
	
	if ( $n_freqs % 2 ) {
	
		$median = $sorted_data[ int ( $n_freqs / 2 ) ];
	} 
	
	# If even
	
	else {
	
		$median = (
		
			$sorted_data[ ( int ( $n_freqs / 2 ) ) - 1 ] +
			$sorted_data[ int ( $n_freqs / 2 ) ]
			
		) / 2;
	
	}
	
	return $median;

}

sub read_regions {

	# Open the repeat file for reading

	my $region_file = $opt->{ 'regions' };

	open ( my $region_in , "<" , $region_file ) or die "$!";

	my $tmp_header;
	
	while (<$region_in>) {

		chomp;
		
		if ( $_ =~ m/^\>(.+)/ ) {
			
			$tmp_header = $1;
		
		}
		
		elsif ( defined $tmp_header ) {
		
			$region_table->{ $tmp_header } .= $_;
		
		}

	}

}

sub read_repeats {

	# Open the repeat file for reading

	my $repeat_file = $opt->{ 'repeats' };

	open ( my $repeat_in , "<" , $repeat_file ) or die "$!";

	my $tmp_header;
	
	while (<$repeat_in>) {

		chomp;
		
		if ( $_ =~ m/^\>(.+)/ ) {
			
			$tmp_header = $1;
		
		}
		
		elsif ( defined $tmp_header ) {
		
			$repeat_table->{ $tmp_header } .= $_;
		
		}

	}

}

sub read_coverage {

	# Open the coverage files for reading

	foreach my $coverage_file ( @{ $opt->{ 'coverage' } } ) {
	
		my $coverage_in;
		
		if ( $coverage_file =~ m/\.gz$/ ) {
		
			open ( $coverage_in , "zcat $coverage_file |" ) or die "$!";
		
		}
		
		else {
		
			open ( $coverage_in , "<" , $coverage_file ) or die "$!";
		
		}
		
		say "Reading coverage from $coverage_file ...";
		
		my $tmp_header;
		
		while (<$coverage_in>) {
		
			chomp;
			
			if ( $_ =~ m/^\>(.+)/ ) {
			
				$tmp_header = $1;
			
			}
			
			else {
			
				$coverage_table->{ $tmp_header } .= $_;
			
			}
		
		}
	
	}
	
}

sub read_hetero {

	# Open the heterozygous genotype files for reading

	foreach my $hetero_file ( @{ $opt->{ 'hetero' } } ) {
	
		my $hetero_in;
		
		if ( $hetero_file =~ m/\.gz$/ ) {
		
			open ( $hetero_in , "zcat $hetero_file |" ) or die "$!";
		
		}
		
		else {
		
			open ( $hetero_in , "<" , $hetero_file ) or die "$!";
		
		}
		
		say "Reading heterozygous genotypes from $hetero_file ...";
		
		while (<$hetero_in>) {
		
			chomp;
			
			if ( $_ =~ m/^(\S+)\s(\d+)/ ) {
			
				substr ( $coverage_table->{ $1 } , $2 - 1 , 1 , "0" );
			
			}
		
		}
	
	}
	
}

sub read_methylation {

	my %seen_contig;

	my %seen_region;

	# Check for repeat or not

	my $check_in_repeat = $opt->{ 'repeat_location' };

	my $region_code = {
		'1' => 'intergenic' ,
		'2' => 'intron' ,
		'3' => '3_UTR' ,
		'4' => 'exon' ,
		'5' => '5_UTR' ,
		'6' => 'cds' ,
	};
	
	# Open the methylation frequency files for reading

	foreach my $methylation_file ( @{ $opt->{ 'methylation' } } ) {

		my $methylation_in;
		
		if ( $methylation_file =~ m/\.gz$/ ) {
		
			open ( $methylation_in , "zcat $methylation_file |" ) or die "$!";
		
		}
		
		else {
		
			open ( $methylation_in , "<" , $methylation_file ) or die "$!";
		
		}
		
		say "Reading methylation frequencies from $methylation_file ...";
		
		# Discard the header line
		
		my $header = <$methylation_in>;
		chomp $header;
		
		# Loop over the file
		
		while (<$methylation_in>) {
		
			chomp;
			
			my @data = split ( /\s+/ , $_ );
			
			# chromosome      start   end     num_motifs_in_group     called_sites    called_sites_methylated methylated_frequency    group_sequence
			# ctg1    77      93      3       15      0       0.000   GCTTCCGCTTGTGTCGCTCCACGTGCC
			
			my ( $contig , $start , $stop , $num_motifs , $called_sites , $freq , $seq ) =
			@data[ 0 , 1 , 2 , 3 , 4 , 6 , 7 ];
			
			$freq *= 100;
			
			my $avg_dna_depth = $called_sites / $num_motifs;
			
			next if $avg_dna_depth < $min_dna_coverage;
			
			unless ( defined $seen_contig{ $contig } ) {
			
				$seen_contig{ $contig } = 1;
				
				say "On contig $contig n=" , scalar @contigs;
				
				push @contigs , $contig;
			
			}
			
			# Temporarily recode the sequence
			
			my $tmp_seq = $seq;
			
			$tmp_seq =~ s/CG/M/g;
		
			# Cut away sequence before the first CG (or M)
			
			$tmp_seq =~ s/^[^M]+M/M/;
		
			# Recode it back
		
			$tmp_seq =~ s/M/CG/g;
		
			my $seq_length = length $tmp_seq;
			
			# Get the exact positions for every methylated CG site
			
			my @positions;
			
			for ( my $i = 0; $i < $seq_length - 1; $i++ ) {
			
				my $tmp_pos = $start + $i;
				
				my $base = substr ( $tmp_seq , $i , 1 );
				my $next_base = substr ( $tmp_seq , $i + 1, 1 );
			
				if ( $base eq "C" and $next_base eq "G" ) {
				
					# Check if the position should be concluded or not
					# i.e. has not been flagged for low coverage or heterozygous
					# genotypes at this location
					
					if ( defined $coverage_table->{ $contig } ) {
					
						my $symbol = substr ( $coverage_table->{ $contig } , $tmp_pos , 1 );
						
						if ( $symbol ) {
						
							push @positions , $tmp_pos;
						
						}
						
						else {
						
							say "CpG site at $contig $tmp_pos discarded ...";
						
						}
					
					}
					
					else {
				
						push @positions , $tmp_pos;
				
					}
				
				}
			
			}
			
			# Check every position if it occurs in a repeat
			
			my @valid_positions;
			
			if ( $check_in_repeat == 1 ) {
			
				if ( defined $repeat_table->{ $contig } ) {
			
					my $seq = $repeat_table->{ $contig };
				
					foreach my $pos ( @positions ) {
				
						my $repeat_base = substr ( $repeat_table->{ $contig } , $pos - 1 , 1 );
						
						if ( $repeat_base ) {

							push @valid_positions , $pos;
							
						}
						
					}
				
				}
				
			}
			
			if ( $check_in_repeat == 0 ) {
			
				if ( defined $repeat_table->{ $contig } ) {
			
					my $seq = $repeat_table->{ $contig };
				
					foreach my $pos ( @positions ) {
				
						my $repeat_base = substr ( $repeat_table->{ $contig } , $pos - 1 , 1 );
						
						unless ( $repeat_base ) {

							push @valid_positions , $pos;
							
						}
						
					}
				
				}
				
				else {
				
					foreach my $pos ( @positions ) {
				
						push @valid_positions , $pos;
				
					}
				
				}
				
			}
			
			if ( $check_in_repeat == 2 ) {

				foreach my $pos ( @positions ) {
			
					push @valid_positions , $pos;
				
				}
			
			}
			
			foreach my $pos ( @valid_positions ) {
				
				if ( defined $contig and defined $pos and defined $avg_dna_depth and defined $freq ) {
			
					my $base = substr ( $region_table->{ $contig } , $pos - 1 , 1 );
			
					my $region = $region_code->{ $base };
			
					say $s_out "$contig\t$pos\t$region\t$avg_dna_depth\t$freq";
			
				}
			
			}
		
		}

	}

}

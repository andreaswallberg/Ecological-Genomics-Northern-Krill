#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

use Getopt::Long;
use Data::Dumper;

# Option table and defaults

my $opt = {
    
    # Input files

    'stats' => [] ,
    
    # Fields to parse
    
    'fields' => [] ,
   
   'output' => "distance.csv" ,
   
};

GetOptions(
    $opt ,
    
    'genes|gtf|gff=s' ,
	
	'stats|stat|s=s{,}' ,
	'fields|field|f=s{,}' ,
	
	'window=s' ,
        
    'output|out|o=s' ,
   
    # Booleans
    
    'verbose|v' ,
    'help|h' ,

);

# Read genes
# ----------

init_genes( $opt );
read_stats( $opt );
print_output( $opt );

sub init_genes {

    my ( $opt ) = @_;
    
    my $gene_table = {};
    
    $opt->{ 'app' }->{ 'gene_table' } = $gene_table;
    
    my $window = $opt->{ 'window' };
    
    if ( defined $opt->{ 'genes' } and -e $opt->{ 'genes' } ) {
        
        my $gtf_file = $opt->{ 'genes' };
        
        say "Reading gene coordinates from $gtf_file ...";

        open ( my $gtf_in , "<" , $gtf_file ) or die "$!";
        
        while (<$gtf_in>) {
        
			chomp;

			next if $_ =~ m/^\#/;
			next unless $_ =~ m/\sgene\s/;

            my @data = split ( /\t/ , $_ );
            
            if ( $data[ 2 ] eq "gene" ) {
            
                my ( $chrom , $start , $stop , $orientation ) =
                @data[ 0 , 3 , 4 , 6 ];
                
				push @{ $gene_table->{ $chrom } } , [
					( int ( $start / $window ) ),
					( int ( $stop / $window ) ) ,
					$start ,
					$stop ,
					$orientation ,
				];
            
            }
			
        }
        
    }

}

sub read_stats {

    my ( $opt ) = @_;
    
    my $gene_table = $opt->{ 'app' }->{ 'gene_table' };
    
    my $stat_table = {};
    
    $opt->{ 'app' }->{ 'stat_table' } = $stat_table;
    
    my $window = $opt->{ 'window' };
    
    my @fields = @{ $opt->{ 'fields' } };
    
    foreach my $stat_file ( @{ $opt->{ 'stats' } } ) {
    
		open ( my $in , "<" , $stat_file ) or die "$!";
    
		say "Reading $window bp window stats from $stat_file ...";
    
		my $header_string = <$in>;
		chomp $header_string;
		
		my @headers = split ( /\t/ , $header_string );
		
		my @selected_headers = @headers[ @fields ];
        
        $opt->{ 'headers' } = \@selected_headers;
        
        while (<$in>) {
        
			chomp;

			next if $_ =~ m/^\#/;

			if ( $_ =~ m/^(\S+)\s(\d+)\s/ ) {
			
				my ( $chrom , $start ) = ( $1 , $2 );
			
				# last if $chrom eq "seq_s_1000";
			
				if ( defined $gene_table->{ $chrom } ) {
			
					# say "$chrom\t$start\t@selected_data";
            
					my $mid_start = $start + $window / 2;
					my $stop = $start + $window;
            
					my $i = 0;
            
					my $min_distance;
					my $min_gene_i;
					my $relation;
            
					foreach my $gene ( @{ $gene_table->{ $chrom } } ) {
            
						my $gene_start = $gene->[ 2 ];
						my $gene_stop = $gene->[ 3 ];
						
						# Window starts within the gene
						
						if ( $start >= $gene_start and $start <= $gene_stop ) {
						
							$min_distance = undef;
							$min_gene_i = undef;
							
							last;
						
						}
						
						# Window stops within the gene
						
						elsif ( $stop >= $gene_start and $stop <= $gene_stop ) {
						
							$min_distance = undef;
							$min_gene_i = undef;
							
							last;
						
						}
						
						# Window spans the gene
						
						elsif ( $start <= $gene_start and $stop >= $gene_stop ) {
						
							$min_distance = undef;
							$min_gene_i = undef;
							
							last;
						
						}
						
						elsif ( $start + $window < $gene_start ) {
						
							my $tmp_distance = $gene_start - $mid_start;
						
							if ( defined $min_distance ) {
							
								if ( $tmp_distance < $min_distance ) {
								
									$min_distance = $tmp_distance;
									$min_gene_i = $i;
									$relation = "upstream";
								
								}
							
							}
							
							else {
							
								$min_distance = $tmp_distance;
								$min_gene_i = $i;
								$relation = "upstream";
								
							}
						
						}
						
						elsif ( $start > $gene_stop ) {
						
							my $tmp_distance = $mid_start - $gene_stop;
						
							if ( defined $min_distance ) {
							
								if ( $tmp_distance < $min_distance ) {
								
									$min_distance = $tmp_distance;
									$min_gene_i = $i;
									$relation = "downstream";
									
								}
							
							}
							
							else {
							
								$min_distance = $tmp_distance;
								$min_gene_i = $i;
								$relation = "downstream";
							
							}
						
						}
						
						$i++;
					
					}
					
					if ( defined $min_distance ) {
					
						my @data = split ( /\t/ , $_ );
						
						my @selected_data = @data[ @fields ];
						
						my $min_distance_i = int ( $min_distance / $window );
						
						my $orientation = $gene_table->{ $chrom }->[ $min_gene_i ]->[ 4 ];
						
						if ( $orientation eq "-" ) {
						
							if ( $relation eq "upstream" ) {
							
								$relation = "downstream";
							
							}
							
							elsif ( $relation eq "downstream" ) {
							
								$relation = "upstream";
							
							}
						
						}
						
						for ( my $i = 0; $i < @selected_data; $i++ ) {
						
							my $header = $selected_headers[ $i ];
							my $field_data = $selected_data[ $i ];
						
							if ( defined $field_data ) {
						
								push @{ $stat_table->{ $relation }->[ $min_distance_i ]->{ $header }->{ 'vals' } } , $field_data;
						
								$stat_table->{ $relation }->[ $min_distance_i ]->{ $header }->{ 'sum' } += $field_data;
								$stat_table->{ $relation }->[ $min_distance_i ]->{ $header }->{ 'n' }++;
						
							}
						
						}
					
					}
            
				}
            
            }
			
        }
        
    }

}

sub print_output {

    my ( $opt ) = @_;
    
    my $gene_table = $opt->{ 'app' }->{ 'gene_table' };
    
    my $stat_table = $opt->{ 'app' }->{ 'stat_table' };
    
    my $window = $opt->{ 'window' };
        
    my @headers = @{ $opt->{ 'headers' } };
    
    open ( my $out , ">" , $opt->{ 'output' } ) or die "$!";

    print $out "DISTANCE";
    
    foreach my $header ( @headers ) {
    
		print $out "\t${header}_N\t${header}_LOW\t${header}_MEAN\t${header}_HIGH";
    
    }
    
    say $out "";
    
	my $max_upstream_distance = @{ $stat_table->{ 'upstream' } };
	
	for ( my $i = $max_upstream_distance; $i >= 1; $i-- ) {
	
		print $out "-" , $i * $window + $window / 2;
		
		foreach my $header ( @headers ) {
		
			if ( defined $stat_table->{ "upstream" }->[ $i ]->{ $header }->{ 'n' } ) {
		
				my $sum =
				$stat_table->{ "upstream" }->[ $i ]->{ $header }->{ 'sum' };
				
				my $n =
				$stat_table->{ "upstream" }->[ $i ]->{ $header }->{ 'n' };
		
				my $vals_ref =
				$stat_table->{ "upstream" }->[ $i ]->{ $header }->{ 'vals' };
		
				say "Bootstrapping upstream $i $header n=" , scalar @{ $vals_ref };
		
				my ( $lower , $upper ) = bootstrap_vals( $vals_ref );
		
				my $mean = $sum / $n;
				
				print $out "\t$n\t$lower\t$mean\t$upper";
		
			}
			
			else {
			
				print $out "\t\t";
			
			}
		
		}
		
		say $out "";
	
	}
	
	my $max_downstream_distance = @{ $stat_table->{ 'downstream' } };
	
	for ( my $i = 1; $i < $max_downstream_distance; $i++ ) {
	
		print $out $i * $window + $window / 2;
		
		foreach my $header ( @headers ) {
		
			if ( defined $stat_table->{ "downstream" }->[ $i ]->{ $header }->{ 'n' } ) {
		
				my $sum =
				$stat_table->{ "downstream" }->[ $i ]->{ $header }->{ 'sum' };
				
				my $n =
				$stat_table->{ "downstream" }->[ $i ]->{ $header }->{ 'n' };
		
				my $vals_ref =
				$stat_table->{ "downstream" }->[ $i ]->{ $header }->{ 'vals' };
		
				say "Bootstrapping downstream $i $header n=" , scalar @{ $vals_ref };
		
				my ( $lower , $upper ) = bootstrap_vals( $vals_ref );
		
				my $mean = $sum / $n;
				
				print $out "\t$n\t$lower\t$mean\t$upper";
		
			}
			
			else {
			
				print $out "\t\t";
			
			}
		
		}
		
		say $out "";
	
	}
    
}

sub bootstrap_vals {

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

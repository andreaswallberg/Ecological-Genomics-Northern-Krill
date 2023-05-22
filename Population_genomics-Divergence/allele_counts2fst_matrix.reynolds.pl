#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

use Getopt::Long;
use Data::Dumper;

# Option table and defaults

my $opt = {

	'group_table' => {} ,

	'group' => {} ,
	
	'input' => [] ,
	
	'window' => 0 ,
	
	'no_print_snps' => 0 ,
	
    'region' => {} ,
	
};

GetOptions(
    $opt ,

    # Input
    # -----

    # Allele file
    'input|in=s{,}' ,
    
    # Keep SNPs file
    'keep|k=s' ,
    
    'seqs=s' ,
    'coverage=s' ,
    
    'group|groups|g=s{,}' ,
    
    'window|windows=i' ,
    
    'region=s{,}' ,
    
    # Output
    # ------
    
    # Output name (string)
    'output|out|o=s' ,
    
    'verbose|v' ,
    
    'no_print_snps|i' ,

);

init_groups( $opt );
init_coverage( $opt );
init_keep( $opt );

parse_alleles( $opt );
print_matrix( $opt );

sub init_groups {
    
    my ( $opt ) = @_;

    foreach my $group ( sort keys %{ $opt->{ 'group' } } ) {
    
        my @samples = split ( /\,/ , $opt->{ 'group' }->{ $group } );
    
        foreach my $sample ( @samples ) {
    
            $opt->{ 'app' }->{ 'sample_table' }->{ $sample } = $group;
            
			push @{
				$opt->{ 'app' }->{ 'group_table' }->{ $group }->{ 'samples' }
			} , $sample;
            
        }
        
        push @{ $opt->{ 'groups' } } , $group;
    
    }
    
    my @groups = @{ $opt->{ 'groups' } };
    
    my $contrast_table = {};
    
    for ( my $i = 0; $i < @groups - 1; $i++ ) {
    
        my $group1_name = $groups[ $i ];
    
        for ( my $j = $i + 1; $j < @groups; $j++ ) {

			my $group2_name = $groups[ $j ];
                                    
            my $contrast =
            ${group1_name} . "_vs_" . ${group2_name};
            
            $contrast_table->{ $contrast } = 1;
            
        }
    
    }
    
    $opt->{ 'app' }->{ 'contrast_table' } = $contrast_table;
	
}

sub init_coverage {
    
    my ( $opt ) = @_;
    
    my $coverage_table = {};
    
    $opt->{ 'coverage_table' } = $coverage_table;
    
    if ( defined $opt->{ 'seqs' } and -e $opt->{ 'seqs' } ) {
    
        open ( my $in , "<" , $opt->{ 'seqs' } ) or die "$!";
    
        say "Reading sequences from $opt->{ 'seqs' } ...";
    
        while (<$in>) {
        
            chomp;
            
            if ( $_ =~ m/^(\S+)/ ) {
            
                $opt->{ 'get_chrom' }->{ $1 } = 1;
            
            }
        
        }
    
    }
    
    if ( defined $opt->{ 'coverage' } and -e $opt->{ 'coverage' } ) {
        
        my $coverage_file = $opt->{ 'coverage' };
        
        say "Reading coverage stats from $coverage_file ...";

        my $coverage_in;
        
        if ( $coverage_file =~ m/\.gz/ ) {
        
            open ( $coverage_in , "zcat $coverage_file |" ) or die "$!";
        
        }
        
        else {
        
            open ( $coverage_in , "<" , $coverage_file ) or die "$!";
        
        }

        my $tmp_header = undef;
        
        while (<$coverage_in>) {
        
			chomp;

			if ( $_ =~ m/^\>(.+)/ ) {
			
                $tmp_header = undef;
			
                if ( defined $opt->{ 'get_chrom' }->{ $1 } ) {
                
                    $tmp_header = $1;
                
                }
			
			}
			
			elsif ( defined $tmp_header ) {
			
                $coverage_table->{ $tmp_header } .= $_;
			
			}
			
        }
        
    }
    
}

sub init_keep {

    my ( $opt ) = @_;
    
    if ( defined $opt->{ 'keep' } ) {
    
		my $keep_table = {};
		
		$opt->{ 'app' }->{ 'keep_table' } = $keep_table;
    
		my $keep_file = $opt->{ 'keep' };
		
		say "Reading SNP positions from $keep_file";
		
		open ( my $keep_in , "<" , $keep_file ) or die "$!";
		
		while (<$keep_in>) {
		
			chomp;
			
			next if $_ =~ m/^#/;
			
			if ( $_ =~ m/^(\S+)\s(\d+)/ ) {
			
				$keep_table->{ $1 }->{ $2 } = 1;
				
			}
		
		}
	
	}

}

sub parse_alleles {
    
    my ( $opt ) = @_;
    
    my @samples;
    
    my $group_table = $opt->{ 'app' }->{ 'group_table' };
    
    my $coverage_table = $opt->{ 'coverage_table' };
    
    my $keep_table = $opt->{ 'app' }->{ 'keep_table' };
    
    my $positions = 0;
    my $last_chrom = "";
    my $chrom_i = 0;
    
    my $output = "";
    
    if ( defined $opt->{ 'output' } ) {
    
		$output = $opt->{ 'output' };
    
    }
    
    else {
    
		$output = $opt->{ 'input' }->[ 0 ];
    
    }
    
    my $contrast_table = $opt->{ 'app' }->{ 'contrast_table' };
    
    my @sorted_contrasts = sort keys %{ $contrast_table };
    
    my $out;
    my $num_out;
    
    unless ( $opt->{ 'no_print_snps' } ) {
    
        open ( $out , ">" , $output . ".fst.reynolds.snps.csv" ) or die 
        "Unable to create SNP FST output: $!";
    
        open ( $num_out , ">" , $output . ".fst.reynolds.snps.num_denom.csv" ) or die 
        "Unable to create SNP FST output: $!";
    
    }
    
    my $win_out;
    
    if ( $opt->{ 'window' } ) {
    
		my $window = $opt->{ 'window' };
    
        open ( $win_out , ">" , $output . ".fst.reynolds.windows.${window}bp.csv" ) or die 
        "Unable to create pairwise FST matrix: $!";
        
        say $win_out "CHROM\tPOS\t" , join (
            "\t" , map { "N_${_}\tFST_${_}" } @sorted_contrasts
        );
    
    }

    foreach my $alleles_file ( @{ $opt->{ 'input' } } ) {
    
		if ( -e $alleles_file ) {
			
			my $alleles_in;
			
			if ( $alleles_file =~ m/\.gz$/ ) {
			
				open ( $alleles_in , "gunzip -c $alleles_file |" ) or die
				"Unable to open ALLELES file $alleles_file: $!";
			
			}
			
			else {
			
				open ( $alleles_in , "<" , $alleles_file ) or die
				"Unable to open ALLELES file $alleles_file: $!";

			}
			
            my @groups;
            
            my $header = <$alleles_in>;
            chomp $header;
            
            my @headers = split ( /\t/ , $header );
            
            my @sample_names = @headers[ 5 .. $#headers ];
            
            my $i = 0;
            
            foreach my $sample ( @sample_names ) {
            
                if ( defined $opt->{ 'app' }->{ 'sample_table' }->{ $sample } ) {
            
                    my $group =
                    $opt->{ 'app' }->{ 'sample_table' }->{ $sample };
                
                    push @{ $group_table->{ $group }->{ 'ids' } } , $i;
            
                }
            
                $i++;
            
            }
            
            @groups = sort keys %{ $group_table };
            
            my $last_chrom = "";
            my $coverage_seq;
            
            my $window_table = {};
                
            while (<$alleles_in>) {
                    
                chomp;
                    
                if ( $_ =~ m/^(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(.+)/ ) {
                
                    my ( $chrom , $pos , $ref , $alt , $gts , $data ) =
                    ( $1 , $2 , $3 , $4 , $5 , $6 );
                    
                    # last if $chrom eq "seq_s_2";
                    
                    if ( $chrom ne $last_chrom ) {
                    
						if ( keys $coverage_table ) {
                    
							next unless defined $coverage_table->{ $chrom };
                    
						}
                    
                        $coverage_seq = $coverage_table->{ $chrom };
                        
                        if ( keys %{ $window_table } and $opt->{ 'window' } ) {
                        
                            print_windows(
                                $opt ,
                                $window_table ,
                                $contrast_table ,
                                $win_out ,
                            );
                        
                        }
                        
                        $window_table = {};
                    
                    }
                    
                    if ( defined $keep_table ) {
                    
                        next unless defined $keep_table->{ $chrom }->{ $pos };
                    
                    }
                    
                    if ( keys $coverage_table ) {
                    
						my $symbol =
						substr ( $coverage_seq , $pos - 1 , 1 );
                    
						next unless defined $opt->{ 'region' }->{ $symbol };

                    }
                    
                    # say "$chrom\t$pos\t$symbol";
                    
                    # last if $pos > 100_000;

                    my @allele_counts = split ( /\t/ , $data );
					
					# Pairwise FST and frequency differences
					# --------------------------------------

					if ( defined $out ) {
					
                        print $out "$chrom\t$pos\t$ref\t$alt";
					
					}
					
					for ( my $i = 0; $i < @groups - 1; $i++ ) {

						my $group1_name = $groups[ $i ];

                        my @group_alleles = @allele_counts[
                            @{ $group_table->{ $group1_name }->{ 'ids' } }
                        ];
                        
                        my ( $n1 , $n1_0 , $n1_1 , $p1_0 , $p1_1 );
                        
                        for ( my $i = 0; $i < @group_alleles; $i += 2 ) {

                            if ( defined $group_alleles[ $i ] and defined $group_alleles[ $i + 1 ] ) {
                        
                                $n1 += $group_alleles[ $i ];
                                $n1 += $group_alleles[ $i + 1 ];
                                
                                $n1_0 += $group_alleles[ $i ];
                                
                                $n1_1 += $group_alleles[ $i + 1 ];
                        
                            }
                            
                        }
                        
                        $p1_0 = $n1_0 / $n1;
                        $p1_1 = $n1_1 / $n1;
                        
                        my $group1_size = $n1;
                        
						for ( my $j = $i + 1; $j < @groups; $j++ ) {

							my $group2_name = $groups[ $j ];
                            
                            my @group_alleles = @allele_counts[
                                @{ $group_table->{ $group2_name }->{ 'ids' } }
                            ];
                            
                            my ( $n2 , $n2_0 , $n2_1 , $p2_0 , $p2_1 );
                            
                            for ( my $i = 0; $i < @group_alleles; $i += 2 ) {
                            
                                if ( defined $group_alleles[ $i ] and defined $group_alleles[ $i + 1 ] ) {
                            
                                    $n2 += $group_alleles[ $i ];
                                    $n2 += $group_alleles[ $i + 1 ];
                                    
                                    $n2_0 += $group_alleles[ $i ];
                                    
                                    $n2_1 += $group_alleles[ $i + 1 ];
                            
                                }
                                
                            }
                            
                            $p2_0 = $n2_0 / $n2;
                            $p2_1 = $n2_1 / $n2;
                            
                            my $group2_size = $n2;
							
							my $p_sum;
							my $p_diff;
							my $p_prod;
							
							# If alleles have been observed for both groups
							
							if ( $n1 and $n2 ) {
							
								# If there no is variation between them,
								# move to next contrast
							
								if (
								
									( $p1_0 == 0 or $p1_0 == 1 )
									and $p1_0 == $p2_0
								) {

								
                                    if ( defined $out ) {
								
                                        print $out "\t";
                                    
                                    }
                                    
									next;
								
								}
							
                                # Reynolds FST
							
								$p_sum += ( $p1_0 ** 2 + $p2_0 ** 2 );
								$p_diff += ( $p1_0 - $p2_0 ) ** 2;
								$p_prod += ( $p1_0 * $p2_0 );
							
								$p_sum += ( $p1_1 ** 2 + $p2_1 ** 2 );
								$p_diff += ( $p1_1 - $p2_1 ) ** 2;
								$p_prod += ( $p1_1 * $p2_1 );

								my $n = $n1 + $n2;
								
								my $numerator =
								
								( 0.5 * $p_diff ) -
								( 1 / ( 2 * ( 2 * $n - 1 ) ) ) *
								( 2 - $p_sum );
								
								my $denominator = 1 - $p_prod;
					
								# Group1 vs Group2
					
								$opt->{ 'app' }->{ 'contrast' }->
								{ $group1_name }->{ $group2_name }->
								{ 'numerator' } += $numerator;
								
								$opt->{ 'app' }->{ 'contrast' }->
								{ $group1_name }->{ $group2_name }->
								{ 'denominator' } += $denominator;
					
								# Group2 vs Group1
					
								$opt->{ 'app' }->{ 'contrast' }->
								{ $group2_name }->{ $group1_name }->
								{ 'numerator' } += $numerator;
								
								$opt->{ 'app' }->{ 'contrast' }->
								{ $group2_name }->{ $group1_name }->
								{ 'denominator' } += $denominator;
								
								if ( defined $num_out ) {
								
									say $num_out "$chrom\t$pos\t$group1_name\t$group2_name\t$numerator\t$denominator";
								
								}
								
								if ( $opt->{ 'window' } ) {
								
                                    my $window_i = int ( $pos / $opt->{ 'window' } );
                                    
                                    my $contrast =
                                    ${group1_name} . "_vs_" . ${group2_name};                                    
                                    
                                    $window_table->
                                    { $chrom }->[ $window_i ]->
                                    { $contrast }->
                                    { 'numerator' } += $numerator;
                                    
                                    $window_table->
                                    { $chrom }->[ $window_i ]->
                                    { $contrast }->
                                    { 'denominator' } += $denominator;
                                    
                                    $window_table->
                                    { $chrom }->[ $window_i ]->
                                    { $contrast }->
                                    { 'n' }++;
								
								}
								
								if ( defined $out ) {
								
                                    # Weir-Cockerham FST
                                    
                                    my $r = 2;
                                
                                    my $n0_tot = $n1_0 + $n2_0;
                                    my $n1_tot = $n1_1 + $n2_1;
                                    
                                    my $p1;
                                    my $p2;
                                    
                                    if ( $n0_tot >= $n1_tot ) {
                                    
                                        $p1 = $n1_0 / $n1;
                                        $p2 = $n2_0 / $n2;
                                    
                                    }
                                    
                                    else {
                                    
                                        $p1 = $n1_1 / $n1;
                                        $p2 = $n2_1 / $n2;
                                    
                                    }
                                    
                                    my $fst;

                                    if ( $p1 == $p2 ) {
                                        
                                        $fst = 0;
                                        
                                    }

                                    else {
                                    
                                        # Average across populations
                                        
                                        my $n_avg = ( $n1 + $n2 ) / 2;
                                        my $p_avg =
                                        ( $p1 * $n1 + $p2 * $n2 ) / ( $r * $n_avg );

                                        # Sample size variance
                                        
                                        my $s = (
                                            ( $n1 * ( $p1 - $p_avg ) ** 2 ) /
                                            ( ( $r - 1 ) * $n_avg )
                                        ) + (
                                            ( $n2 * ( $p2 - $p_avg ) ** 2 ) /
                                            ( ( $r - 1 ) * $n_avg )
                                        );
                                        
                                        # Standard deviation
                                        
                                        my $dev = sqrt(
                                            ( ( $n1 - $n_avg ) ** 2 +
                                            ( $n2 - $n_avg ) ** 2 ) / $r );
                                        
                                        # Squared coefficient of variation of sample sizes
                                        
                                        my $c = $dev / $n_avg;
                                        my $c_to_2 = $c**2; 
                                                                
                                        $fst = ( $s -
                                            ( 1 / ( 2 * $n_avg - 1 ) ) *
                                            (
                                                ( $p_avg * ( 1 - $p_avg ) ) -
                                                ( ( ( $r - 1 ) / $r ) * $s )
                                            )
                                        ) / (
                                                ( 1 - ( 2 * $n_avg * $c_to_2 ) /
                                                    ( ( 2 * $n_avg - 1 ) * $r )
                                                ) *
                                                ( $p_avg * ( 1 - $p_avg ) ) +
                                                ( 1 + ( 2 * $n_avg * ( $r - 1 ) * $c_to_2 ) /
                                                    ( ( 2 * $n_avg - 1 ) * $r )
                                                ) * ( $s / $r )
                                        );
                                        
                                        $fst = 0 if $fst < 0;
                                
                                    }
                                    
                                    $fst = 0 if $fst < 0;
                                
                                    $fst = sprintf( "%.4f" , $fst );
                                
                                    my $tot_size = $group1_size + $group2_size;
                                    
                                    my $tot_prop = ( $n1 + $n2 ) / $tot_size;
                                    $tot_prop = sprintf( "%.4f" , $tot_prop );
                                    
                                    my $n1_prop = $n1 / $group1_size;
                                    $n1_prop = sprintf( "%.4f" , $n1_prop );
                                    
                                    my $n2_prop = $n2 / $group2_size;
                                    $n2_prop = sprintf( "%.4f" , $n2_prop );
                                    
                                    print $out "\t${group1_name}/${group2_name}:$fst:";
                                    print $out "$tot_size:$tot_prop:$n1_prop:$n2_prop";
                                    
                                    $n1_0 = sprintf( "%.4f" , $n1_0 );
                                    $n1_1 = sprintf( "%.4f" , $n1_1 );
                                    $p1_0 = sprintf( "%.4f" , $p1_0 );
                                    $p1_1 = sprintf( "%.4f" , $p1_1 );
                                    
                                    $n2_0 = sprintf( "%.4f" , $n2_0 );
                                    $n2_1 = sprintf( "%.4f" , $n2_1 );
                                    $p2_0 = sprintf( "%.4f" , $p2_0 );
                                    $p2_1 = sprintf( "%.4f" , $p2_1 );
                                    
                                    print $out "|$group1_name,$n1,$n1_0,$n1_1,$p1_0,$p1_1|$group2_name,$n2,$n2_0,$n2_1,$p2_0,$p2_1";
                                    
                                    # say "$group1_name|$n1|$n1_0|$n1_1|$p1_0|$p1_1|$group2_name|$n2|$n2_0|$n2_1|$p2_0|$p2_1";
								
								}
								
							}
							
							else {
							
                                print $out "\t";
							
							}

						}
						
					}
					
					say $out "" if defined $out;
					
					$positions++;
					
					if ( $chrom ne $last_chrom ) {
					
						$chrom_i++;
					
						say "$chrom ($chrom_i) $positions snps parsed so far ...";
					
					}
					
					$last_chrom = $chrom;
					
				}
				
			}
			
            if ( keys %{ $window_table } and $opt->{ 'window' } ) {
            
                print_windows(
                    $opt ,
                    $window_table ,
                    $contrast_table ,
                    $win_out ,
                );
                
                $window_table = {};
            
            }
			
		}
		
		else {
		
            die "No such file $alleles_file: $!";
		
		}
		
	}
	
	say "$last_chrom ($chrom_i) $positions snps parsed in total";
    
}

sub print_matrix {
    
    my ( $opt ) = @_;
    
    my $output = "";
    
    if ( defined $opt->{ 'output' } ) {
    
		$output = $opt->{ 'output' };
    
    }
    
    else {
    
		$output = $opt->{ 'input' }->[ 0 ];
    
    }
    
    my $contrast_table =
    $opt->{ 'app' }->{ 'contrast' };
    
    # PHYLIP
    
    open ( my $out , ">" , $output . ".fst.reynolds.phy" ) or die 
    "Unable to create pairwise FST matrix: $!";

	my @groups = @{ $opt->{ 'groups' } };
    
    say $out scalar @groups;
    
    # say Dumper ( $contrast_table );
    
    foreach my $group1_name ( @groups ) {
        
        print $out $group1_name;
        
        foreach my $group2_name ( @groups ) {
            
            print $out "\t";
            
            if ( $group1_name eq $group2_name ) {
                
                print $out 0;
                
            }

            elsif (
            
                defined $contrast_table->{ $group1_name }->{ $group2_name }->
                { 'denominator' }
            
            ) {
                
                # say "$group1_name\t$group2_name";
                
                print $out sprintf( "%.4f",
                $contrast_table->{ $group1_name }->{ $group2_name }->
                { 'numerator' } /
                $contrast_table->{ $group1_name }->{ $group2_name }->
                { 'denominator' } );
            
            }
            
        }
        
        say $out "";
        
    }
    
    close $out;
    
    # CSV
    
    open ( $out , ">" , $output . ".fst.reynolds.csv" ) or die 
    "Unable to create pairwise FST matrix: $!";
    
    say $out "\t" , join ( "\t" , @groups );
    
    foreach my $group1_name ( @groups ) {

        print $out $group1_name;

        foreach my $group2_name ( @groups ) {
            
            print $out "\t";
            
            if ( $group1_name eq $group2_name ) {
                
                print $out 0;
                
            }

            elsif (
            
                defined $contrast_table->{ $group1_name }->{ $group2_name }->
                { 'denominator' }
            
            ) {
                
                print $out sprintf( "%.4f", 
                $contrast_table->{ $group1_name }->{ $group2_name }->
                { 'numerator' } /
                $contrast_table->{ $group1_name }->{ $group2_name }->
                { 'denominator' } );
            
            }
            
        }
        
        say $out "";
        
    }    
    
}

sub print_windows {
    
    my ( $opt , $window_table , $contrast_table , $out ) = @_;

    my @sorted_contrasts = sort keys %{ $contrast_table };
    
    my @sorted_chromosomes = sort keys %{ $window_table };
    
    foreach my $chrom ( @sorted_chromosomes ) {
    
        for ( my $i = 0; $i < @{ $window_table->{ $chrom } }; $i++ ) {
    
            print $out "$chrom\t" , $i * $opt->{ 'window' };
    
            foreach my $contrast ( @sorted_contrasts ) {
    
                my $fst = "";
                my $n = "";
                
                if ( defined $window_table->{ $chrom }->[ $i ]-> { $contrast } ) {
                    
                    $fst = sprintf( "%.4f", 
                    $window_table->{ $chrom }->[ $i ]-> { $contrast }->
                    { 'numerator' } /
                    $window_table->{ $chrom }->[ $i ]-> { $contrast }->
                    { 'denominator' } );
                
                    $n = $window_table->{ $chrom }->[ $i ]-> { $contrast }->
                    { 'n' };
                
                }
                
                print $out "\t$n\t$fst";
            
            }
            
            say $out "";
            
        }

    }
    
}

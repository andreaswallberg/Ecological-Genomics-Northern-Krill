#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)
# See note for the BioPerl subroutine tajima_D_counts_simple

use strict;
use warnings;
use 5.010;

use Getopt::Long;
use Data::Dumper;

# Option table and defaults

my $opt = {
    
    # Input files
    
    'vcf' => [] ,

    # Groups of samples

    'groups' => [] ,

    'group' => {} ,
    
    'window' => [] ,
    
    'region' => {} ,
    
};

GetOptions(
    $opt ,

    'group|groups|g=s{,}' ,
    
    'genes|gtf|gff=s' ,
    
    'coverage|c=s' ,
    
    'seqs|s=s' ,
    
    'vcf|input=s{,}' ,
    
    'window|windows=i{,}' ,
    
    'region=s{,}' ,
    
    'output|out|o=s' ,
   
    # Booleans
    
    'verbose|v' ,
    'help|h' ,

);

# Possibly print help
# -------------------

#pod2usage( { -exitval => 1, -verbose => 2, -output => \*STDOUT } ) if not
#scalar @ARGV;

=head1 SYNOPSIS

 watterson_pi.pl [options]

 Options and switches:

=cut

# Possibly print help
# -------------------

pod2usage( { -exitval => 1, -verbose => 2, -output => \*STDOUT } ) if
$opt->{ 'help' };

# Set defaults
# ------------

init_defaults( $opt );

# Create groups
# -------------

init_groups( $opt );

# Read coverage data
# ------------------

init_coverage( $opt );

# Read genes
# ------------------

init_genes( $opt );

# Parse phased data in vcf format
# ----------------------------------

init_vcf( $opt );

# Print results
# -------------

foreach my $symbol ( keys %{ $opt->{ 'region' } } ) {

    my $region = $opt->{ 'region' }->{ $symbol };

    next unless defined $region;
    
    print_stats( $opt , $symbol , $region );
    
    foreach my $window ( @{ $opt->{ 'window' } } ) {
    
        print_stats_windows( $opt , $symbol , $region , $window );

    }
    
}

sub init_defaults {
    
    my ( $opt ) = @_;
    
    if ( defined $opt->{ 'vcf' } and not defined $opt->{ 'output' } ) {
        
        $opt->{ 'output' } = $opt->{ 'vcf' }->[ 0 ];
        
    }
    
}

sub init_groups {
    
    my ( $opt ) = @_;

    foreach my $group ( keys %{ $opt->{ 'group' } } ) {
    
        # Parse and generate groups from the CSV file
        
        my $group_file = $opt->{ 'group' }->{ $group };

        open ( my $in , "<" , $group_file ) or die "$!";
        
        say "Reading group file $group_file ..." if $opt->{ 'verbose' };
        
        while (<$in>) {
        
            chomp;
            
            my $sample = $_;
            
            $opt->{ 'app' }->{ 'sample_table' }->{ $sample } = $group;

            push @{
                $opt->{ 'app' }->{ 'group_table' }->{ $group }->{ 'samples' }
            } , $sample;
        
        }
	
	}
	
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
            
                unless ( defined $opt->{ 'get_chrom' }->{ $1 } ) {
            
                    $opt->{ 'get_chrom' }->{ $1 } = 1;
                    
                    push @{ $opt->{ 'get_chroms' } } , $1;
            
                }
            
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

sub init_genes {

    my ( $opt ) = @_;
    
    my $gene_table = {};
    
    $opt->{ 'app' }->{ 'gene_table' } = $gene_table;
    
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
                
                foreach my $window ( @{ $opt->{ 'window' } } ) {
                
                    push @{ $gene_table->{ $window }->{ $chrom } } , [
                        ( int ( $start / $window ) ),
                        ( int ( $stop / $window ) ) ,
                        $start ,
                        $stop ,
                        $orientation ,
                    ];
                    
                    if ( $opt->{ 'verbose' } ) {
                    
                        # say "Adding gene at $chrom $start $stop ($orientation)";
                    
                    }
                
                }
            
            }
			
        }
        
    }

}

sub init_vcf {
    
    my ( $opt ) = @_;
       
    # Read VCF
    # --------
    
    foreach my $vcf_file ( @{ $opt->{ 'vcf' } } ) {
        
        parse_vcf( $opt , $vcf_file );
        
    }
    
}

sub parse_vcf {
    
    my ( $opt , $vcf_file ) = @_;
    
    my $coverage_table = $opt->{ 'coverage_table' };
    
    my @windows = @{ $opt->{ 'window' } };
        
    if ( -e $vcf_file ) {
        
        my $vcf_in;
        
        if ( $vcf_file =~ m/\.gz/ ) {
        
            open ( $vcf_in , "zcat $vcf_file |" ) or die "$!";
        
        }
        
        else {
        
            open ( $vcf_in , "<" , $vcf_file ) or die "$!";
        
        }
        
        if ( $opt->{ 'verbose' } ) {
            
            say "==> Reading from $vcf_file";
            
        }
        
        my @sample_names;
        
        my $group_table = {};
        
        $opt->{ 'group_table' } = $group_table;
        
        my @groups;
        
        while (<$vcf_in>) {
            
            chomp;
            
            if ( $_ =~ m/^\#\#/ ) {
            
                next;
            
            }
            
            elsif ( $_ =~ m/^#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT\s+(.+)/ ) {

                @sample_names = ( split ( /\s+/ , $1 ) );

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
                
                foreach my $group ( @groups ) {
                
                    my $samples = @{ $group_table->{ $group }->{ 'ids' } };
        
                    my $chromosomes = 2 * $samples;
                    
                    my $a_N = 0;
                    my $a2 = 0;
                    
                    for ( my $i = 1; $i < $chromosomes; $i++ ) {
                        
                        $a_N += 1 / $i;
                        $a2 += ( 1 / $i**2 );
                        
                    }

                    $group_table->{ $group }->{ 'a_N' } = $a_N;
                    $group_table->{ $group }->{ 'a_2' } = $a2;
                
                }
                
                last;
                
            }
            
        }
        
        my $last_chrom = "";
        
        my $coverage_seq;
        
        my $sample_seq_table;
        my $sample_win_table;
        my $variation_table;
        
        while (<$vcf_in>) {
            
            chomp;
            
            if ( $_ =~m/^(\S+)\s+(\d+)\s+\S+\s+(\S)\s+(\S)\s+\S+\s+.+\S+\s+GT\S*\s+(.+)/ ) {
                
                my ( $chrom , $pos , $ref , $alt , $data ) = ( $1  , $2 , $3 , $4 , $5 );
                
                # last if $chrom eq "seq_s_30";
                
                # say "$chrom\t$pos";
                
                # last if $pos > 100_000;
                
                if ( $chrom ne $last_chrom ) {
                
                    next unless defined $coverage_table->{ $chrom };
                
                    $coverage_seq = $coverage_table->{ $chrom };
                    
                    # last if $chrom eq "seq_s_2";
                    
                    say "$chrom\t$pos";
                
                }
                
                my $symbol =
                substr ( $coverage_seq , $pos - 1 , 1 );
                
                next unless defined $opt->{ 'region' }->{ $symbol };
                
                if ( defined $opt->{ 'region_table' }->{ $symbol } ) {
                
                    $sample_seq_table =
                    $opt->{ 'region_table' }->
                    { $symbol }->{ 'sample_seq_table' };
                    
                    $sample_win_table =
                    $opt->{ 'region_table' }->
                    { $symbol }->{ 'sample_win_table' };
                    
                    $variation_table =
                    $opt->{ 'region_table' }->
                    { $symbol }->{ 'variation_table' };
                
                }
                
                else {
                
                    $sample_seq_table = {};
                    
                    $opt->{ 'region_table' }->
                    { $symbol }->{ 'sample_seq_table' } =
                    $sample_seq_table;
                    
                    $sample_win_table = {};
                    
                    $opt->{ 'region_table' }->
                    { $symbol }->{ 'sample_win_table' } =
                    $sample_win_table;
                    
                    $variation_table = {};
                    
                    $opt->{ 'region_table' }->
                    { $symbol }->{ 'variation_table' } =
                    $variation_table;
                
                }
                
                my $tmp_win_table = {};
                
                foreach my $window ( @windows ) {
                
                    $tmp_win_table->{ $window } = int ( $pos / $window );
                
                }
                
                my @gts = split ( /\t/ , $data );
                
                foreach my $group ( @groups ) {
                    
                    my @group_gts = @gts[
						@{ $group_table->{ $group }->{ 'ids' } }
					];

					my $group_allele_table = {};

					my $i = 0;
					
					foreach my $gt ( @group_gts ) {
					
                        if ( $gt =~ m/^0\W0/ ) {
                        
                            $group_allele_table->{ $ref } += 2;
                            
                            $sample_seq_table->{ 'ALL' }->{ $group }->[ $i ] .= $ref;
                            $sample_seq_table->{ 'ALL' }->{ $group }->[ $i + 1 ] .= $ref;
                        
                            $sample_seq_table->{ $chrom }->{ $group }->[ $i ] .= $ref;
                            $sample_seq_table->{ $chrom }->{ $group }->[ $i + 1 ] .= $ref;
                        
                            foreach my $window ( @windows ) {
                        
                                $sample_win_table->{ $window }->{ $chrom }->
                                [ $tmp_win_table->{ $window } ]->{ $group }->[ $i ] .= $ref;
                                
                                $sample_win_table->{ $window }->{ $chrom }->
                                [ $tmp_win_table->{ $window } ]->{ $group }->[ $i + 1 ] .= $ref;
                            
                            }
                        
                        }
                        
                        elsif ( $gt =~ m/^1\W1/ ) {
                        
                            $group_allele_table->{ $alt } += 2;
                        
                            $sample_seq_table->{ 'ALL' }->{ $group }->[ $i ] .= $alt;
                            $sample_seq_table->{ 'ALL' }->{ $group }->[ $i + 1 ] .= $alt;
                            
                            $sample_seq_table->{ $chrom }->{ $group }->[ $i ] .= $alt;
                            $sample_seq_table->{ $chrom }->{ $group }->[ $i + 1 ] .= $alt;
                            
                            foreach my $window ( @windows ) {
                            
                                $sample_win_table->{ $window }->{ $chrom }->
                                [ $tmp_win_table->{ $window } ]->{ $group }->[ $i ] .= $alt;
                                
                                $sample_win_table->{ $window }->{ $chrom }->
                                [ $tmp_win_table->{ $window } ]->{ $group }->[ $i + 1 ] .= $alt;
                        
                            }
                        
                        }
                        
                        elsif ( $gt =~ m/^0\W1/ or $gt =~ m/^1\W0/ ) {
                        
                            $group_allele_table->{ $ref }++;
                            $group_allele_table->{ $alt }++;
                        
                            # Do not keep track on haplotype for the heterozygous
                            # genotypes here
                        
                            $sample_seq_table->{ 'ALL' }->{ $group }->[ $i ] .= $ref;
                            $sample_seq_table->{ 'ALL' }->{ $group }->[ $i + 1 ] .= $alt;
                        
                            $sample_seq_table->{ $chrom }->{ $group }->[ $i ] .= $ref;
                            $sample_seq_table->{ $chrom }->{ $group }->[ $i + 1 ] .= $alt;
                            
                            foreach my $window ( @windows ) {
                            
                                $sample_win_table->{ $window }->{ $chrom }->
                                [ $tmp_win_table->{ $window } ]->{ $group }->[ $i ] .= $ref;
                                
                                $sample_win_table->{ $window }->{ $chrom }->
                                [ $tmp_win_table->{ $window } ]->{ $group }->[ $i + 1 ] .= $alt;
                        
                            }
                        
                        }
                        
                        else {
                        
                            $group_allele_table->{ '?' }++;
                            $group_allele_table->{ '?' }++;
                        
                            # Do not keep track on haplotype for the heterozygous
                            # genotypes here
                        
                            $sample_seq_table->{ 'ALL' }->{ $group }->[ $i ] .= "?";
                            $sample_seq_table->{ 'ALL' }->{ $group }->[ $i + 1 ] .= "?";
                        
                            $sample_seq_table->{ $chrom }->{ $group }->[ $i ] .= "?";
                            $sample_seq_table->{ $chrom }->{ $group }->[ $i + 1 ] .= "?";
                            
                            foreach my $window ( @windows ) {
                            
                                $sample_win_table->{ $window }->{ $chrom }->
                                [ $tmp_win_table->{ $window } ]->{ $group }->[ $i ] .= "?";
                                
                                $sample_win_table->{ $window }->{ $chrom }->
                                [ $tmp_win_table->{ $window } ]->{ $group }->[ $i + 1 ] .= "?";
                        
                            }
                        
                        }
                        
                        $i += 2;
					
					}
                    
					# For pi
					# ------
                    
                    if ( defined $group_allele_table->{ '?' } ) {
                    
						delete $group_allele_table->{ '?' };
                    
                    }
                    
                    # Keep track on variable SNPs in the group
                    
                    if ( keys %{ $group_allele_table } == 2 ) {
                    
						# say "$group\t@group_alleles";
                    
						# Number of samples with data
						
						$group_allele_table->{ $ref } = 0 unless defined
						$group_allele_table->{ $ref };
						
						$group_allele_table->{ $alt } = 0 unless defined
						$group_allele_table->{ $alt };
						
						my $samples = $group_allele_table->{ $ref };
						$samples = $group_allele_table->{ $alt };
						
						$samples /= 2;
                    
                        # Register the number of samples
                        
                        $variation_table->{ 'chrom' }->{ 'ALL' }->
                        { $group }->{ 'samples' } += $samples;
                        
                        $variation_table->{ 'chrom' }->{ 'ALL' }->
                        { $group }->{ 'samples_n' }++;
                    
                        $variation_table->{ 'chrom' }->{ 'ALL' }->
                        { $group }->{ 'variable' }++;
                        
                        $variation_table->{ 'chrom' }->{ $chrom }->
                        { $group }->{ 'samples' } += $samples;
                        
                        $variation_table->{ 'chrom' }->{ $chrom }->
                        { $group }->{ 'samples_n' }++;
                    
                        $variation_table->{ 'chrom' }->{ $chrom }->
                        { $group }->{ 'variable' }++;
                        
                        foreach my $window ( @windows ) {
                        
                            $variation_table->{ 'win' }->{ $window }->{ $chrom }->
                            [ $tmp_win_table->{ $window } ]->
                            { $group }->{ 'samples' } += $samples;
                            
                            $variation_table->{ 'win' }->{ $window }->{ $chrom }->
                            [ $tmp_win_table->{ $window } ]->
                            { $group }->{ 'samples_n' }++;
                        
                            $variation_table->{ 'win' }->{ $window }->{ $chrom }->
                            [ $tmp_win_table->{ $window } ]->
                            { $group }->{ 'variable' }++;
                    
                        }
                    
                    }
                    
                }
                
                $last_chrom = $chrom;
        
            }

        }
            
        close $vcf_in;
        
    }
    
    else {
        
        die "Unable to locate VCF formatted input file $vcf_file: $!";
        
    }
    
}

sub print_stats {
    
    my ( $opt , $symbol , $region ) = @_;
    
    my $sample_seq_table =
    $opt->{ 'region_table' }->
    { $symbol }->{ 'sample_seq_table' };
    
    my $variation_table =
    $opt->{ 'region_table' }->
    { $symbol }->{ 'variation_table' }->{ 'chrom' };
    
	my $coverage_table = $opt->{ 'coverage_table' };

	my $covered_table = {};
	
	my $length_table = {};
	
	foreach my $chrom ( @{ $opt->{ 'get_chroms' } } ) {
	
        next if $chrom eq "ALL";
	
        next unless defined $variation_table->{ $chrom };
	
        my $chrom_seq = $coverage_table->{ $chrom };
	
        $length_table->{ $chrom } = length $chrom_seq;
        $length_table->{ 'ALL' } += $length_table->{ $chrom };
	
        $covered_table->{ $chrom } = 0;
	
        while ( $chrom_seq =~ m/(${symbol}+)/g ) {
	
            my $start = $-[0];
            my $stop = $+[0];
            
            $stop--;
            
            my $covered_length = $stop - $start + 1;
	
            $covered_table->{ $chrom } += $covered_length;
            $covered_table->{ 'ALL' } += $covered_length;
	
        }
        
        $covered_table->{ $chrom } = 0 unless defined
        $covered_table->{ $chrom };
	
	}
	
    $covered_table->{ 'ALL' } = 0 unless defined
    $covered_table->{ 'ALL' };
	
	my @groups = sort keys %{ $opt->{ 'app' }->{ 'group_table' } };

	open (
		my $out ,
		">" ,
		$opt->{ 'output' } . ".wattersons_theta_pi_tajimas_D.${region}.csv" ,
	);
	
	print $out "CHROM\tLENGTH\tCOVERED\tCOVERED_PROP";
    
    foreach my $group ( @groups ) {
    
        print $out "\t${group}_N\t${group}_VARIABLE\t${group}_THETA\t${group}_PI\t${group}_TD";
    
    }
    
    say $out "";
    
    my @sorted_chromosomes = @{ $opt->{ 'get_chroms' } };
    
    unshift @sorted_chromosomes , "ALL";
    
	my $group_table = $opt->{ 'group_table' };
    
    foreach my $chrom ( @sorted_chromosomes ) {
        
		my $chrom_length = $length_table->{ $chrom };
        my $covered_length = $covered_table->{ $chrom };
        
        next unless defined $variation_table->{ $chrom };
        
        print $out "$chrom\t$chrom_length\t$covered_length\t" , 
        $covered_length / $chrom_length;
        
        foreach my $group ( @groups ) {

            if ( defined
            
                $variation_table->{ $chrom }->
                { $group }->{ 'variable' }
            
            ) {

				my $avg_chrom_length = 0;
				my $avg_chrom_length_n = 0;

                my $variable = 0;
                my $watterson = "";
    
                my $pi = 0;
                my $tajimasd = 0;
            
                $variable =
                $variation_table->{ $chrom }->
                { $group }->{ 'variable' };
                                
                my $samples = @{ $group_table->{ $group }->{ 'ids' } };
                
                my $a_N = $group_table->{ $group }->{ 'a_N' };
                my $a2 = $group_table->{ $group }->{ 'a_2' };
                          
                $watterson = $variable / ( $covered_length * $a_N );
                          
                # Pi & Tajimas' D
                # ---------------
                                                
                my $nr_samples = @{ $sample_seq_table->{ $chrom }->{ $group } };
                
                for ( my $i = 0; $i < $nr_samples; $i++ ) {
                
                    $sample_seq_table->{ $chrom }->{ $group }->[ $i ] =
                    uc
                    $sample_seq_table->{ $chrom }->{ $group }->[ $i ];
            
                    $sample_seq_table->{ $chrom }->{ $group }->[ $i ] =~
                    s/[NRYSWJNBDHV\.\?\-]/Z/g;
                
                }
                
                my $alleles_freq = 1 / $nr_samples;
                
                for (
                
                    my $i = 0;
                    $i < $nr_samples - 1;
                    $i++
                
                ) {

					my $seq = $sample_seq_table->{ $chrom }->{ $group }->[ $i ];

					my $seq_length = length $seq;

					my $n_unknown = $seq =~ tr/Z//;
					
					# say "$seq_length\t$n_unknown";
					
					for (
                    
                        my $j = $i + 1;
                        $j < $nr_samples;
                        $j++
                    
                    ) {
					
						my $tmp_seq = $seq;
					
						# say "Comparing sample $i to $j for sequence $chrom ...";
					
						my $other_alleles_freq = 1 / $nr_samples;
					
						my $other_seq = $sample_seq_table->{ $chrom }->{ $group }->[ $j ];

						my $tmp_other_seq = $other_seq;
						
						my $length_other_seq = length $tmp_other_seq;
						my $other_n_unknown = $tmp_other_seq =~ tr/Z//;
						
						# say "$length_other_seq\t$other_n_unknown";
						
						my $diff = $tmp_seq ^ $tmp_other_seq;

						( my $mask = $diff ) =~ tr{\x00}{\xff}c;

						$tmp_seq &= $mask;
						$tmp_other_seq &= $mask;
							
						my $n_differences = $tmp_seq =~ tr/[ACGT]//;
						my $n_unknown_diff_other = $tmp_other_seq =~ tr/Z//;
							
                        my $tot = $covered_length;
							
						# my $tot = $seq_length - $n_unknown - $n_unknown_diff_other;
							
						if ( $tot ) {
							
							my $differences =
							$n_differences - $n_unknown_diff_other;
                        
							$avg_chrom_length += $tot;
							$avg_chrom_length_n++;
							
							$pi +=
							(

								$alleles_freq *
								$other_alleles_freq *
								( $differences / $tot )

							);
							
							# say "Differences: $differences " , $differences / $tot;
                        
                        }
                    
                    }
                    
                    # exit;
                    
                }
                
                $pi *= 2;
                $pi *= ( $nr_samples / ( $nr_samples -1 ) );
                
                # $avg_chrom_length /= $avg_chrom_length_n;
                
                $tajimasd = "";
                
                if ( $nr_samples and $variable and $pi and $a_N and $a2 ) {
                
                    # Pi: 0.04944 Theta: 2.16841 Tajima's D: -2.21623 SegSites: 18
                
# 					$tajimasd = tajima_D_counts_simple(
# 						7 ,
# 						18 ,
# 						0.04944,
# 						$a_N ,
# 						$a2 ,
# 					);
                
					$tajimasd = tajima_D_counts_simple(
						$nr_samples ,
						$variable ,
						$pi * $covered_length ,
						$a_N ,
						$a2 ,
					);
                
                }
                
				print $out "\t$nr_samples\t$variable\t$watterson\t$pi\t$tajimasd";
                    
            }
                
        }
        
        say $out "";
        
    }

}

sub print_stats_windows {
    
    my ( $opt , $symbol , $region , $window ) = @_;
    
    my $sample_win_table =
    $opt->{ 'region_table' }->
    { $symbol }->{ 'sample_win_table' }->{ $window };
    
    my $variation_table =
    $opt->{ 'region_table' }->
    { $symbol }->{ 'variation_table' }->{ 'win' }->{ $window };
    
	my $coverage_table = $opt->{ 'coverage_table' };
	
	my $gene_table = $opt->{ 'app' }->{ 'gene_table' };

	my $length_table = {};
	my $covered_table = {};
	
	foreach my $chrom ( keys %{ $variation_table } ) {
	
        my $coverage_seq = $coverage_table->{ $chrom };
	
        my $length = length $coverage_seq;

        my $i = 0;
        
        for ( my $pos = 0; $pos < $length; $pos += $window ) {
        
            my $win_length = $window;
            
            if ( $pos + $win_length > $length ) {
            
                $win_length = $length - $pos;
            
            }
        
            my $coverage_subseq = substr ( $coverage_seq , $pos , $win_length );
            
            # my $covered_length = 
            
            while ( $coverage_subseq =~ m/(${symbol}+)/g ) {

                my $start = $-[0];
                my $stop = $+[0];
                
                $stop--;
                
                my $covered_length = $stop - $start + 1;

                $covered_table->{ $chrom }->[ $i ] += $covered_length;

            }
            
            $covered_table->{ $chrom }->[ $i ] = 0 unless defined
            $covered_table->{ $chrom }->[ $i ];
            
            $length_table->{ $chrom }->[ $i ] = $win_length;
        
            $i++;
        
        }
	
	}
	
	my @groups = sort keys %{ $opt->{ 'app' }->{ 'group_table' } };

	open (
		my $out ,
		">" ,
		$opt->{ 'output' } . ".wattersons_theta_pi_tajimas_D.window_${window}.${region}.csv" ,
	);
	
	print $out "CHROM\tPOS\tGLOBAL_POS\tLENGTH\tCOVERED\tCOVERED_PROP";
    
    foreach my $group ( @groups ) {
    
        print $out "\t${group}_N\t${group}_VARIABLE\t${group}_THETA\t${group}_PI\t${group}_TD";
    
    }
    
    say $out "";
    
    my @sorted_chromosomes = @{ $opt->{ 'get_chroms' } };
	
	my $group_table = $opt->{ 'group_table' };
    
    my $dist_table = {};
    
    my $global_pos = 0;
    
    foreach my $chrom ( @sorted_chromosomes ) {
        
        next unless defined $covered_table->{ $chrom };
        
        my $pos = 0;
        
        for ( my $win_i = 0; $win_i < @{ $covered_table->{ $chrom } }; $win_i ++ ) {
        
            my $chrom_length = $length_table->{ $chrom }->[ $win_i ];
            my $covered_length = $covered_table->{ $chrom }->[ $win_i ];
        
            print $out "$chrom\t" , "$pos\t$global_pos\t" , "$chrom_length\t$covered_length\t" , 
            $covered_length / $chrom_length;
            
            my $min_dist;
            my $min_dist_i;
            my $min_dist_relation;
            
            my $min_dist_valid = 0;
            
            if (
                
                defined $gene_table->{ $window }->{ $chrom } and
                ( $region eq "intergenic" or $region eq "any" )
            
            ) {
            
                my $avg_pos = $pos + $chrom_length / 2;
            
                foreach my $gene_ref (
                
                    @{ $gene_table->{ $window }->{ $chrom } }
                
                ) {
                
                    my $start = $gene_ref->[ 2 ];
                    my $stop = $gene_ref->[ 3 ];
                    
                    my $dist;
                    my $dist_relation;
                    
                    if ( $avg_pos < $start ) {
                    
                        $dist = $start - $avg_pos;
                        $dist_relation = "upstream";
                    
                    }
                    
                    elsif ( $avg_pos > $stop ) {
                    
                        $dist = $avg_pos - $stop;
                        $dist_relation = "downstream";
                    
                    }
                    
                    if ( defined $dist ) {
                    
                        if ( not defined $min_dist or $dist < $min_dist ) {
                    
                            $min_dist = $dist;
                            $min_dist_relation = $dist_relation;
                            
                            if (
                            
                                $gene_ref->[ 4 ] eq "+" or
                                $gene_ref->[ 4 ] eq "-"
                            
                            ) {
                        
                                $min_dist_valid = 1;
                        
                            }
                            
                            else {
                            
                                $min_dist_valid = 0;
                            
                            }
                            
                            if ( $gene_ref->[ 4 ] eq "-" ) {
                            
                                if ( $min_dist_relation eq "upstream" ) {
                                
                                    $min_dist_relation = "downstream";
                                
                                }
                                
                                elsif ( $min_dist_relation eq "downstream" ) {
                                
                                    $min_dist_relation = "upstream";
                                
                                }
                            
                            }
                        
                        }
                    
                    }
                
                }
                
                if ( defined $min_dist and $min_dist_valid ) {
                
                    $min_dist_i = int ( $min_dist / $window );
                    
                }
            
            }
            
            $pos += $chrom_length;
            $global_pos += $chrom_length;
        
            foreach my $group ( @groups ) {

                my $nr_samples = '';
                my $variable = '';
                my $watterson = '';
                my $pi = '';
                my $tajimasd = '';
                
                if ( $sample_win_table->{ $chrom }->[ $win_i ]->{ $group } and $covered_length ) {
            
                    $variation_table->{ $chrom }->[ $win_i ]->
                    { $group }->{ 'variable' } = 0
                    unless defined
                    $variation_table->{ $chrom }->[ $win_i ]->
                    { $group }->{ 'variable' };
                    
                    $variation_table->{ $chrom }->[ $win_i ]->
                    { $group }->{ 'samples' } = 0
                    unless defined
                    $variation_table->{ $chrom }->[ $win_i ]->
                    { $group }->{ 'samples' };
                    
                    $variation_table->{ $chrom }->[ $win_i ]->
                    { $group }->{ 'samples_n' } = 0
                    unless defined
                    $variation_table->{ $chrom }->[ $win_i ]->
                    { $group }->{ 'samples_n' };

                    my $avg_chrom_length = 0;
                    my $avg_chrom_length_n = 0;
                    
                    $variable =
                    $variation_table->{ $chrom }->[ $win_i ]->
                    { $group }->{ 'variable' };
                            
                    my $samples = @{ $group_table->{ $group }->{ 'ids' } };
                        
                    my $a_N = $group_table->{ $group }->{ 'a_N' };
                    my $a2 = $group_table->{ $group }->{ 'a_2' };
                                
                    $watterson = $variable / ( $covered_length * $a_N );
                                
                    # Pi & Tajimas' D
                    # ---------------
                    
                    my $nr_samples =
                    @{ $sample_win_table->{ $chrom }->[ $win_i ]->{ $group } };
                    
                    for ( my $i = 0; $i < $nr_samples; $i++ ) {
                    
                            $sample_win_table->{ $chrom }->[ $win_i ]->{ $group }->[ $i ] =
                            uc
                            $sample_win_table->{ $chrom }->[ $win_i ]->{ $group }->[ $i ];
                    
                            $sample_win_table->{ $chrom }->[ $win_i ]->{ $group }->[ $i ] =~
                            s/[NRYSWJNBDHV\.\?\-]/Z/g;
                    
                    }
                    
                    my $alleles_freq = 1 / $nr_samples;
                    
                    for (
                    
                        my $i = 0;
                        $i < $nr_samples - 1;
                        $i++
                    
                    ) {

                        my $seq = $sample_win_table->{ $chrom }->[ $win_i ]->{ $group }->[ $i ];

                        my $seq_length = length $seq;

                        my $n_unknown = $seq =~ tr/Z//;
                        
                        # say "$seq_length\t$n_unknown";
                        
                        for (
                        
                            my $j = $i + 1;
                            $j < $nr_samples;
                            $j++
                        
                        ) {
                        
                            my $tmp_seq = $seq;
                        
                            # say "Comparing sample $i to $j for sequence $chrom ...";
                        
                            my $other_alleles_freq = 1 / $nr_samples;
                        
                            my $other_seq = $sample_win_table->{ $chrom }->[ $win_i ]->{ $group }->[ $j ];

                            my $tmp_other_seq = $other_seq;
                            
                            my $length_other_seq = length $tmp_other_seq;
                            my $other_n_unknown = $tmp_other_seq =~ tr/Z//;
                            
                            # say "$length_other_seq\t$other_n_unknown";
                            
                            my $diff = $tmp_seq ^ $tmp_other_seq;

                            ( my $mask = $diff ) =~ tr{\x00}{\xff}c;

                            $tmp_seq &= $mask;
                            $tmp_other_seq &= $mask;
                                
                            my $n_differences = $tmp_seq =~ tr/[ACGT]//;
                            my $n_unknown_diff_other = $tmp_other_seq =~ tr/Z//;
                                
                            my $tot = $covered_length;
                                
                            # my $tot = $seq_length - $n_unknown - $n_unknown_diff_other;
                                
                            if ( $tot ) {
                                
                                my $differences =
                                $n_differences - $n_unknown_diff_other;
                            
                                $avg_chrom_length += $tot;
                                $avg_chrom_length_n++;
                                
                                $pi = 0 unless $pi;
                                
                                $pi +=
                                (

                                    $alleles_freq *
                                    $other_alleles_freq *
                                    ( $differences / $tot )

                                );
                                
                                # say "Differences: $differences " , $differences / $tot;
                            
                            }
                        
                        }
                        
                        # exit;
                        
                    }
                    
                    $pi *= 2;
                    $pi *= ( $nr_samples / ( $nr_samples -1 ) );
                    
                    # $avg_chrom_length /= $avg_chrom_length_n;
                    
                    if ( $nr_samples and $variable and $pi and $a_N and $a2 ) {
                    
                        # Pi: 0.04944 Theta: 2.16841 Tajima's D: -2.21623 SegSites: 18
                    
    # 					$tajimasd = tajima_D_counts_simple(
    # 						7 ,
    # 						18 ,
    # 						0.04944,
    # 						$a_N ,
    # 						$a2 ,
    # 					);
                    
                        $tajimasd = tajima_D_counts_simple(
                            $nr_samples ,
                            $variable ,
                            $pi * $covered_length ,
                            $a_N ,
                            $a2 ,
                        );
                        
                        if ( defined $min_dist_i ) {
                        
                            $dist_table->{ $min_dist_relation }->[ $min_dist_i ]->{ $group }->{ 'variable' } += $variable;
                            $dist_table->{ $min_dist_relation }->[ $min_dist_i ]->{ $group }->{ 'watterson' } += $watterson;
                            $dist_table->{ $min_dist_relation }->[ $min_dist_i ]->{ $group }->{ 'pi' } += $pi;
                            $dist_table->{ $min_dist_relation }->[ $min_dist_i ]->{ $group }->{ 'tajimasd' } += $tajimasd;
                            $dist_table->{ $min_dist_relation }->[ $min_dist_i ]->{ $group }->{ 'n' }++;

                        }
                    
                    }
                
                }
                
                print $out "\t$nr_samples\t$variable\t$watterson\t$pi\t$tajimasd";
                
            }
            
            say $out "";
        
        }
        
    }
    
    if ( keys %{ $dist_table } and $region eq "intergenic" or $region eq "any" ) {
    
        open (
            my $out ,
            ">" ,
            $opt->{ 'output' } . ".wattersons_theta_pi_tajimas_D.window_${window}.${region}.distance_from_gene.csv" ,
        );
        
        print $out "DISTANCE\tDISTANCE_MID";
        
        foreach my $group ( @groups ) {
        
            print $out "\t${group}_N\t${group}_VARIABLE\t${group}_THETA\t${group}_PI\t${group}_TD";
        
        }
        
        say $out "";
        
        if ( defined $dist_table->{ 'upstream' } ) {
        
            for ( my $i = $#{ $dist_table->{ 'upstream' } }; $i >= 0; $i-- ) {
            
                print $out "-" , $i * $window , "\t-" , $i * $window + $window / 2;
            
                foreach my $group ( @groups ) {
                
                    my $watterson = "";
                    my $pi = "";
                    my $tajimasd = "";
                    my $n = "";
                    my $variable = "";
                    
                    if ( defined $dist_table->{ 'upstream' }->[ $i ]->{ $group }->{ 'n' } ) {
                    
                        $n =
                        $dist_table->{ 'upstream' }->[ $i ]->{ $group }->{ 'n' };
                    
                        $variable =
                        $dist_table->{ 'upstream' }->[ $i ]->{ $group }->{ 'variable' };
                    
                        $watterson =
                        $dist_table->{ 'upstream' }->[ $i ]->{ $group }->{ 'watterson' };
                        
                        $pi =
                        $dist_table->{ 'upstream' }->[ $i ]->{ $group }->{ 'pi' };

                        $tajimasd =
                        $dist_table->{ 'upstream' }->[ $i ]->{ $group }->{ 'tajimasd' };
                        
                        $watterson /= $n;
                        $pi /= $n;
                        $tajimasd /= $n;
                        
                    }
                    
                    print $out "\t${n}\t$variable\t$watterson\t$pi\t$tajimasd";
                
                }
                
                say $out "";
            
            }
        
        }
        
        if ( defined $dist_table->{ 'downstream' } ) {
        
            for ( my $i = 0; $i <= $#{ $dist_table->{ 'downstream' } }; $i++ ) {
            
                print $out $i * $window , "\t" , $i * $window + $window / 2;
            
                foreach my $group ( @groups ) {
                
                    my $watterson = "";
                    my $pi = "";
                    my $tajimasd = "";
                    my $n = "";
                    my $variable = "";
                    
                    if ( defined $dist_table->{ 'downstream' }->[ $i ]->{ $group }->{ 'n' } ) {
                    
                        $n =
                        $dist_table->{ 'downstream' }->[ $i ]->{ $group }->{ 'n' };
                        
                        $variable =
                        $dist_table->{ 'downstream' }->[ $i ]->{ $group }->{ 'variable' };
                        
                        $watterson =
                        $dist_table->{ 'downstream' }->[ $i ]->{ $group }->{ 'watterson' };
                        
                        $pi =
                        $dist_table->{ 'downstream' }->[ $i ]->{ $group }->{ 'pi' };

                        $tajimasd =
                        $dist_table->{ 'downstream' }->[ $i ]->{ $group }->{ 'tajimasd' };
                        
                        $watterson /= $n;
                        $pi /= $n;
                        $tajimasd /= $n;
                        
                    }
                    
                    print $out "\t${n}\t$variable\t$watterson\t$pi\t$tajimasd";
                
                }
                
                say $out "";
            
            }
        
        }
    
    }

}

# The subroutine originates from BioPerl and was written by Jason Stajich

# https://metacpan.org/release/CJFIELDS/BioPerl-1.6.924/source/Bio/PopGen/Statistics.pm

# It has minor modifications and differs by reusing already calculated $a_N ($a1) and $a2

sub tajima_D_counts_simple {
    
    my ( $n , $seg_sites , $pi , $a1 , $a2 ) = @_;
        
    my $b1 = ( $n + 1 ) / ( 3 * ( $n - 1) );
    
    my $b2 = ( 2 * ( $n ** 2 + $n + 3 ) ) / ( ( 9 * $n ) * ( $n - 1 ) );
    
    my $c1 = $b1 - ( 1 / $a1 );
    
    my $c2 = $b2 - ( ( $n + 2 ) / ( $a1 * $n ) )+ ( $a2 / $a1 ** 2 );
    
    my $e1 = $c1 / $a1;
    my $e2 = $c2 / ( $a1**2 + $a2 );
    
    my $sum = ( $e1 * $seg_sites ) + ( ( $e2 * $seg_sites ) * ( $seg_sites - 1 ) );
    
    if ( $sum <= 0 ) {
    
		return "";
    
    }
    
    my $denom = sqrt ( $sum );

    return "" if $denom == 0;
    
    my $D = ( $pi - ( $seg_sites / $a1 ) ) / $denom;
    
    return $D;
    
}

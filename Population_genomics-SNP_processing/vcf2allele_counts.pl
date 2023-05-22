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
	'snp_table' => {} ,

	'groups' => {} ,
	
	'input' => [] ,
	
	'gt' => 0 ,
	'gl' => 0 ,
	'gp' => 0 ,
	'pl' => 0 ,
	
	'haploid' => 0 ,
	
};

GetOptions(
    $opt ,

    # Input
    # -----

    # VCF file
    'input|in|vcf=s{,}' ,

    # Groups
    'groups|g=s{,}' ,
    
    # Output
    # ------
    
    # Output name (string)
    'output|out|o=s' ,
    
    'verbose|v' ,
    
    'gt' ,
    'gl' ,
    'gp' ,
    'pl' ,

);

if ( @{ $opt->{ 'input' } } ) {

    parse_groups ( $opt );
    parse_vcf( $opt );

}

sub parse_groups {

	my ( $opt ) = @_;
	
	foreach my $group ( sort keys %{ $opt->{ 'groups' } } ) {
	
		my $group_file = $opt->{ 'groups' }->{ $group };
		
		if ( -e $group_file ) {
				
			open ( my $group_in , "<" , $group_file ) or die
			"Unable to open group file $group_file: $!";
		
			while (<$group_in>) {
			
				chomp;
				
				next if $_ =~ m/^\#/;
				
				$opt->{ 'sample_table' }->{ $_ } = $group;
				
				$opt->{ 'group_table' }->{ $group }->{ 'sample' } = $_;
			
			}
		
		}
	
	}

}

sub parse_vcf {
    
    my ( $opt ) = @_;
    
    my $group_table = $opt->{ 'group_table' };
    
    my @groups = sort keys %{ $opt->{ 'groups' } };
    
    my $output = "";
    
    if ( defined $opt->{ 'output' } ) {
    
		$output = $opt->{ 'output' };
    
    }
    
    else {
    
		$output = $opt->{ 'input' }->[ 0 ];
    
    }
    
    if ( $opt->{ 'gt' } ) {
    
        open_out_file( $opt , $output , "GT" , \@groups );
    
    }
    
    if ( $opt->{ 'gl' } ) {
    
        open_out_file( $opt , $output , "GL" , \@groups );
    
    }
    
    if ( $opt->{ 'gp' } ) {
    
        open_out_file( $opt , $output , "GP" , \@groups );
    
    }
    
    if ( $opt->{ 'pl' } ) {
    
        open_out_file( $opt , $output , "PL" , \@groups );
    
    }
    
    foreach my $vcf_file ( @{ $opt->{ 'input' } } ) {
    
		if ( -e $vcf_file ) {
			
			my $vcf_in;
			
			if ( $vcf_file =~ m/\.gz$/ ) {
			
				open ( $vcf_in , "gunzip -c $vcf_file |" ) or die
				"Unable to open VCF file $vcf_file: $!";
			
			}
			
			else {
			
				open ( $vcf_in , "<" , $vcf_file ) or die
				"Unable to open VCF file $vcf_file: $!";

			}
			
			# Get sample positions in the VCF file
			
            get_sample_positions( $opt , $group_table , $vcf_in );

            # Get the allele counts
            
            get_alleles( $opt , $group_table , $vcf_in );
			
		}
		
	}
    
}

sub open_out_file {

    my ( $opt , $output , $stat , $groups_ref ) = @_;
    
    open ( my $out , ">" , $output . ".allele_counts.${stat}.csv" ) or die 
    "Unable to create output: $!";
    
    $opt->{ 'out_table' }->{ $stat } = $out;

    print $out "CHROM\tPOS\tREF\tALT\tN_WITH_GT\t";
    
    if ( $opt->{ 'haploid' } ) {
    
        print $out join ( "\t" , "TOT" , @{ $groups_ref } );
    
    }
    
    else {
    
        print $out join ( "\t" , map { "$_\t$_" } ( "TOT" , @{ $groups_ref } ) );
    
    }
    
    say $out "";
    
    open ( my $f_out , ">" , $output . ".allele_freqs.${stat}.csv" ) or die 
    "Unable to create output: $!";
    
    $opt->{ 'out_freq_table' }->{ $stat } = $f_out;

    print $f_out "CHROM\tPOS\tREF\tALT\tN_WITH_GT\t";
    
    if ( $opt->{ 'haploid' } ) {
    
        print $f_out join ( "\t" , "TOT" , @{ $groups_ref } );
    
    }
    
    else {
    
        print $f_out join ( "\t" , map { "$_\t$_" } ( "TOT" , @{ $groups_ref } ) );
    
    }
    
    say $f_out "";

}

sub get_sample_positions {

    my ( $opt , $group_table , $vcf_in ) = @_;

    while (<$vcf_in>) {
            
        chomp;
        
        if ( $_ =~ m/^\#CHROM.+FORMAT\s+(.+)/ ) {
        
            my @samples = split ( /\s+/ , $1 );
            
            foreach my $group ( sort keys %{ $opt->{ 'groups' } } ) {
            
                $group_table->{ $group }->{ 'samples' } = [];
                $group_table->{ $group }->{ 'samples_i' } = [];
            
            }
            
            my $i = 0;
        
            foreach my $sample ( @samples ) {
            
                if ( defined $opt->{ 'sample_table' }->{ $sample } ) {
                
                    my $group = $opt->{ 'sample_table' }->{ $sample };
                    
                    push @{
                        $group_table->{ $group }->{ 'samples' }
                    } , $sample;
                    
                    push @{
                        $group_table->{ $group }->{ 'samples_i' }
                    } , $i;
                
                    say "Adding sample $sample [$i] to group $group ...";
                
                }
            
                $i++;
            
            }
            
            last;
        
        }
        
    }

}

sub get_alleles {

    my ( $opt , $group_table , $vcf_in ) = @_;

    my @groups = sort keys %{ $opt->{ 'groups' } };
    my $groups_ref = \@groups;
    
    my $get_gt = $opt->{ 'gt' };
    my $get_gl = $opt->{ 'gl' };
    my $get_gp = $opt->{ 'gp' };
    my $get_pl = $opt->{ 'pl' };

    my $format_table = {};
    
    my $positions = 0;
    my $last_chrom = "";
    my $chrom_i = 0;
    
    while (<$vcf_in>) {
    
        chomp;
        
        if (
        
            $_ =~ m/^(\S+)\s+(\d+)\s+\S+\s+(\S)\s+(\S)\s+\S+\s\S+\s\S+\;\S+\s+(GT\:\S+)\s+(.+)/ or
            $_ =~ m/^(\S+)\s+(\d+)\s+\S+\s+(\S)\s+(\S)\s+\S+\s\S+\s\S+\s+(GT)\s+(.+)/
            
        ) {
                
            my ( $chrom , $pos , $ref , $alt , $format , $data ) =
            ( $1 , $2 , $3 , $4 , $5 , $6 );
            
            unless ( keys %{ $format_table } ) {
            
                my $i = 0;
                
                foreach my $field ( split ( /\:/ , $format ) ) {
                
                    if ( $field eq "GT" ) {
                    
                        $format_table->{ 'GT' } = $i;
                    
                    }
                    
                    elsif ( $field eq "GL" ) {
                    
                        $format_table->{ 'GL' } = $i;
                    
                    }
                    
                    elsif ( $field eq "GP" ) {
                    
                        $format_table->{ 'GP' } = $i;
                    
                    }
                    
                    elsif ( $field eq "PL" ) {
                    
                        $format_table->{ 'PL' } = $i;
                    
                    }
                
                    $i++;
                
                }
            
            }
            
            my @genotypes = split ( /\s+/ , $data );
            
            my $genotypes_ref = \@genotypes;
            
            if ( $get_gt and defined $format_table->{ 'GT' } ) {

                my ( $n_with_gt , $counts_ref , $freqs_ref ) = get_counts_gt(
                    $opt ,
                    $group_table ,
                    $groups_ref ,
                    $genotypes_ref ,
                    $format_table->{ 'GT' } ,
                );
            
                my $out = $opt->{ 'out_table' }->{ 'GT' };
            
                say $out "$chrom\t$pos\t$ref\t$alt\t$n_with_gt\t" ,
                join ( "\t" , @{ $counts_ref } );
            
                my $f_out = $opt->{ 'out_freq_table' }->{ 'GT' };
                
                say $f_out "$chrom\t$pos\t$ref\t$alt\t$n_with_gt\t" ,
                join ( "\t" , @{ $freqs_ref } );
            
            }
            
            if ( $get_gl and defined $format_table->{ 'GL' } ) {
            
                 my ( $n_with_gt , $counts_ref , $freqs_ref ) = get_counts_gl(
                    $opt ,
                    $group_table ,
                    $groups_ref ,
                    $genotypes_ref ,
                    $format_table->{ 'GL' } ,
                );
                
                my $out = $opt->{ 'out_table' }->{ 'GL' };
            
                say $out "$chrom\t$pos\t$ref\t$alt\t$n_with_gt\t" ,
                join ( "\t" , @{ $counts_ref } );
            
                my $f_out = $opt->{ 'out_freq_table' }->{ 'GL' };
                
                say $f_out "$chrom\t$pos\t$ref\t$alt\t$n_with_gt\t" ,
                join ( "\t" , @{ $freqs_ref } );
            
            }
            
            if ( $get_gp and defined $format_table->{ 'GP' } ) {
            
                 my ( $n_with_gt , $counts_ref , $freqs_ref ) = get_counts_gp(
                    $opt ,
                    $group_table ,
                    $groups_ref ,
                    $genotypes_ref ,
                    $format_table->{ 'GP' } ,
                );
                
                my $out = $opt->{ 'out_table' }->{ 'GP' };
            
                say $out "$chrom\t$pos\t$ref\t$alt\t$n_with_gt\t" ,
                join ( "\t" , @{ $counts_ref } );
                
                my $f_out = $opt->{ 'out_freq_table' }->{ 'GP' };
                
                say $f_out "$chrom\t$pos\t$ref\t$alt\t$n_with_gt\t" ,
                join ( "\t" , @{ $freqs_ref } );
            
            }
            
            $positions++;
            
            if ( $chrom ne $last_chrom ) {
            
                $chrom_i++;
            
                say "$chrom ($chrom_i) $positions snps parsed so far ...";
            
            }
            
            $last_chrom = $chrom;
            
        }
    
    }
    
    say "$last_chrom ($chrom_i) $positions snps parsed in total";

}

sub get_counts_gt {

    my ( $opt , $group_table , $groups_ref , $genotypes_ref , $field_i ) = @_;

    my @counts;
    my @freqs;
    
    my $tot_n_0 = 0;
    my $tot_n_1 = 0;
    
    my $n_with_gt = 0;
    
    foreach my $group_name ( @{ $groups_ref } ) {
	
        my $n_0 = 0;
        my $n_1 = 0;
	
        foreach my $group_gt (
            
            @{ $genotypes_ref }[
                @{ $group_table->{ $group_name }->{ 'samples_i' } }
            ]
            
        ) {
        
            my @group_fields = split ( /\:/ , $group_gt );
            
            my $gt = $group_fields[ $field_i ];
        
            # Diploid
        
            if ( $gt =~ m/^(\d)\W(\d)$/ ) {
            
                my ( $a1 , $a2 ) = ( $1 , $2 );
                
                if ( $a1 == 0 ) {
                
                    $n_0++;
                
                }
            
                elsif ( $a1 == 1 ) {
                
                    $n_1++;
                
                }
                
                if ( $a2 == 0 ) {
                
                    $n_0++;
                
                }
            
                elsif ( $a2 == 1 ) {
                
                    $n_1++;
                
                }
                
                $n_with_gt++;
            
            }
            
            # Haploid
            
            elsif ( $group_gt =~ m/^(\d)$/ ) {
            
                my ( $a1 ) = ( $1 );
                
                if ( $a1 == 0 ) {
                
                    $n_0++;
                
                }
            
                elsif ( $a1 == 1 ) {
                
                    $n_1++;
                
                }
                
                $n_with_gt++;
            
            }
        
        }
        
        if ( $n_0 or $n_1 ) {
        
            push @counts , ( $n_0 , $n_1 );
        
            my $sum = $n_0 + $n_1;
            
            push @freqs , (
                sprintf( "%.5f", $n_0 / $sum ) ,
                sprintf( "%.5f", $n_1 / $sum )
            );
        
        }
        
        else {
        
            push @counts , ( '' , '' );
            push @freqs , ( '' , '' );
        
        }
	
        $tot_n_0 += $n_0;
        $tot_n_1 += $n_1;
	
    }
    
    if ( $tot_n_0 or $tot_n_1 ) {
    
        unshift @counts , ( $tot_n_0 , $tot_n_1 );
    
        my $sum = $tot_n_0 + $tot_n_1;
        
        unshift @freqs , (
            sprintf( "%.5f", $tot_n_0 / $sum ) ,
            sprintf( "%.5f", $tot_n_1 / $sum )
        );
    
    }
    
    else {
    
        unshift @counts , ( '' , '' );
        unshift @freqs , ( '' , '' );
    
    }
	
	return ( $n_with_gt , \@counts , \@freqs );

}

sub get_counts_gl {

    my ( $opt , $group_table , $groups_ref , $genotypes_ref , $field_i ) = @_;

    my @counts;
    my @freqs;
    
    my $tot_n_0 = 0;
    my $tot_n_1 = 0;
    
    my $n_with_gt = 0;
    
    foreach my $group_name ( @{ $groups_ref } ) {
	
        my $n_0 = 0;
        my $n_1 = 0;
	
        foreach my $group_gt (
            
            @{ $genotypes_ref }[
                @{ $group_table->{ $group_name }->{ 'samples_i' } }
            ]
            
        ) {
        
            my @group_fields = split ( /\:/ , $group_gt );
            
            my $gt = $group_fields[ $field_i ];
                
            my @gs;
        
            my $g_sum;
            
            # say "$group_gt\t$gt";
            
            if ( $gt =~ m/\d/ ) {
            
                foreach my $gl ( split ( /\,/ , $gt ) ) {
                
                    my $gl_trans = 10 ** $gl;
                
                    push @gs , $gl_trans;
                    
                    $g_sum += $gl_trans;
                    
                }
                
            }
            
            if ( $g_sum ) {
                
                foreach my $gl_trans ( @gs ) {
                
                    $gl_trans /= $g_sum;

                }
                
                $n_with_gt++;
            
                # Diploid
                
                if ( @gs == 3 ) {
                
                    $n_0 += 2 * $gs[ 0 ];
                    $n_0 += 1 * $gs[ 1 ];
                    $n_1 += 1 * $gs[ 1 ];
                    $n_1 += 2 * $gs[ 2 ];

                }
                
                # Haploid
                
                elsif ( @gs == 2 ) {
                
                    $n_0 += 1 * $gs[ 0 ];
                    $n_1 += 1 * $gs[ 1 ];

                }
            
            }
            
        }
            
        if ( $n_0 or $n_1 ) {
        
            $n_0 = sprintf( "%.5f", $n_0 );
            $n_1 = sprintf( "%.5f", $n_1 );
        
            push @counts , ( $n_0 , $n_1 );
        
            my $sum = $n_0 + $n_1;
            
            push @freqs , (
                sprintf( "%.5f", $n_0 / $sum ) ,
                sprintf( "%.5f", $n_1 / $sum )
            );
        
        }
        
        else {
        
            push @counts , ( '' , '' );
            push @freqs , ( '' , '' );
        
        }
        
        $tot_n_0 += $n_0;
        $tot_n_1 += $n_1;
	
    }
    
    if ( $tot_n_0 or $tot_n_1 ) {
    
        $tot_n_0 = sprintf( "%.5f", $tot_n_0 );
        $tot_n_1 = sprintf( "%.5f", $tot_n_1 );
    
        unshift @counts , ( $tot_n_0 , $tot_n_1 );
    
        my $sum = $tot_n_0 + $tot_n_1;
        
        unshift @freqs , (
            sprintf( "%.5f", $tot_n_0 / $sum ) ,
            sprintf( "%.5f", $tot_n_1 / $sum )
        );
    
    }
    
    else {
    
        unshift @counts , ( '' , '' );
        unshift @freqs , ( '' , '' );
    
    }
	
	return ( $n_with_gt , \@counts , \@freqs );

}

sub get_counts_gp {

    my ( $opt , $group_table , $groups_ref , $genotypes_ref , $field_i ) = @_;

    my @counts;
    my @freqs;
    
    my $tot_n_0 = 0;
    my $tot_n_1 = 0;
    
    my $n_with_gt = 0;
    
    foreach my $group_name ( @{ $groups_ref } ) {
	
        my $n_0 = 0;
        my $n_1 = 0;
	
        foreach my $group_gt (
            
            @{ $genotypes_ref }[
                @{ $group_table->{ $group_name }->{ 'samples_i' } }
            ]
            
        ) {
        
            my @group_fields = split ( /\:/ , $group_gt );
            
            my $gt = $group_fields[ $field_i ];
            
            my @gps = split ( /\,/ , $gt );
        
            my $g_sum;
            
            map { $g_sum += $_ } @gps;
                        
            if ( $g_sum =~ m/\d/ ) {
                
                $n_with_gt++;
            
                # Diploid
                
                if ( @gps == 3 ) {
                
                    $n_0 += 2 * $gps[ 0 ];
                    $n_0 += 1 * $gps[ 1 ];
                    $n_1 += 1 * $gps[ 1 ];
                    $n_1 += 2 * $gps[ 2 ];

                }
                
                # Haploid
                
                elsif ( @gps == 2 ) {
                
                    $n_0 += 1 * $gps[ 0 ];
                    $n_1 += 1 * $gps[ 1 ];

                }
            
            }
            
        }
            
        if ( $n_0 or $n_1 ) {
        
            $n_0 = sprintf( "%.5f", $n_0 );
            $n_1 = sprintf( "%.5f", $n_1 );
        
            push @counts , ( $n_0 , $n_1 );
            
            my $sum = $n_0 + $n_1;
            
            push @freqs , (
                sprintf( "%.5f", $n_0 / $sum ) ,
                sprintf( "%.5f", $n_1 / $sum )
            );
        
        }
        
        else {
        
            push @counts , ( '' , '' );            
            push @freqs , ( '' , '' );
        
        }
        
        $tot_n_0 += $n_0;
        $tot_n_1 += $n_1;
	
    }
    
    if ( $tot_n_0 or $tot_n_1 ) {
    
        $tot_n_0 = sprintf( "%.5f", $tot_n_0 );
        $tot_n_1 = sprintf( "%.5f", $tot_n_1 );
    
        unshift @counts , ( $tot_n_0 , $tot_n_1 );
    
        my $sum = $tot_n_0 + $tot_n_1;
        
        unshift @freqs , (
            sprintf( "%.5f", $tot_n_0 / $sum ) ,
            sprintf( "%.5f", $tot_n_1 / $sum )
        );
    
    }
    
    else {
    
        unshift @counts , ( '' , '' );
        unshift @freqs , ( '' , '' );
    
    }
	
	return ( $n_with_gt , \@counts , \@freqs );

}

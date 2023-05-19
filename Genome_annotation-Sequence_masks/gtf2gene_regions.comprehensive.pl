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
    
    # Sequence length table
    
    'length' => [] ,
    
    # Sequence length table (bed format)
    
    'bed' => [] ,
    
    # GTF file with gene coordinates
    
    'gtf' => [] ,
    
    # Accessibility mask in FASTA format
    
    'coverage' => [] ,
    
};

# Region/symbols

my $region_table = {

    'uncovered' => 0 ,
    'inaccessible' => 0 ,
    
    'intergenic' => 1 ,
    
    'gene' => 2 ,
    'transcript' => 2 ,
    
    '3_utr' => 3 ,
    'three_prime_utr' => 3 ,
    
    'exon' => 4 ,
    
    '5_utr' => 5 ,
    'five_prime_utr' => 5 ,
    
    'cds' => 6 ,
    'start_codon' => 6 ,
    'stop_codon' => 6 ,

};

my $symbol_table = {

    '0' => 'inaccessible' ,
    
    '1' => 'intergenic' ,
    
    '2' => 'intron' ,
    
    '3' => 'three_prime_utr' ,
    
    '4' => 'exon' ,
    
    '5' => 'five_prime_utr' ,
    
    '6' => 'cds' ,

};

$opt->{ 'app' }->{ 'region_table' } = $region_table;
$opt->{ 'app' }->{ 'symbol_table' } = $symbol_table;

GetOptions(
    $opt ,

    'length|l=s{,}' ,
    'bed|b=s{,}' ,
    
    'gtf|gff|g=s{,}' ,
    
    'coverage|cov|c=s{,}' ,
    
    'output|out|o=s' ,
   
    # Booleans
    
    'verbose|v' ,
    'help|h' ,

);

if ( @{ $opt->{ 'gtf' } } ) {

    unless ( defined $opt->{ 'output' } ) {

        $opt->{ 'output' } = $opt->{ 'gtf' }->[ 0 ];

    }

    # Get the lengths files

    get_lengths( $opt );

    # Parse the GTF/GFF files

    read_genes( $opt );

    # Print data

    print_data( $opt , "complete" );
    
    # Mask data
    
    mask_data( $opt );
    
    # Print data

    print_data( $opt , "masked_inaccessible" );

}

sub get_lengths {

    my ( $opt ) = @_;

    my $seq_table = {};
    my @sequences;
    
    $opt->{ 'app' }->{ 'seq_table' } = $seq_table;
    $opt->{ 'app' }->{ 'sequences' } = \@sequences;
    
    foreach my $length_file ( @{ $opt->{ 'length' } } ) {

        open ( my $length_in , "<" , $length_file ) or die "$!";

        say "Reading lengths from $length_file ...";

        my $header = <$length_in>;

        while (<$length_in>) {

            chomp;
            
            if ( $_ =~ m/^\S+\t([^\t]+)\t(\d+)/ ) {
            
                my $sequence = $1;
            
                unless ( defined $seq_table->{ $sequence } ) {
            
                    $seq_table->{ $sequence }->{ 'length' } = $2;
                    $seq_table->{ $sequence }->{ 'seq' } = 1 x $2;
                
                    push @sequences , $sequence;
            
                }
            
            }

        }
        
    }
    
    foreach my $length_file ( @{ $opt->{ 'bed' } } ) {

        open ( my $length_in , "<" , $length_file ) or die "$!";

        say "Reading lengths from $length_file ...";

        my $header = <$length_in>;

        while (<$length_in>) {

            chomp;
            
            next if $_ =~ m/^\#/;
            
            if ( $_ =~ m/^([^\t]+)\t[^\t]+\t(\d+)/ ) {
            
                my ( $sequence , $length ) = ( $1 , $2 + 1 );
            
                unless ( defined $seq_table->{ $sequence } ) {
            
                    $seq_table->{ $sequence }->{ 'length' } = $length;
                    $seq_table->{ $sequence }->{ 'seq' } = 1 x $length;
                
                    push @sequences , $sequence;
                
                }
            
            }

        }
        
    }

}

sub read_genes {

    my ( $opt ) = @_;

    my $seq_table = $opt->{ 'app' }->{ 'seq_table' };
    
    foreach my $gtf_file ( @{ $opt->{ 'gtf' } } ) {

        say "Reading gene coordinates from $gtf_file ...";

        open ( GTF_IN , "<" , $gtf_file ) or die "$!";

        while (<GTF_IN>) {

            chomp;
            
            next if $_ =~ m/^\#/;
            next unless $_ =~ m/^\w/;
            
            my @data = split ( /\s/ , $_ );
            
            my ( $sequence , $region , $start , $stop ) =
            @data[ 0 , 2 , 3 , 4 ];
            
            $region = lc $region;
            
            if (
            
                defined $region_table->{ $region } and
                defined $seq_table->{ $sequence }->{ 'seq' }
            
            ) {
            
                my $symbol = $region_table->{ $region };

                # Swap coordinates around if they are in opposite order
            
                if ( $start > $stop ) {
                
                    ( $start , $stop ) = ( $stop , $start );
                
                }
                
                # Convert to 0-based coordinates
                
                $start--;
                $stop--;
            
                my $region_length = $stop - $start + 1;
            
                $seq_table->{ $sequence }->{ 'region' }->{ $region }->{ 'n_orig' }++;
                $seq_table->{ $sequence }->{ 'region' }->{ $region }->{ 'length_orig' } += $region_length;
            
                my $current_region = substr (
                    $seq_table->{ $sequence }->{ 'seq' } ,
                    $start ,
                    $region_length ,
                );
                
                my $n_intergenic = $current_region =~ tr/1//;
                
                if ( $n_intergenic == $region_length ) {
            
                    if ( $opt->{ 'verbose' } ) {
            
                        say "Fast insertion of $region ($symbol) sequence along $sequence ...";
                        say "\t$start - $stop: $region_length vs $n_intergenic";
                    
                    }
                    
                    substr (
                        $seq_table->{ $sequence }->{ 'seq' } ,
                        $start ,
                        $region_length ,
                        $symbol x $region_length ,
                    );
                
                }
                
                else {
                
                    if ( $opt->{ 'verbose' } ) {
                
                        say "Careful insertion of $region ($symbol) sequence along $sequence ...";
                        say "\t$start - $stop: $region_length vs $n_intergenic";
                
                    }
                
                    for ( my $i = $start; $i <= $stop; $i++ ) {
                    
                        my $base = substr (
                            $seq_table->{ $sequence }->{ 'seq' } ,
                            $i ,
                            1 ,
                        );
                        
                        if ( $symbol == 4 and $base == 3 ) {
                        
                            next;
                        
                        }
                        
                        elsif ( $base < $symbol ) {
                        
                            substr (
                                $seq_table->{ $sequence }->{ 'seq' } ,
                                $i ,
                                1 ,
                                $symbol ,
                            );
                        
                        }
                        
                        elsif ( $symbol == 3 and $base == 4 ) {
                        
                            substr (
                                $seq_table->{ $sequence }->{ 'seq' } ,
                                $i ,
                                1 ,
                                $symbol ,
                            );
                        
                        }
                    
                    }
                
                }

            }

        }
    
    }

}

sub print_data {

    my ( $opt , $tag ) = @_;

    my $seq_table = $opt->{ 'app' }->{ 'seq_table' };
        
    my $region_table = $opt->{ 'app' }->{ 'region_table' };
    my $symbol_table = $opt->{ 'app' }->{ 'symbol_table' };
    
    my @sorted_symbols = sort keys %{ $symbol_table };
    
    my @sequences = @{ $opt->{ 'app' }->{ 'sequences' } };
    
    my $table_out_file = $opt->{ 'output' } . ".regions.${tag}.table.csv";
    
    open ( my $t_out , ">" , $table_out_file ) or die "";
    
    print $t_out "SEQUENCE\tLENGTH\tGENES";
    
    foreach my $symbol ( @sorted_symbols ) {
    
        my $region = $symbol_table->{ $symbol };
        
        print $t_out
        "\t${region}_LENGTH\t${region}_PROP\t${region}_N\t${region}_AVG";
    
    }
    
    say $t_out "";
    
    my $fasta_out_file = $opt->{ 'output' } . ".regions.${tag}.fasta";
    
    open ( my $f_out , ">" , $fasta_out_file ) or die "$!";
    
    my $stat_table = {};
    
    foreach my $sequence ( @sequences ) {
    
        my $seq = $seq_table->{ $sequence }->{ 'seq' };
        
        my $seq_length = length $seq;
        
        say $f_out ">$sequence\n" , $seq;
        
        my $n_genes = 0;
        
        if (
        
            defined
            $seq_table->{ $sequence }->{ 'region' }->{ 'gene' }->{ 'n_orig' }
        
        ) {
        
            $n_genes =
            $seq_table->{ $sequence }->{ 'region' }->{ 'gene' }->{ 'n_orig' };
        
        }
        
        $stat_table->{ 'genes' } += $n_genes;
        $stat_table->{ 'length' } += $seq_length;
        
        print $t_out "$sequence\t$seq_length\t$n_genes";
        
        foreach my $symbol ( @sorted_symbols ) {
        
            my $region = $symbol_table->{ $symbol };
        
            my $symbol_length = 0;
            my $symbol_n = 0;
        
            while ( $seq =~ m/(${symbol}+)/g ) {
        
                $symbol_length += length $1;
                $symbol_n++;
        
            }
            
            my $avg = "";
            
            if ( $symbol_n ) {
            
                $avg = ${symbol_length} / $symbol_n;
            
            }
            
            $stat_table->{ 'symbol' }->{ $symbol }->{ 'n' } += $symbol_n;
            $stat_table->{ 'symbol' }->{ $symbol }->{ 'length' } += $symbol_length;
            
            my $symbol_prop = $symbol_length / $seq_length;
            
            print $t_out "\t${symbol_length}\t$symbol_prop\t$symbol_n\t$avg";
        
        }
        
        say $t_out "";
    
    }
    
    my $global_out_file =
    $opt->{ 'output' } . ".regions.${tag}.table.GLOBAL.csv";
    
    open ( my $g_out , ">" , $global_out_file ) or die "";
    
    print $g_out "SEQUENCE\tLENGTH\tGENES";
    
    foreach my $symbol ( @sorted_symbols ) {
    
        my $region = $symbol_table->{ $symbol };
        
        print $g_out
        "\t${region}_LENGTH\t${region}_PROP\t${region}_N\t${region}_AVG";
    
    }
    
    say $g_out "";
    
    $stat_table->{ 'genes' } = 0 unless defined
    $stat_table->{ 'genes' };
    
    $stat_table->{ 'length' } = 0 unless defined
    $stat_table->{ 'length' };
    
    print $g_out
        "ALL\t" ,
        $stat_table->{ 'length' } , "\t" ,
        $stat_table->{ 'genes' };
    
    foreach my $symbol ( @sorted_symbols ) {
    
        my $region = $symbol_table->{ $symbol };
        
        my $symbol_n =
        $stat_table->{ 'symbol' }->{ $symbol }->{ 'n' };
        
        my $symbol_length =
        $stat_table->{ 'symbol' }->{ $symbol }->{ 'length' };
        
        my $avg = "";
        
        if ( $symbol_n ) {
        
            $avg = ${symbol_length} / $symbol_n;
        
        }
        
        my $symbol_prop = $symbol_length / $stat_table->{ 'length' };
        
        print $g_out "\t${symbol_length}\t$symbol_prop\t$symbol_n\t$avg";
        
    }
        
    say $g_out "";
    
}

sub mask_data {

    my ( $opt ) = @_;

    my $seq_table = $opt->{ 'app' }->{ 'seq_table' };
    
    foreach my $coverage_file ( @{ $opt->{ 'coverage' } } ) {
        
        my $coverage_table = {};
        
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
            
                if ( defined $seq_table->{ $1 } ) {
                
                    $tmp_header = $1;
                
                }
            
            }
            
            elsif ( defined $tmp_header ) {
            
                $coverage_table->{ $tmp_header } .= $_;
            
            }
            
        }
        
        foreach my $sequence ( keys %{ $coverage_table } ) {
        
            my $cov_seq = $coverage_table->{ $sequence };
        
            while ( $cov_seq =~ m/(0+)/g ) {
            
                my $start = $-[ 0 ];
                my $stop = $+[ 0 ];
            
                $stop--;
            
                my $length = $stop - $start + 1;
            
                substr (
                    $seq_table->{ $sequence }->{ 'seq' } ,
                    $start ,
                    $length ,
                    0 x $length ,
                );
            
            }
        
        }
        
    }
    
}

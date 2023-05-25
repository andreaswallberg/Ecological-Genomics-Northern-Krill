#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

use Getopt::Long;
use Data::Dumper;

# Option table and defaults

my $opt = {

	'fst' => [] ,
	'ad' => [] ,
	'xp' => [] ,
	'div' => [] ,
	'div_com' => [] ,
	'div_region' => [] ,
	'z' => [] ,

	'types' => [ 'fst' , 'ad' , 'xp' , 'div' , 'div_com' , 'div_region' , 'z' ] ,
	
	'out' => "consolidated.csv" ,
	
};

GetOptions(
    $opt ,

    # Input
    # -----

    'chromosome_list=s' ,
    
    'fst|fst_windows=s{,}' ,
    'ad|ancestral_derived_windows|ancestral_derived=s{,}' ,
    'xp|xp_windows=s{,}' ,
    'div|diversity_windows=s{,}' ,
    'div_com|diversity_comparative_windows=s{,}' ,
    'div_region|diversity_region_windows=s{,}' ,
    'z|z_windows=s{,}' ,
    
    # Output
    # ------
    
    # Output name (string)
    'out|output|o=s' ,
    
    'verbose|v' ,

);

my $data_table = {};
$opt->{ 'app' }->{ 'data_table' } = $data_table;

read_chromosome_list();
read_all_windows();
print_output();

sub read_chromosome_list {

    my $chromosome_list_file = $opt->{ 'chromosome_list' };
    
    my ( $in ) = get_filehandle( $chromosome_list_file , 'chromosome_list' );

    while (<$in>) {
    
        chomp;
    
        if ( $_ =~ m/^(\S+)/ ) {
        
            my $chrom = $1;
            
            unless ( defined $data_table->{ 'universal' }->{ $chrom } ) {
    
                push @{
                    $data_table->{ 'universal' }->{ 'chromosomes' }
                } , $chrom;
                
                $data_table->{ 'universal' }->{ $chrom } = {};
    
            }
        
        }
    
    }
    
}

sub read_all_windows {

    foreach my $type ( 'fst' , 'z' , 'div_region' ) {

        read_windows(

            {
                'type' => $type ,
                'get_header' => 1 ,
                'start_field' => 2 ,
            } ,
            
        );
    
    }

    foreach my $type ( 'ad' , 'xp' , 'div' , 'div_com' ) {
        
        read_windows(
        
            {
                'type' => $type ,
                'get_header' => 1 ,
                'start_field' => 3 ,
            } ,
        
        );

    }

}

sub read_windows {

    my ( $config_table ) = @_;
    
    my $type = $config_table->{ 'type' };
    
    return unless defined $opt->{ $type };
    
    my $start_field = $config_table->{ 'start_field' };
    
    my $data_table = $opt->{ 'app' }->{ 'data_table' };

    foreach my $file ( @{ $opt->{ $type } } ) {
    
        my ( $in ) = get_filehandle( $file , $type );
        
        my $headers_ref;
        
        if ( $config_table->{ 'get_header' } ) {
        
            ( $headers_ref ) = get_headers ( $in , $start_field ); 
    
        }
    
        $data_table->{ $type }->{ 'header' } = $headers_ref;
    
        if ( $start_field == 2 ) {
    
            while (<$in>) {
            
                chomp;
                
                if ( $_ =~ m/^(\S+)\s(\d+)\s(.+)/ ) {
                
                    my ( $chrom , $start , $data ) = ( $1 , $2 , $3 );
                    
                    register_data( $data_table , $type , $chrom , $start , $data );
                
                }
            
            }
        
        }
        
        elsif ( $start_field == 3 ) {
    
            while (<$in>) {
            
                chomp;
                
                if ( $_ =~ m/^(\S+)\s(\d+)\s\d+\s(.+)/ ) {
                
                    my ( $chrom , $start , $data ) = ( $1 , $2 , $3 );
                    
                    register_data( $data_table , $type , $chrom , $start , $data );

                }
            
            }
        
        }
    
    }

}

sub get_filehandle {

    my ( $file , $type ) = @_;
    
    my $in;
    
    if ( $file =~ m/\.gz$/ ) {
    
        open ( $in , "zcat $file |" ) or die "Unable to read $file: $!";
    
    }
    
    else {
    
        open ( $in , "<" , $file ) or die "Unable to read $file: $!";
    
    }
    
    say "Reading $type $file ...";
    
    return ( $in );
    
}

sub get_headers {

    my ( $in , $cut ) = @_;

    my $header = <$in>;
    chomp $header;

    my @tmp = split ( /\t/ , $header );
    
    my @headers = @tmp[ $cut .. $#tmp ];
    
    return ( \@headers );
    
}

sub register_data {

    my ( $data_table , $type , $chrom , $start , $data ) = @_;

    if ( defined $data_table->{ 'universal' }->{ $chrom } ) {
    
        if ( $start =~ m/0$/ ) {
        
            $start += 1;
        
        }
        
        unless ( defined $data_table->{ 'universal' }->{ $chrom }->{ $start } ) {
        
            $data_table->{ 'universal' }->{ $chrom }->{ $start } = 1;
        
        }
        
        if ( $data =~ m/\S/ ) {
        
            $data_table->{ $type }->{ 'data' }->{ $chrom }->{ $start } = $data;

        }
    
    }
    
}

sub print_output {
    
    my $data_table = $opt->{ 'app' }->{ 'data_table' };
    
    my $out_file = $opt->{ 'out' };
    
    open ( my $out , ">" , $out_file ) or die "$!";
    
    print $out "CHROM\tSTART";
    
    my @types = grep { defined $data_table->{ $_ }->{ 'data' } } @{ $opt->{ 'types' } };
    
    foreach my $type ( @types ) {
    
        print $out "\t$type\t" , join ( "\t" , @{ $data_table->{ $type }->{ 'header' } } );
    
    }
    
    say $out "";
    
    my @chromosomes = grep {
        keys %{ $data_table->{ 'universal' }->{ $_ } }
    } @{ $data_table->{ 'universal' }->{ 'chromosomes' } };
    
    foreach my $chrom ( @chromosomes ) {
    
        my @starts = sort { $a <=> $b } keys %{ $data_table->{ 'universal' }->{ $chrom } };
        
        foreach my $start ( @starts ) {
        
            print $out "$chrom\t$start";
            
            foreach my $type ( @types ) {
            
                my $n_fields = @{ $data_table->{ $type }->{ 'header' } };
            
                my $string = "\t" x $n_fields;
                
                if ( defined $data_table->{ $type }->{ 'data' }->{ $chrom }->{ $start } ) {
                
                    $string = "\t" . $data_table->{ $type }->{ 'data' }->{ $chrom }->{ $start };
                
                }
                
                print $out "\t${type}${string}";
            
            }
            
            say $out "";
    
        }
    
    }

    close $out;
    
}

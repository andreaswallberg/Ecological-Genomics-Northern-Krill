#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $pol_table = {};

foreach my $gt_file ( @ARGV ) {

    open ( my $in , "<" , $gt_file ) or die "$!";

    say STDERR $gt_file;
    
    my $header = <$in>;
    
    while (<$in>) {

        chomp;
        
        my @data = split ( /\t/ , $_ );
        
        my $n_pol;
        
        # Skip positions where the minor allele is a singleton
        
        if ( $data[ 5 ] < 2 or $data[ 6 ] < 2 ) {
        
            next;
        
        }
        
        for ( my $i = 7; $i < @data; $i += 2 ) {
        
            # say "$data[$i]\t$data[$i+1]";
        
            if ( $data[ $i ] == 0 or $data[ $i + 1 ] == 0 ) {
            
                next;
            
            }
            
            else {
            
                $n_pol++;
            
            }
        
        }
        
        next unless defined $n_pol;
        
        $pol_table->{ 'n_pol' }->{ $n_pol }++;
        $pol_table->{ 'n' }++;
        
        # last if $pol_table->{ 'n' } >= 100_000;

    }

}

my @n_pols = sort { $a <=> $b } keys %{ $pol_table->{ 'n_pol' } };

foreach my $n_pol ( @n_pols ) {

    say "$n_pol\t" , $pol_table->{ 'n_pol' }->{ $n_pol } , "\t" ,
    $pol_table->{ 'n_pol' }->{ $n_pol } / $pol_table->{ 'n' };

}

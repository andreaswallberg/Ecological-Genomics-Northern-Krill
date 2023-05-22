#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $locus_table = {};

while (<>) {

    chomp;
    
    # XLOC_000022.1.p1	XP_027237621.1	77.9	810	172	4	11	820	3	805	0.0e+00	1239.9
    # XLOC_000022.1.p1	XP_037775203.1	76.7	810	182	4	11	820	3	805	0.0e+00	1220.7
    
    if ( $_ =~ m/^(\S+)\s.+\s(\S+)$/ ) {
    
        my ( $transcript , $score ) = ( $1 , $2 );
    
		# Assumes all accessions belong to the same locus are indexed following a dot.
		# What comes before the dot is the locus tag.
    
        if ( $transcript =~ m/^([^\.]+)\./ ) {
        
            my $locus = $1;
            
            if ( defined $locus_table->{ $locus } and $score > $locus_table->{ $locus }->{ 'score' } ) {
            
                $locus_table->{ $locus }->{ 'score' } = $score;
                $locus_table->{ $locus }->{ 'transcript' } = $transcript;
                
            }
            
            elsif ( not defined $locus_table->{ $locus } ) {
            
                $locus_table->{ $locus }->{ 'score' } = $score;
                $locus_table->{ $locus }->{ 'transcript' } = $transcript;
            
            }
        
        }
    
    }

}

foreach my $locus ( sort keys %{ $locus_table } ) {

    say "$locus\t" , $locus_table->{ $locus }->{ 'transcript' } , "\t" ,
    $locus_table->{ $locus }->{ 'score' };

}

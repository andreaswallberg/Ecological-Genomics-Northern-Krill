#!/usr/bin/perl

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $og_table = {};
my @ogs;

my $i = 1;

while (<>) {

    chomp;
    
    my ( $og , $genes_string ) = ( "SWIFTORTHO_" . $i , $_ );

    push @ogs , $og;

    my @data = split ( /\s/ , $genes_string );
    
    foreach my $gene_string ( @data ) {
    
        if ( $gene_string =~ m/([^\|]+)\|(\S+)/ ) {
    
            my ( $species , $gene ) = ( $1 , $2 );
    
            $og_table->{ 'og' }->{ $og }->{ $species }++;
            
            if ( $gene =~ m/mnor/ ) {
            
                unless ( defined $og_table->{ 'og_rep' }->{ $og } ) {
                
                    $og_table->{ 'og_rep' }->{ $og } = $gene;
                
                }
            
            }
            
            unless ( defined $og_table->{ 'species' }->{ $species } ) {
            
                $og_table->{ 'species' }->{ $species } = 1;
            
            }
            
        }
    
    }
    
    $i++;

}

delete $og_table->{ 'species' }->{ 'mnor' };

my @sorted_species = ( sort keys %{ $og_table->{ 'species' } } );

say "Desc\tFamily ID\tmnor\t" , join ( "\t" , @sorted_species );

foreach my $og ( @ogs ) {

    my $desc = "NA";
    
    unless ( defined $og_table->{ 'og' }->{ $og }->{ 'mnor' } ) {
    
        $og_table->{ 'og' }->{ $og }->{ 'mnor' } = 0;
    
    }
    
    my $n_mnor = $og_table->{ 'og' }->{ $og }->{ 'mnor' };
    
    my $tot = 0;
    
    my @vals;
    
    foreach my $species ( @sorted_species ) {
    
        $og_table->{ 'og' }->{ $og }->{ $species } = 0 unless defined
        $og_table->{ 'og' }->{ $og }->{ $species };
        
        $tot += $og_table->{ 'og' }->{ $og }->{ $species };
        
        push @vals , $og_table->{ 'og' }->{ $og }->{ $species };
    
    }
    
    my $avg = median( @vals );

    if ( defined $og_table->{ 'og_rep' }->{ $og } ) {
    
        $desc = $og_table->{ 'og_rep' }->{ $og };
    
    }
    
    my $ratio = "";
    
    if ( $avg ) {
    
        $ratio = sprintf("%.2f" , $n_mnor / $avg );
    
    }
    
    $avg = sprintf("%.2f" , $avg );
    
    $desc .= "|$avg|$ratio";
    
    print "$desc\t$og\t$n_mnor";
    
    foreach my $species ( @sorted_species ) {
    
        print "\t" , $og_table->{ 'og' }->{ $og }->{ $species };
        
    }
    
    say "";
    
}

sub median {

    my @vals = sort {$a <=> $b} @_;
    
    my $len = @vals;
    
    if ( $len % 2 ) {
        return $vals[int($len/2)];
    }
    
    else {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
    
}

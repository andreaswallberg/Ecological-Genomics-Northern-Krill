#!/usr/bin/perl -w

use strict;
use warnings;
use 5.010;

my $seq_table = {};

my @seqs;

my $out_file = shift @ARGV;

open ( my $out , ">" , $out_file ) or die "$!";

foreach my $depth_file ( @ARGV ) {

    my $in;

    my $i = 1;
    
    if ( $depth_file =~ m/\.gz$/ ) {
    
        open ( $in , "zcat $depth_file |" ) or die "$!";
    
    }
    
    else {

        open ( $in , "<" , $depth_file ) or die "$!";

    }
    
    while (<$in>) {
        
        if ( $_ =~ m/^(\S+)\s+(.+)/ ) {
        
            unless ( defined $seq_table->{ $1 } ) {
        
                push @seqs, $1;
        
            }
            
            say "$depth_file\t$i\t$1";
            
            $seq_table->{ $1 } .= $_;
        
        }
        
        $i++;

    }
    
}

foreach my $seq ( @seqs ) {

    print $out $seq_table->{ $seq };

}

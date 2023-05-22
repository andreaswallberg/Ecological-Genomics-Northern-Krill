#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $seq_table = {};

my $tmp_header;

my @headers;

while (<>) {

    chomp;
    
    if ( $_ =~ m/^\>(.+)/ ) {
    
        $tmp_header = $1;
    
    }
    
    else {
    
        unless ( defined $seq_table->{ $tmp_header } ) {
        
            push @headers , $tmp_header;
        
        }
    
        $seq_table->{ $tmp_header } .= $_;
    
    }

}

foreach my $tmp_header ( @headers ) {

    my $seq = $seq_table->{ $tmp_header };
    
    my $seq_length = length $seq;
    
    my $n = $seq =~ tr/Nn//;
    
    say STDERR "$tmp_header\t$seq_length\t$n";
    
    if ( ( $n / $seq_length ) < 0.8 ) {
    
        say $tmp_header;
    
    }

}

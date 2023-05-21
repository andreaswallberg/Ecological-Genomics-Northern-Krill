#!/usr/bin/perl -w

use strict;
use warnings;
use 5.010;

my $out_file = shift @ARGV;

open ( my $out , ">" , $out_file ) or die "$!";

foreach my $depth_file ( @ARGV ) {

    open ( my $in , "<" , $depth_file ) or die "$!";

    my $last_seq = "";
    
    my @depths;
    my $sample = 0;
    
    print $out "#SAMPLE";
    
    while (<$in>) {
        
        chomp;
        
        if ( $_ =~ m/^\#SAMPLE\s(.+)/ ) {
        
            print $out " $1"; 
        
        }
        
        elsif ( $_ =~ m/^(\S+)\s+(.+)/ ) {
        
            if ( $1 ne $last_seq ) {
            
                if ( @depths ) {
            
                    say $out $last_seq , "\t" , join ( "\t" , @depths );
                
                    @depths = ();
                    
                    $sample = 0;
                
                }
                
                else {
                
                    say $out "";
                
                }
                
                say "Computing for $1 ...";
            
            }
            
            $sample++;
            
            my @tmp_depths = split ( /\t/ , $2 );
            
            if ( @depths and @tmp_depths != @depths ) {
            
                die "Unequal number of depth data points at $1 for sample $sample!";
            
            }
            
            for ( my $i = 0; $i < @tmp_depths; $i++ ) {
            
                $depths[ $i ] += $tmp_depths[ $i ];
            
            }
            
            $last_seq = $1;
        
        }

    }
    
    if ( @depths ) {

        say $out $last_seq , "\t" , join ( "\t" , @depths );

        @depths = ();

    }
    
}


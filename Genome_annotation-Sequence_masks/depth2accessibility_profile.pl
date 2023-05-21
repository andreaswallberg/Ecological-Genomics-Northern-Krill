#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $min_depth = shift @ARGV;
my $max_depth = shift @ARGV;

my $out_file = shift @ARGV;

my @all_depths;
my $all_n = 0;

open ( my $d_fa_out , ">" , $out_file . ".depth_profile.fa" ) or die "$!";
open ( my $d_out , ">" , $out_file . ".depth_profile.csv" ) or die "$!";
open ( my $d_table_out , ">" , $out_file . ".depth_profile.table.csv" ) or die "$!";

say $d_table_out "SEQ\tLENGTH\tLENGTH_COVERED\tPROP_COVERED";

my $global_tmp_pos = 0;

foreach my $depth_file ( @ARGV ) {
    
    my $in;
    
    if ( $depth_file =~ m/\.gz$/ ) {
    
        open ( $in , "zcat $depth_file|" ) or die "$!";
    
    }
    
    else {
    
        open ( $in , "<" , $depth_file ) or die "$!";
    
    }
    
    while (<$in>) {
        
        chomp;
        
        if ( $_ =~ m/^\#/ ) {
        
            say $d_out "$_";
        
        }
        
        elsif ( $_ =~ m/^(\S+)\s(.+)/ ) {

            my $i = 1;
            
            my @depths = split ( /\t/ , $2 );
            
            my $depth_profile = "";
            
            my $length = @depths;
            my $lengths_covered = 0;
            
            foreach my $depth ( @depths ) {
            
                if ( $depth >= $min_depth and $depth <= $max_depth ) {
            
                    $lengths_covered++;
                    $depth_profile .= "1";
                
                }
                
                else {
                
                    $depth_profile .= "0";
                
                }
                
                $i++;
                
            }
            
            say $d_fa_out ">$1\n" , $depth_profile;
            
            say $d_out "$1\t" , join ( "\t" , split ( // , $depth_profile ) );
            
            say $d_table_out "$1\t$length\t$lengths_covered\t" , $lengths_covered / $length;
        
        }
        
    }

}

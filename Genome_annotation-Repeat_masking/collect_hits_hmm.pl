#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $repeat_table = {};

my $out_file = shift @ARGV;
my $complete = shift @ARGV;

foreach my $hmm_tbl_file ( @ARGV ) {

    open ( my $in , "<" , $hmm_tbl_file ) or die "$!";

    while (<$in>) {

        chomp;
        
        next if $_ =~ m/^\#/;
        
        # seq_s_35267_38642_43749_LTR_LTR_Finder#LTR/Copia|aa2         -          AP_1731              -            1.5e-10   40.4   1.4   3.9e-10   39.1   1.4   1.6   1   1   0   1   1   1   1 -
        
        if ( $_ =~ m/^(\S+)\s+\-\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)/ ) {
        
            my ( $repeat , $domain , $e , $score ) =
            ( $1 , $2 , $3 , $4 , $5 );
            
            next if $e >= 0.01;
            
            my $domain_spec = "";
            
            if ( $domain =~ m/([^\_]+)\_(\S+)/ ) {
            
                $domain = $1;
                $domain_spec = "_" . $2;
            
            }
            
            $repeat =~ s/\|[^\|]+$//;
        
            if ( $domain =~ m/^GAG[^\_]/ ) {
            
                $domain = "GAG";
            
            }
            
            if (
            
                defined $repeat_table->{ $repeat }->{ $domain } and
                $e <
                $repeat_table->{ $repeat }->{ $domain }->{ 'e' }
            
            ) {
            
                $repeat_table->{ $repeat }->{ $domain }->{ 'e' } = $e;
                $repeat_table->{ $repeat }->{ $domain }->{ 'info' } = "${domain}$domain_spec,$e,$score";
                
            
            }
            
            else {
            
                $repeat_table->{ $repeat }->{ $domain }->{ 'e' } = $e;
                $repeat_table->{ $repeat }->{ $domain }->{ 'info' } = "${domain}$domain_spec,$e,$score";
            
            }
        
        }
        
        else {
        
            say $_;
        
        }

    }

}

my @repeats = keys %{ $repeat_table };

my @sorted_repeats = sort {
    ( scalar keys %{ $repeat_table->{ $b } } ) <=>
    ( scalar keys %{ $repeat_table->{ $a } } )
} @repeats;

my @counts;

open ( my $out , ">" , $out_file . ".hits.csv" ) or die "$!";
open ( my $c_out , ">" , $out_file . ".hits.complete.csv" ) or die "$!";
open ( my $i_out , ">" , $out_file . ".hits.incomplete.csv" ) or die "$!";

say $out "NR\tREPEAT\tN\tDOMAINS\tDOMAINS_SCORES";

my $i = 1;

foreach my $repeat ( @sorted_repeats ) {

    my @sorted_domains = sort keys %{ $repeat_table->{ $repeat } };

    my $n = @sorted_domains;
    
    my $sorted_domains_string = join ( "|" , @sorted_domains );
    
    my $scores;
    
    foreach my $domain ( @sorted_domains ) {
    
        $scores .= "|" . $repeat_table->{ $repeat }->{ $domain }->{ 'info' };
    
    }
    
    $scores =~ s/^.//;
    
    say $out "$i\t$repeat\t$n\t$sorted_domains_string\t$scores";

    $counts[$n]++;

    if ( $n >= $complete ) {
    
		say $c_out $repeat;
    
    }
    
    $i++;
    
}

open ( my $d_out , ">" , $out_file . ".hits.dist.csv" ) or die "$!";

say $d_out "NR\tN";

for ( my $i = $#counts; $i > 0; $i-- ) {

    $counts[ $i ] = 0 unless defined $counts[ $i ];

    say $d_out $i , "\t" , $counts[ $i ];

}

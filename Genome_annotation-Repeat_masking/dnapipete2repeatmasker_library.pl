#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $min_contig = shift @ARGV;
my $min_target = shift @ARGV;

open ( my $anno_in , "<" , $ARGV[ 0 ] ) or die "$!";

my $anno_table = {};

while (<$anno_in>) {

    chomp;
    
    # comp_TRINITY_DN10046_c1_g1_i1	205	0.975609756097561	DNA3-11_Mad	DNA/PIF-Harbinger	861	[320-518]	0.22996515679442509
    
    my @data = split ( /\t/ , $_ );
    
    my ( $contig , $tmp_min_contig , $repeat , $class , $tmp_min_target ) =
    @data[ 0 , 2 , 3 , 4 , 7 ];
    
    next if $class eq "Simple_repeat";
    next if $class eq "Low_complexity";
    next if $class eq "Satellite";
    
    if ( $tmp_min_contig >= $min_contig and $tmp_min_target >= $min_target ) {
    
        $anno_table->{ $contig } = "$class $repeat";
    
    }

}

open ( my $fasta_in , "<" , $ARGV[ 1 ] ) or die "$!";

open (
    my $fasta_out ,
    ">" ,
    $ARGV[ 1 ] . ".repeatmasker.contig_${min_contig}.template_${min_target}.lib"
) or die "$!";

my $tmp_header = undef;
my $tmp_class = undef;

my $contig;

my $stat_table = {};

while (<$fasta_in>) {

    chomp;

    if ( $_ =~ m/^\>(.+)/ ) {
    
        $tmp_header = undef;
        $tmp_class = undef;
    
        my $header = $1;
        $contig = $header;
        
        $contig =~ s/^.+\_comp/comp/;
    
        if ( defined $anno_table->{ $contig } ) {
        
            $tmp_header = "$header#" . $anno_table->{ $contig };
            
            $tmp_class = $anno_table->{ $contig };
            $tmp_class =~ s/\W.+//;
            
            say $fasta_out ">$tmp_header";

            $stat_table->{ $tmp_class }->{ 'n' }++;
            $stat_table->{ 'ANY' }->{ 'n' }++;
        
        }
    
    }
    
    elsif ( defined $tmp_header ) {
    
        say $fasta_out "$_";
        
        if ( defined $tmp_class ) {
        
            $stat_table->{ $tmp_class }->{ 'length' } += length $_;
            $stat_table->{ 'ANY' }->{ 'length' } += length $_;
        
        }
    
    }

}

open (
    my $stats_out ,
    ">" ,
    $ARGV[ 1 ] . ".repeatmasker.contig_${min_contig}.template_${min_target}.lib.stats.csv"
) or die "$!";

say $stats_out "CLASSIFICATION\tPROP\tLENGTH\tN\tAVG_LENGTH";

foreach my $tmp_class ( sort keys %{ $stat_table } ) {

    my $n = $stat_table->{ $tmp_class }->{ 'n' };
    my $length = $stat_table->{ $tmp_class }->{ 'length' };

    my $prop = $length / $stat_table->{ 'ANY' }->{ 'length' };
    
    next unless defined $length;
    
    say $stats_out "$tmp_class\t$prop\t$length\t$n\t" , $length / $n;
    
}

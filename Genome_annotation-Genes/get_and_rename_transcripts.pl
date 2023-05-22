#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $tag = shift @ARGV;
my $keep_file = shift @ARGV;

open ( my $in , "<" , $keep_file ) or die "$!";

my $locus_table = {};

my $header = <$in>;

while (<$in>) {

    chomp;

    my @data = split ( /\t/ , $_ );
    
    $locus_table->{ $data[ 1 ] }->{ 'transcript' } = $data[ 0 ];
    
    $locus_table->{ $data[ 1 ] }->{ 'string' } =
    "score:" . $data[ 3 ] . "|" . "target:" . $data[ 2 ] . "|" .
    "n_species:" . $data[ 4 ] . "|" . "species:" . $data[ 5 ];
    
    
}

my $gxf_file = shift @ARGV;

open ( my $out , ">" , $gxf_file . ".kept" ) or die "$!";

open ( my $gxf_in , "<" , $gxf_file ) or die "$!";

my $print = 0;

my $transcript;
my $locus;

my $gene_id;
my $transcript_id;

while (<$gxf_in>) {

    chomp;
    
    my $line = $_;
    
    if ( $_ =~ m/^\#/ ) {
    
        say $out $line;
        next;
    
    }
    
    elsif ( $_ =~ m/\stranscript\s.+transcript_id\s\"([^\"]+)\.p\d+\"\;\sgene_id\s\"([^\"]+)\"/ ) {
    
        $print = 0;
    
        ( $transcript , $locus ) = ( $1 , $2 );
        
        my $locus_num = $locus;
        $locus_num =~ s/.+\_//;
        
        $gene_id = "${tag}.$locus_num";
        $transcript_id = "${tag}.${locus_num}.1";
        
        $line =~ s/\"$transcript[^\"]+\"/\"$transcript_id\"/;
        $line =~ s/\"$locus\"/\"$gene_id\"/;
    
        if ( defined $locus_table->{ $locus } and $transcript eq $locus_table->{ $locus }->{ 'transcript' } ) {
        
            $print = 1;
            
            my $info = $locus_table->{ $locus }->{ 'string' };
            
            $line .= "; " . "info \"$info\";";
        
        }
    
    }
    
    # seq_c_100039    transdecoder    mRNA    12331   43892   .       +       .       ID=COMPARATIVE_00000001.p1;Parent=XLOC_000001;
    
    elsif ( $line =~ m/\sgene\s.+\sID\=([^\;]+)/ ) {
    
        ( $locus ) = $1;
        
        my $locus_num = $locus;
        $locus_num =~ s/.+\_//;
        
        my $gene_id = "${tag}.$locus_num";
        
        if ( defined $locus_table->{ $locus } ) {

            $line =~ s/\"$locus\"/\"$gene_id\"/;
        
            say $out $line;
            next;
        
        }
    
    }
    
    elsif ( $line =~ m/\smRNA\s.+\sID=([^\;]+)\.p\d+\;Parent\=([^\;]+);/ ) {
           
        $print = 0;

        ( $transcript , $locus ) = ( $1 , $2 );
        
        my $locus_num = $locus;
        $locus_num =~ s/.+\_//;
        
        $gene_id = "${tag}.$locus_num";
        $transcript_id = "${tag}.${locus_num}.1";
        
        $line =~ s/\=$transcript[^\"]+\;/\=$transcript_id\;/;
        $line =~ s/\=$locus\;/\=$gene_id\;/;
        
        if ( defined $locus_table->{ $locus } and $transcript eq $locus_table->{ $locus }->{ 'transcript' } ) {
        
            $print = 1;
            
            my $info = $locus_table->{ $locus }->{ 'string' };
            
            $line .= ";Info=$info;";
        
        }
    
    }
    
    if ( $print ) {
    
        $line =~ s/\"$transcript[^\"]+\"/\"$transcript_id\"/;
        $line =~ s/\"$locus\"/\"$gene_id\"/;
        
        $line =~ s/\=$transcript[^\"]+\;/\=$transcript_id\;/;
        $line =~ s/\=$locus\;/\=$gene_id\;/;
    
        say $out $line;
    
    }

}

#!/usr/bin/perl

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

use Getopt::Long;
use Data::Dumper;

# Main table for command line options, defaults and data structures

my $opt = {

	# Minimum quality to reprint read mappings
	
	'minimum_quality' => 60 ,

	# Edge distance to reprint read mappings
	
	'maximum_distance' => 20000 ,
	
	'read_length' => 150 ,
	
	# ( 134 + 150 ) / 2 = 142
	
};

GetOptions(

	$opt ,

	# Input BAM

	'bam|b=s' ,
	
	'minimum_quality|min_quality|min_qual|mq=i' ,
	
	'maximum_distance|max_distance|max_dist|md=i' ,
	
	'read_length|rl=s' ,

);

my $seq_length_table = {};

my $all_seq_table = {};
my $keep_seq_table = {};

my $all_barcode_table = {};
my $keep_barcode_table = {};
my $dist_table = [];

my $minimum_quality = $opt->{ 'minimum_quality' };
my $maximum_distance = $opt->{ 'maximum_distance' };

my $all_mappings_n = 0;
my $keep_mappings_n = 0;

my $tot_length;
my $tot_keep_length;

scan_data( $opt );

sub scan_data {

    my ( $opt ) = @_;
    
    my $bam_file = $opt->{ 'bam' };
    
    unless ( -f $bam_file ) {
        
        die "Can not read from BAM file " , $bam_file;
    
    }
    
    open ( my $bam_in , "samtools view -H $bam_file |" ) or die "$!";
    
    my $bam_out_file = "${bam_file}.analysed.edge_${maximum_distance}.qual_${minimum_quality}.bam";
    
    open ( my $bam_out , "| samtools view -O BAM -o $bam_out_file -" ) or die "$!";
    
    while (<$bam_in>) {
    
        chomp;
        
        say $bam_out $_;
        
        if ( $_ =~ m/^\@SQ\t+SN:([^\t]+)\tLN:(\S+)/ ) {
        
            $seq_length_table->{ $1 }->{ 'length' } = $2;
            
        }
    
    }
    
    close $bam_in;
    
    open ( $bam_in , "samtools view $bam_file |" ) or die "$!";
    
    my $last_seq = "";
    
    my $seq_length = "";
    
    while (<$bam_in>) {
    
        chomp;
        
        if ( $_ =~ m/^\S+\t\S+\t(\S+)\t+(\d+)\t(\d+)\t.+\tBX:Z:([^\t]+)\t/ ) {
        
            if ( $1 ne $last_seq ) {
            
                if ( $last_seq ) {
            
                    say $last_seq , "\t" ,
                        $seq_length , "\t" ,
                        $all_mappings_n , "\t" , $keep_mappings_n;
                
                }
                
                $seq_length = $seq_length_table->{ $1 }->{ 'length' };
            
            }
            
            my $dist = $2;
            
            my $stop_dist = $seq_length - $2;
            
            if ( $stop_dist < $dist ) {
            
                $dist = $stop_dist;
            
            }
            
            if ( $dist <= $maximum_distance and $3 >= $minimum_quality ) {

                say $bam_out $_;
            
                $keep_mappings_n++;
            
            }
            
            $all_mappings_n++;
            
            $last_seq = $1;
        
        }
    
    }
    
    close $bam_in;
    close $bam_out;
    
    system ( "samtools index $bam_out_file" );
    system ( "samtools flagstat $bam_out_file > ${bam_out_file}.flagstat.csv" );
    system ( "samtools idxstats $bam_out_file > ${bam_out_file}.idxstats.csv" );
    
}

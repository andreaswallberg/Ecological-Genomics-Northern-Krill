#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $chosen_file = shift @ARGV;

my $chosen_table = {};

open ( my $chosen_in , "<" , $chosen_file ) or die "$!";

while (<$chosen_in>) {

    chomp;
    
    my ( $gene , $transcript ) = split ( /\t/ , $_ );

    $chosen_table->{ $transcript } = $gene;
    
}

my $gxf_file = shift @ARGV;

open ( my $gxf_in , "<" , $gxf_file ) or die "$!";

my $print = 0;

my $transcript;
my $gene;
my $old_id;

my $gene_id;
my $transcript_id;

my $tag_table = {};

while (<$gxf_in>) {

    chomp;
    
    my $line = $_;
       
    if ( $_ =~ m/^\#/ ) {
    
        say $line;
        next;
    
    }
    
    # seq_c_100039	transdecoder	transcript	42721	43512	.	+	.	transcript_id "REFERENCE_00000001"; gene_id "XLOC_000001"; oId "STRG.ref_mixed.103811.1.p1"; tss_id "TSS1"; num_samples "1";
    
    # transcript_id "REFERENCE_BEST_00009514"; gene_id "XLOC_005711"; oId "STRG.ref_mixed.44637.1.p1"; tss_id "TSS8312"; num_samples "1"; info "score:793|target:ref_mixed|n_species:1|species:ref_mixed";
    
    elsif ( $line =~ m/\stranscript\s.+transcript_id\s\"([^\"]+)\"\;\sgene_id\s\"([^\"]+)\"\;\soId\s\"([^\"]+)\"\;\s/ ) {
    
        $print = 1;
    
        ( $transcript , $gene , $old_id ) = ( $1 , $2 , $3 );

        if ( $transcript =~ m/^(\S+)\.(\d+)/ ) {
        
            my ( $trans , $num ) = ( $1 , $2 );
        
            if ( defined $chosen_table->{ $transcript } ) {

                my $tag = "GENE";
            
                if ( $old_id =~ m/STRG\.(\d+)/ ) {
            
                    $tag = "REF_STRG_1_" . $1;
            
                }
                
                elsif ( $old_id =~ m/COMPARATIVE_SPALN\.(\d+\.\d+)\.\d+/ ) {
            
                    $tag = "COM_SPAL_" . $1;
            
                }
                
                elsif ( $old_id =~ m/mRNA\.ref\.(\d+\.\d+)/ ) {
            
                    $tag = "REF_TRIN_" . $1;
            
                }
                
                $tag =~ s/\W/\_/g;
        
                my $gene_tag = $tag . "_" . $gene;
                
                $tag_table->{ $gene }->{ 'tag' } = $gene_tag;
                $tag_table->{ $gene }->{ 'chosen' } = $transcript;
                
                say STDERR "Building $gene_tag";
        
            }
            
        }
    
    }

}

close $gxf_in;

open ( $gxf_in , "<" , $gxf_file ) or die "$!";

my $gene_tag;
my $trans_tag;

while (<$gxf_in>) {

    chomp;
    
    my $line = $_;
       
    if ( $_ =~ m/^\#/ ) {
    
        say $line;
        next;
    
    }
    
    # seq_c_100039	transdecoder	transcript	42721	43512	.	+	.	transcript_id "REFERENCE_00000001"; gene_id "XLOC_000001"; oId "STRG.ref_mixed.103811.1.p1"; tss_id "TSS1"; num_samples "1";
    
    # transcript_id "REFERENCE_BEST_00009514"; gene_id "XLOC_005711"; oId "STRG.ref_mixed.44637.1.p1"; tss_id "TSS8312"; num_samples "1"; info "score:793|target:ref_mixed|n_species:1|species:ref_mixed";
    
    elsif ( $line =~ m/\stranscript\s.+transcript_id\s\"([^\"]+)\"\;\sgene_id\s\"([^\"]+)\"\;\soId\s\"([^\"]+)\"\;\s/ ) {
    
        $print = 0;
    
        ( $transcript , $gene , $old_id ) = ( $1 , $2 , $3 );

        if ( $transcript =~ m/^(\S+)\.(\d+)/ ) {
        
            my ( $trans , $num ) = ( $1 , $2 );
            
            if ( defined $tag_table->{ $gene }->{ 'tag' } ) {
            
                $print = 1;
                
                $gene_tag = $tag_table->{ $gene }->{ 'tag' };
                $trans_tag = $tag_table->{ $gene }->{ 'tag' } . "." . $num;
            
                if ( defined $chosen_table->{ $transcript } ) {
                
                    $line .= " chosen_isoform \"1\";";
                    
                }
            
            }
            
            else {
            
                die "No such gene:\n$_";
            
            }
            
        }
    
    }
    
    if ( $print ) {
    
        $line =~ s/transcript_id\s\"[^\"]+\"\;/transcript_id "${trans_tag}"\;/;
        $line =~ s/gene_id\s\"[^\"]+\"\;/gene_id "${gene_tag}"\;/;
        
        say $line;
    
    }

}

#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;
use Data::Dumper;

my $stop_plus_table = {

    "TAA" => "TAA" ,
    "TGA" => "TGA" ,
    "TAG" => "TAG" ,

};

my $stop_minus_table = {

    "TTA" => "TAA" ,
    "TCA" => "TGA" ,
    "CTA" => "TAG" ,

};

my $seq_table = {};

my $fasta_file = $ARGV[ 0 ];

open ( my $fasta_in , "<" , $fasta_file ) or die "$!";

say STDERR "Reading FASTA sequence from $fasta_file ...";

my $tmp_header;

while (<$fasta_in>) {

    chomp;
    
    if ( $_ =~ m/^\>(.+)/ ) {
    
        $tmp_header = $1;
    
    }
    
    else {
    
        $seq_table->{ $tmp_header } = uc $_;
    
    }

}

my $gff_file = $ARGV[ 1 ];

open ( my $gff_in , "<" , $gff_file ) or die "$!";

say STDERR "Reading GFF coordinates from $gff_file ...";

my $last_table = {};

my $in_3prime_partial = 0;

my $three_prime_partial_table = {};

while (<$gff_in>) {

    chomp;
    
    if ( $_ =~ m/\sgene\s/ ) {
    
        $in_3prime_partial = 0;
        
        say STDERR "Resetting 3-prime partial or internal status ...";
    
    }
    
    if ( $_ =~ m/^\S+.+\s+gene\s+\d+\s+\d+\s+\S+\s+\S\s+.+type\%3A(3prime_partial|internal)\%20/ ) {
    
        $in_3prime_partial = 1;
    
        say STDERR "Gene is 3-prime partial or internal ...";
    
    }
    
    elsif ( $in_3prime_partial and $_ =~ m/^(\S+).+\s+CDS\s+(\d+)\s+(\d+)\s+\S+\s+(\S)\s+(\S+)\s+.+Parent=(\S+)\.\d+\.p\d+/ ) {
    
        my ( $header , $start , $stop , $orient , $phase , $gene ) = ( $1 , $2 , $3 , $4 , $5 , $6 );
        
        my $length = $stop - $start + 1;
        
        my $modulo_length = $length - $phase;

        if ( $modulo_length % 3 ) {
        
            say STDERR "Skipping CDS with out of phase CDS!";
        
            next;
        
        }
        
        if ( defined $seq_table->{ $header } and $orient eq "+" ) {
        
            my $gene_seq = substr ( $seq_table->{ $header } , $start - 1 , $length );
            
            my $gene_length = length $gene_seq;
            
            my $last_codon = substr ( $gene_seq , $gene_length - 3 , 3 );
        
            if ( not defined $stop_plus_table->{ $last_codon } ) {
        
                $three_prime_partial_table->{ $gene }->{ 'header' } = $header;
                $three_prime_partial_table->{ $gene }->{ 'orient' } = $orient;
                $three_prime_partial_table->{ $gene }->{ 'start' } = $start;
                $three_prime_partial_table->{ $gene }->{ 'stop' } = $stop;
                $three_prime_partial_table->{ $gene }->{ 'last_codon' } = $last_codon;
            
            }
        
        }
        
        elsif ( defined $seq_table->{ $header } and $orient eq "-" ) {
        
            my $gene_seq = substr ( $seq_table->{ $header } , $start - 1 , $length );
            
            my $gene_length = length $gene_seq;
            
            my $last_codon = substr ( $gene_seq , 0 , 3 );
        
            if ( not defined $stop_minus_table->{ $last_codon } ) {
        
                    $three_prime_partial_table->{ $gene }->{ 'header' } = $header;
                    $three_prime_partial_table->{ $gene }->{ 'orient' } = $orient;
                    $three_prime_partial_table->{ $gene }->{ 'start' } = $start;
                    $three_prime_partial_table->{ $gene }->{ 'stop' } = $stop;
                    $three_prime_partial_table->{ $gene }->{ 'last_codon' } = $last_codon;
            
            }
        
        }
    
    }

}

my $i = 0;

foreach my $gene ( keys %{ $three_prime_partial_table } ) {

    my $header = $three_prime_partial_table->{ $gene }->{ 'header' };
    
    if ( $three_prime_partial_table->{ $gene }->{ 'orient' } eq "+" ) {
    
        my $stop = $three_prime_partial_table->{ $gene }->{ 'stop' };

        my $next_codon = substr ( $seq_table->{ $header } , $stop , 3 );

        if ( length $next_codon == 3 ) {

            if ( defined $stop_plus_table->{ $next_codon } ) {

                $three_prime_partial_table->{ $gene }->{ 'stop_codon' } = $stop_plus_table->{ $next_codon };
                $three_prime_partial_table->{ $gene }->{ 'stop_codon_outer' } = $stop + 3;
                
                $last_table->{ $next_codon }++;
        
                $i++;
        
                say "Added stop codon for plus strand gene ($i)";

            }
        
        }
    
    }
    
    if ( $three_prime_partial_table->{ $gene }->{ 'orient' } eq "-" ) {
    
        my $start = $three_prime_partial_table->{ $gene }->{ 'start' };

        next unless $start > 4;
        
        my $next_codon = substr ( $seq_table->{ $header } , $start - 4 , 3 );

        if ( length $next_codon == 3 ) {

            if ( defined $stop_minus_table->{ $next_codon } ) {

                $three_prime_partial_table->{ $gene }->{ 'stop_codon' } = $stop_minus_table->{ $next_codon };
                $three_prime_partial_table->{ $gene }->{ 'stop_codon_outer' } = $start - 3;
        
                $last_table->{ $next_codon }++;
        
                $i++;
        
                say "Added stop codon for minus strand gene ($i)";

            }
        
        }
    
    }

}

close $gff_in;

open ( $gff_in , "<" , $gff_file ) or die "$!";
open ( my $gff_out , ">" , $gff_file . ".added_stop_codons.gff" ) or die "$!";

say STDERR "Updating GFF coordinates for $gff_file ...";

while (<$gff_in>) {

    chomp;
    
    my @data = split ( /\s/ , $_ );
    
    if ( $_ =~ m/\s+(gene|mRNA)\s+/ ) {

        my ( $orient , $gene ) = @data[ 6 , 8 ];
        
        $gene =~ s/^ID\=//;
        $gene =~ s/\;.*//;
        $gene =~ s/\..*//;
        
        if ( $three_prime_partial_table->{ $gene }->{ 'stop_codon' } ) {
        
            if ( $orient eq "+" ) {
            
                $data[ 4 ] = $three_prime_partial_table->{ $gene }->{ 'stop_codon_outer' };
            
            }
            
            if ( $orient eq "-" ) {
            
                $data[ 3 ] = $three_prime_partial_table->{ $gene }->{ 'stop_codon_outer' };
            
            }
            
            if ( $data[ 8 ] =~ m/3A3prime_partial/ ) {
            
                $data[ 8 ] =~ s/3A3prime_partial/3Acomplete/g;
        
            }
            
            elsif ( $data[ 8 ] =~ m/3Ainternal/ ) {
            
                $data[ 8 ] =~ s/3Ainternal/3A5prime_partial/g;
            
            }
        
        }
    
    }
    
    if ( $_ =~ m/\s+(CDS|exon)\s+/ ) {

        my ( $orient , $gene ) = @data[ 6 , 8 ];
        
        $gene =~ s/^ID\=//;
        $gene =~ s/cds\.//;
        $gene =~ s/\;.*//;
        $gene =~ s/\..*//;
        
        if ( $three_prime_partial_table->{ $gene }->{ 'stop_codon' } ) {
        
            if (
            
                $data[ 3 ] eq $three_prime_partial_table->{ $gene }->{ 'start' } and
                $data[ 4 ] eq $three_prime_partial_table->{ $gene }->{ 'stop' }
            
            ) {
        
                if ( $orient eq "+" ) {
                
                    $data[ 4 ] = $three_prime_partial_table->{ $gene }->{ 'stop_codon_outer' };
                
                }
                
                if ( $orient eq "-" ) {
                
                    $data[ 3 ] = $three_prime_partial_table->{ $gene }->{ 'stop_codon_outer' };
                
                }
            
            }
        
        }
    
    }
    
    say $gff_out join ( "\t" , @data );

}

$fasta_file = $ARGV[ 2 ];

open ( $fasta_in , "<" , $fasta_file ) or die "$!";

open ( my $fasta_out , ">" , $fasta_file . ".added_stop_codons.fasta" ) or die "$!";

my $gene_table = {};

my @genes;

my $gene = undef;

while (<$fasta_in>) {

    chomp;
    
    if ( $_ =~ m/^\>([^\.]+)/ ) {
    
        $gene = $1;
        
        unless ( defined $gene_table->{ $gene } ) {
        
            $gene_table->{ $gene }->{ 'name' } = $_;
        
            push @genes , $gene;
    
        }
    
    }
    
    else {
    
        $gene_table->{ $gene }->{ 'seq' } .= $_;
    
    }

}

foreach my $gene ( @genes ) {

    my $name = $gene_table->{ $gene }->{ 'name' };
    my $seq = $gene_table->{ $gene }->{ 'seq' };
    
    if ( defined $three_prime_partial_table->{ $gene }->{ 'stop_codon' } ) {
    
        $seq .= $three_prime_partial_table->{ $gene }->{ 'stop_codon' };
        
    }
    
    say $fasta_out "${name}\n" , $seq;

}

# say Dumper ( $last_table );

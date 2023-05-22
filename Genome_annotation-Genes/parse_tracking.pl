#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $tracking_file = shift @ARGV;

open ( my $in , "<" , $tracking_file ) or die "$!";

my $locus_table = {};
my @loci;

my $species_table = {};

while (<$in>) {

    chomp;
    
    # COMPARATIVE_00000147	XLOC_000049	-	u	-	q2:gene.tepidariorum.6.00042|mRNA.tepidariorum.6.00042|2|0.000000|0.000000|0.000000|195	-	-	q5:gene.longiramus.6.00042|mRNA.longiramus.6.00042|2|0.000000|0.000000|0.000000|195	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
    
    my @data = split ( /\t/ , $_ );
    
    my ( $model , $locus ) = @data[ 0 , 1 ];
    
    for ( my $i = 4; $i < @data; $i++ ) {
    
        my $transcript_string = $data[ $i ];
        
        if ( $transcript_string =~ /^[^\:]+\:\w+\.(\w+)\./ ) {
        
            my $species = $1;
        
            my @transcript_data = split ( /\|/ , $transcript_string );
        
            my $n_exons = $transcript_data[ 2 ];
            my $exons_length = $transcript_data[ -1 ];
            
            if ( $n_exons >= 2 and $exons_length >= 500 ) {
        
                say "$model\t$locus\t$species";
            
                unless ( defined $locus_table->{ $locus } ) {
                
                    push @loci , $locus;
                
                }
            
                push @{ $locus_table->{ $locus }->{ 'species' }->{ $species } } , $model;
                push @{ $species_table->{ $species }->{ 'locus' }->{ $locus } } , $model;
        
            }
        
        }
    
    }

}

my $model_table = {};

foreach my $blast_file ( @ARGV ) {

    open ( my $blast_in , "<" , $blast_file ) or die "$!";

    while (<$blast_in>) {

        chomp;
        
        next if $_ =~ m/^\#/;
        
        if ( $_ =~ m/^(\S+)\.p\d+\s.+\s(\S+)$/ ) {
        
            unless ( defined $model_table->{ $1 } ) {
            
                $model_table->{ $1 } = $2;
            
            }
            
            if ( defined $model_table->{ $1 } and $2 > $model_table->{ $1 } ) {
            
                $model_table->{ $1 } = $2;
            
            }
        
        }

    }

}

my $locus_dist = [];
my $score_dist = [];
my $species_dist = {};

my $keep_locus_table = {};

my $n_loci = @loci;

open ( my $keep_out , ">"  , "${tracking_file}.keep_transcripts.csv" ) or die "$!";

say $keep_out "TRANSCRIPT\tLOCUS\tSPECIES\tSCORE\tN_SPECIES\tSPECIES";

foreach my $locus ( @loci ) {

    my $keep = 0;

    my @all_species = keys %{ $locus_table->{ $locus }->{ 'species' } };
    
    my $n_species = @all_species;
    
    $locus_dist->[ $n_species ]->{ 'ALL' }++;
        
    if ( $n_species >= 2 and defined $locus_table->{ $locus }->{ 'species' }->{ 'superba' } ) {
    
        $keep = 1;
    
    }
    
    elsif ( $n_species >= 2 ) {
    
        if (
        
            defined $locus_table->{ $locus }->{ 'species' }->{ 'azteca' } or
            defined $locus_table->{ $locus }->{ 'species' }->{ 'hawaiensis' } or
            defined $locus_table->{ $locus }->{ 'species' }->{ 'longiramus' } or
            defined $locus_table->{ $locus }->{ 'species' }->{ 'californicus' } or
            defined $locus_table->{ $locus }->{ 'species' }->{ 'californicus' } or
            defined $locus_table->{ $locus }->{ 'species' }->{ 'sinensis' } or
            defined $locus_table->{ $locus }->{ 'species' }->{ 'quadricarinatus' } or
            defined $locus_table->{ $locus }->{ 'species' }->{ 'virginalis' } or
            defined $locus_table->{ $locus }->{ 'species' }->{ 'americanus' } or
            defined $locus_table->{ $locus }->{ 'species' }->{ 'monodon' } or
            defined $locus_table->{ $locus }->{ 'species' }->{ 'vannamei' } or
            defined $locus_table->{ $locus }->{ 'species' }->{ 'magna' } or
            defined $locus_table->{ $locus }->{ 'species' }->{ 'pulex' } or
            defined $locus_table->{ $locus }->{ 'species' }->{ 'vulgare' } or
            defined $locus_table->{ $locus }->{ 'species' }->{ 'monodon' }
    
        ) {
        
            $keep = 1;
        
        }
    
    }
    
    my $max_score = 0;
    my $max_species = "";
    my $max_model = "";
    
    foreach my $species ( @all_species ) {
    
        $species_dist->{ $species }->{ 'PRESENT' }++;
    
        foreach my $model ( @{ $locus_table->{ $locus }->{ 'species' }->{ $species } } ) {

            if ( defined $model_table->{ $model } ) {
        
                my $score = $model_table->{ $model };
                
                if ( $score > 100 and $score > $max_score ) {
                
                    $max_score = $score;
                    $max_species = $species;
                    $max_model = $model;
                
                }
            
            }
        
        }
        
    }

    if ( $max_score ) {
    
        my $score_i = int ( $max_score / 10 );
        
        $score_dist->[ $score_i ]->{ 'ALL' }++;
        $species_dist->{ $max_species }->{ 'ALL' }++;
    
        if ( $keep ) {
    
            $keep_locus_table->{ $locus } = 1;
            $locus_dist->[ $n_species ]->{ 'KEEP' }++;
        
            foreach my $species ( keys %{ $locus_table->{ $locus }->{ 'species' } } ) {
                
                $species_dist->{ $species }->{ 'KEEP_PRESENT' }++;
                
            }
            
            $score_dist->[ $score_i ]->{ 'KEEP' }++;
            $species_dist->{ $max_species }->{ 'KEEP' }++;
            
            say $keep_out "${max_model}\t${locus}\t${max_species}\t${max_score}\t" ,
                "$n_species\t" , join ( "," , @all_species );
        
        }

    }

}

my $n_keep_loci = keys %{ $keep_locus_table };

my @headers = ( 'ALL' , 'KEEP' );

open ( my $dist_out , ">"  , "${tracking_file}.dist.csv" ) or die "$!";

say $dist_out "N_SPECIES\t" , join ( "\t" , map { $_ . "_N_GENES"  } @headers );

for ( my $i = 1; $i < @{ $locus_dist }; $i++ ) {

    print $dist_out "$i";
    
    foreach my $type ( @headers ) {
    
        $locus_dist->[ $i ]->{ $type } = 0 unless defined
        $locus_dist->[ $i ]->{ $type };
        
        print $dist_out "\t" , $locus_dist->[ $i ]->{ $type };
    
    }

    say $dist_out "";
    
}

open ( my $b_dist_out , ">"  , "${tracking_file}.dist_blast.csv" ) or die "$!";

say $b_dist_out "SCORE\t" , join ( "\t" , map { $_ . "_N_GENES"  } @headers );

for ( my $i = 0; $i < @{ $score_dist }; $i++ ) {

    print $b_dist_out $i * 10;
    
    foreach my $type ( @headers ) {
    
        $score_dist->[ $i ]->{ $type } = 0 unless defined
        $score_dist->[ $i ]->{ $type };
        
        print $b_dist_out "\t" , $score_dist->[ $i ]->{ $type };
    
    }

    say $b_dist_out "";
    
}

open ( my $s_dist_out , ">"  , "${tracking_file}.dist_max_species.csv" ) or die "$!";

say $s_dist_out "SPECIES\tN_LOCI\tPRESENT\tPROP_PRESENT\tN_KEEP_LOCI\tKEEP_PRESENT\tPROP_KEEP_PRESENT\t" , join ( "\t" , ( map { $_ . "_MAX" } @headers ) );

foreach my $species ( sort keys %{ $species_dist } ) {

    $species_dist->{ $species }->{ "PRESENT" } = 0 unless defined
    $species_dist->{ $species }->{ "PRESENT" };

    print $s_dist_out $species , "\t" , $n_loci , "\t" ,
    $species_dist->{ $species }->{ "PRESENT" } , "\t" ,
    $species_dist->{ $species }->{ "PRESENT" } / $n_loci;
    
    $species_dist->{ $species }->{ "KEEP_PRESENT" } = 0 unless defined
    $species_dist->{ $species }->{ "KEEP_PRESENT" };
    
    print $s_dist_out "\t" , $n_keep_loci , "\t" ,
    $species_dist->{ $species }->{ "KEEP_PRESENT" } , "\t" ,
    $species_dist->{ $species }->{ "KEEP_PRESENT" } / $n_keep_loci;
    
    foreach my $type ( ( @headers ) ) {
    
        $species_dist->{ $species }->{ $type } = 0 unless defined
        $species_dist->{ $species }->{ $type };
    
        print $s_dist_out "\t" , $species_dist->{ $species }->{ $type };
    
    }

    say $s_dist_out "";
    
}

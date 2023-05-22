#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

use Getopt::Long;
use Data::Dumper;

# Option table and defaults

my $opt = {

	'min_hits' => 10 ,
	'min_identity' => 80 ,
	'min_prop' => 0.8 ,
	'min_length' => 1500 ,
	
	'repeats' => {} ,
	'reference' => [] ,
	'blast' => {} ,
	
	'repeat_table' => {} ,
	'tot_repeat_table' => {} ,
	
	'reference_table' => {} ,
	'tot_reference_table' => {} ,
	
};

GetOptions(
    $opt ,

    # Input
    # -----

    'min_hits|mh=i' ,
    'min_identity|mi=s' ,
    'min_prop|mp=s' ,
    'min_length|ml=i' ,
    
    'repeats|rep|r=s{,}' ,
    'reference|ref|f=s{,}' ,
    
    'repeat_tag|tag|t=s' ,
    
    'blast|b=s{,}' ,
    
    # Output
    # ------
    
    # Output name (string)
    'output|out|o=s' ,
    
    'verbose|v' ,

);

read_seq_lengths( $opt );
read_blast_hits( $opt );
summarize_hits_dists( $opt );
summarize_hits_reference( $opt );
print_repeat_lists( $opt );

sub read_seq_lengths {

    my ( $opt ) = @_;
    
    my $repeat_table = $opt->{ 'repeat_table' };
    my $tot_repeat_table = $opt->{ 'tot_repeat_table' };
    
    my $reference_table = $opt->{ 'reference_table' };
    my $tot_reference_table = $opt->{ 'tot_reference_table' };
    
    my $tag = $opt->{ 'repeat_tag' };
    
    foreach my $repeat_group ( sort keys %{ $opt->{ 'repeats' } } ) {
    
        my $repeat_file = $opt->{ 'repeats' }->{ $repeat_group };
    
        if ( -e $repeat_file ) {
        
            open ( my $in , "<" , $repeat_file ) or die "$!";
    
            say "Reading repeats for $repeat_group from $repeat_file ...";
    
            my $tmp_header;
            
            while (<$in>) {
            
                chomp;
                
                if ( $_ =~ m/^\>(\S+)/ ) {
                
                    my $test_header = $1;
                
                    $tmp_header = undef;
                    
                    if ( defined $tag ) {
                    
                        if ( $test_header =~ m/${tag}/ ) {
                        
                            $tmp_header = $test_header;
                            $tot_repeat_table->{ $repeat_group }->{ 'n' }++;
                
                        }
                
                    }
                    
                    else {
                    
                        $tmp_header = $test_header;
                        $tot_repeat_table->{ $repeat_group }->{ 'n' }++;
                    
                    }
                
                }
                
                elsif ( defined $tmp_header ) {
                
                    my $seq_length = length $_;
                
                    $repeat_table->{ $repeat_group }->
                    { $tmp_header }->{ 'length' } += $seq_length;
                    
                    $tot_repeat_table->{ $repeat_group }->{ 'length' } +=
                    $seq_length;
                
                }
            
            }
            
            say "Read " ,
            $tot_repeat_table->{ $repeat_group }->{ 'n' } ,
            " repeats totalling " ,
            $tot_repeat_table->{ $repeat_group }->{ 'length' } ,
            " bp";
        
        }
    
    }
    
    foreach my $reference ( @{ $opt->{ 'reference' } } ) {
        
        open ( my $in , "<" , $reference ) or die "$!";
        
        say "Reading genome sequences from $reference ...";
        
        my $tmp_header;
        
        while (<$in>) {
        
            chomp;
            
            if ( $_ =~ m/^\>(\S+)/ ) {
            
                $tmp_header = $1;
                
                $tot_reference_table->{ 'n' }++;
            
            }
            
            else {
            
                my $seq_length = length $_;
            
                $reference_table->{ $tmp_header }->{ 'length' } += $seq_length;
                
                $tot_reference_table->{ 'length' } += $seq_length;
            
            }
        
        }
        
        say "Read " ,
        $tot_reference_table->{ 'n' } ,
        " seqs totalling " ,
        $tot_reference_table->{ 'length' } ,
        " bp";
    
    }

}

sub read_blast_hits {

    my ( $opt ) = @_;
    
    my $repeat_table = $opt->{ 'repeat_table' };
    my $tot_repeat_table = $opt->{ 'tot_repeat_table' };
    
    my $reference_table = $opt->{ 'reference_table' };
    my $tot_reference_table = $opt->{ 'tot_reference_table' };
    
    my $min_hits = $opt->{ 'min_hits' };
    my $min_identity = $opt->{ 'min_identity' };
    my $min_length = $opt->{ 'min_length' };
    my $min_prop = $opt->{ 'min_prop' };
    
    foreach my $repeat_group ( sort keys %{ $opt->{ 'blast' } } ) {
    
        my $blast_file = $opt->{ 'blast' }->{ $repeat_group };
    
        my $blast_in;
    
        if ( $blast_file =~ m/\.gz$/ ) {
        
            open ( $blast_in , "zcat $blast_file |" ) or die "$!";
        
        }
        
        else {
        
            open ( $blast_in , "<" , $blast_file ) or die "$!";
        
        }
        
        say "Reading blast hits from $blast_file ...";
        
        my $last_repeat = "";
        
        my @all_hits;
        my $n_hits = 0;
        
        my $repeat_length;
        
        while (<$blast_in>) {
        
            next if $_ =~ m/^\#/;
        
            chomp;
            
            my @data = split ( /\t/ , $_ );
        
            my (
                $repeat ,
                $ref_seq ,
                $identity ,
                $al_length ,
                $ref_start ,
                $ref_stop
            ) = @data[ 0 , 1 , 2 , 3 , 8 , 9 ];
            
            unless ( $repeat eq $last_repeat ) {
            
                if (
                
                    $n_hits >= $min_hits and
                    defined $repeat_table->{ $repeat_group }->{ $last_repeat }
                
                ) {
                
                    register_hits(
                        $opt ,
                        $last_repeat ,
                        $repeat_group ,
                        \@all_hits ,
                    );
                
                }
                
                if (
            
                    defined $repeat_table->{ $repeat_group }->{ $repeat }
                    
                ) {
                
                    $repeat_length =
                    $repeat_table->{ $repeat_group }->{ $repeat }->{ 'length' };
                
                }
                
                @all_hits = ();
                $n_hits = 0;
            
            }
            
            if (
            
                defined $repeat_table->{ $repeat_group }->{ $repeat } and
                $identity >= $min_identity
            
            ) {
            
                my $repeat_prop = $al_length / $repeat_length;
            
                if ( $al_length >= $min_length or $repeat_prop >= $min_prop ) {
            
                    $n_hits++;
                
                    push @all_hits , [
                        $ref_seq ,
                        $identity ,
                        $ref_start ,
                        $ref_stop ,
                    ];
                
                }
            
            }
            
            $last_repeat = $repeat;
            
            if ( defined $tot_reference_table->{ 'tot_n' } ) {
            
                # last if $tot_reference_table->{ 'tot_n' } >= 100;
        
            }
        
        }
        
        if (
        
            $n_hits >= $min_hits and
            defined $repeat_table->{ $repeat_group }->{ $last_repeat }
            
        ) {
        
            register_hits(
                $opt ,
                $last_repeat ,
                $repeat_group ,
                \@all_hits ,
            );
        
        }
    
    }

}

sub register_hits {

    my ( $opt , $repeat , $repeat_group , $all_hits_ref ) = @_;

    my $repeat_table = $opt->{ 'repeat_table' };
    my $tot_repeat_table = $opt->{ 'tot_repeat_table' };
    
    my $reference_table = $opt->{ 'reference_table' };
    my $tot_reference_table = $opt->{ 'tot_reference_table' };
    
    $tot_reference_table->{ 'tot_n' }++;
    
    say "Registering $repeat_group $repeat (" ,
    $tot_reference_table->{ 'tot_n' } , ")";
    
    my $repeat_hit_length = 0;
    
    foreach my $hit_ref ( @{ $all_hits_ref } ) {
    
        my ( $ref_seq , $identity , $ref_start , $ref_stop ) = @{ $hit_ref };
    
        if ( $ref_start > $ref_stop ) {
        
            ( $ref_start , $ref_stop ) = ( $ref_stop , $ref_start );
        
        }
    
        my $hit_length = $ref_stop - $ref_start;
    
        my $ref_length = $reference_table->{ $ref_seq }->{ 'length' };
        
        unless ( defined $reference_table->{ $ref_seq }->{ 'seq' } ) {
        
            $reference_table->{ $ref_seq }->{ 'seq' } = 0 x $ref_length;
        
        }
        
        $reference_table->{ $ref_seq }->{ 'tot_hit' } += $hit_length;
        
        $tot_reference_table->{ 'tot_hit' } += $hit_length;
        
        substr (
            $reference_table->{ $ref_seq }->{ 'seq' } , 
            $ref_start - 1 ,
            $hit_length ,
            "1" x $hit_length ,
        );
        
        $tot_repeat_table->{ $repeat_group }->{ 'identity' }->
        [ int $identity ]++;
        
        $tot_repeat_table->{ $repeat_group }->{ 'hit_length' }->
        [ int ( $hit_length / 10 ) ]++;

        $tot_repeat_table->{ $repeat_group }->{ 'tot_hit_n' }++;
        $tot_repeat_table->{ $repeat_group }->{ 'tot_hit_length' } += $hit_length;
        
        $repeat_hit_length += $hit_length;
        
        $tot_reference_table->{ 'max_hit_length' } = $hit_length unless
        defined
        $tot_reference_table->{ 'max_hit_length' };
        
        if ( $hit_length > $tot_reference_table->{ 'max_hit_length' } ) {
        
            $tot_reference_table->{ 'max_hit_length' } = $hit_length;
        
        }
        
    }
    
    unless (
    
        defined
        $tot_repeat_table->{ $repeat_group }->{ 'seen' }->{ $repeat }
        
    ) {
    
        $tot_repeat_table->{ $repeat_group }->{ 'seen_n' }++;
        
        $tot_repeat_table->{ $repeat_group }->{ 'seen_length' } +=
        $repeat_table->{ $repeat_group }->{ $repeat }->{ 'length' };
        
        my $n_hits = scalar @{ $all_hits_ref };
        
        $tot_reference_table->{ 'max_n_hits' } = $n_hits unless
        defined
        $tot_reference_table->{ 'max_n_hits' };
        
        if ( $n_hits > $tot_reference_table->{ 'max_n_hits' } ) {
        
            $tot_reference_table->{ 'max_n_hits' } = $n_hits;
        
        }
        
        $tot_repeat_table->{ $repeat_group }->{ 'hit_n_dist' }->
        [ int $n_hits / 10 ]++;
    
        $repeat_table->{ $repeat_group }->{ $repeat }->{ 'hits_n' } = $n_hits;
        $repeat_table->{ $repeat_group }->{ $repeat }->{ 'hits_length' } = $repeat_hit_length;
    
        $tot_repeat_table->{ $repeat_group }->{ 'seen' }->{ $repeat } = 1;
    
    }

}

sub summarize_hits_dists {

    my ( $opt ) = @_;

    my $repeat_table = $opt->{ 'repeat_table' };
    my $tot_repeat_table = $opt->{ 'tot_repeat_table' };
    
    my $reference_table = $opt->{ 'reference_table' };
    my $tot_reference_table = $opt->{ 'tot_reference_table' };
    
    my $min_identity = $opt->{ 'min_identity' };
    
    # IDENTITY
    # --------
    
    my $out_file = $opt->{ 'output' } . ".identity.dist.csv";
    
    open ( my $dist_out , ">" , $out_file ) or die "$!";
    
    my @groups = sort keys %{ $tot_repeat_table };
    
    print $dist_out "IDENTITY";
    
    foreach my $repeat_group ( @groups ) {
    
        print $dist_out "\t${repeat_group}_N\t${repeat_group}_PROP";
    
    }
    
    say $dist_out "";
    
    for ( my $i = $min_identity; $i <= 100; $i++ ) {
    
        print $dist_out $i;
    
        foreach my $repeat_group ( keys %{ $tot_repeat_table } ) {
        
            $tot_repeat_table->{ $repeat_group }->{ 'identity' }->[ $i ] = 0
            unless defined
            $tot_repeat_table->{ $repeat_group }->{ 'identity' }->[ $i ];
            
            my $n =
            $tot_repeat_table->{ $repeat_group }->{ 'identity' }->[ $i ];
            
            my $prop = $n /
            $tot_repeat_table->{ $repeat_group }->{ 'tot_hit_n' };
            
            print $dist_out "\t$n\t$prop";
        
        }
        
        say $dist_out "";
    
    }
    
    # LENGTH
    # ------
    
    my $max_hit_length = $tot_reference_table->{ 'max_hit_length' };
    
    $out_file = $opt->{ 'output' } . ".length.dist.csv";
    
    open ( $dist_out , ">" , $out_file ) or die "$!";
    
    print $dist_out "IDENTITY";
    
    foreach my $repeat_group ( @groups ) {
    
        print $dist_out "\t${repeat_group}_N\t${repeat_group}_PROP";
    
    }
    
    say $dist_out "";
    
    for ( my $i = 0; $i <= $max_hit_length / 10 ; $i++) {
    
        print $dist_out $i * 10;
    
        foreach my $repeat_group ( keys %{ $tot_repeat_table } ) {
        
            $tot_repeat_table->{ $repeat_group }->{ 'hit_length' }->[ $i ] = 0
            unless defined
            $tot_repeat_table->{ $repeat_group }->{ 'hit_length' }->[ $i ];
            
            my $n =
            $tot_repeat_table->{ $repeat_group }->{ 'hit_length' }->[ $i ];
            
            my $prop = $n /
            $tot_repeat_table->{ $repeat_group }->{ 'tot_hit_n' };
            
            print $dist_out "\t$n\t$prop";
        
        }
        
        say $dist_out "";
    
    }

    # N HITS
    # ------
    
    my $max_n_hits = $tot_reference_table->{ 'max_n_hits' };
    
    $out_file = $opt->{ 'output' } . ".n_hits.dist.csv";
    
    open ( $dist_out , ">" , $out_file ) or die "$!";
    
    print $dist_out "N_HITS";
    
    foreach my $repeat_group ( @groups ) {
    
        print $dist_out "\t${repeat_group}_N\t${repeat_group}_PROP";
    
    }
    
    say $dist_out "";
    
    for ( my $i = 0; $i <= $max_n_hits / 10 ; $i++) {
    
        print $dist_out $i * 10;
    
        foreach my $repeat_group ( keys %{ $tot_repeat_table } ) {
        
            $tot_repeat_table->{ $repeat_group }->{ 'hit_n_dist' }->[ $i ] = 0
            unless defined
            $tot_repeat_table->{ $repeat_group }->{ 'hit_n_dist' }->[ $i ];
            
            my $n =
            $tot_repeat_table->{ $repeat_group }->{ 'hit_n_dist' }->[ $i ];
            
            my $prop = $n /
            $tot_repeat_table->{ $repeat_group }->{ 'seen_n' };
            
            print $dist_out "\t$n\t$prop";
        
        }
        
        say $dist_out "";
    
    }
    
}

sub summarize_hits_reference {

    my ( $opt ) = @_;

    my $repeat_table = $opt->{ 'repeat_table' };
    my $tot_repeat_table = $opt->{ 'tot_repeat_table' };
    
    my $reference_table = $opt->{ 'reference_table' };
    my $tot_reference_table = $opt->{ 'tot_reference_table' };

    my $out_file = $opt->{ 'output' } . ".hits_summary.csv";
    
    open ( my $out , ">" , $out_file ) or die "$!";
    
    my @groups = sort keys %{ $tot_repeat_table };
    
    say $out "REPEAT_GROUP\tN\tLENGTH\tAVG_LENGTH\t" ,
        "HIT_N\tHIT_LENGTH\tAVG_HIT_LENGTH\t" ,
        "AVG_HITS_PER_REPEAT\tAVG_PROP_LENGTH_PER_REPEAT\t" ,
        "REPEAT_AMPLIFICATION";
    
    my $tot_n = 0;
    my $tot_length = 0;
    
    my $tot_hit_n = 0;
    my $tot_hit_length = 0;
    
    foreach my $repeat_group ( @groups ) {
    
        my $n = $tot_repeat_table->{ $repeat_group }->{ 'seen_n' };
        my $length = $tot_repeat_table->{ $repeat_group }->{ 'seen_length' };
        
        $tot_n += $n;
        $tot_length += $length;
        
        my $avg_length = $length / $n;
        
        my $hit_n = $tot_repeat_table->{ $repeat_group }->{ 'tot_hit_n' };
        my $hit_length = $tot_repeat_table->{ $repeat_group }->{ 'tot_hit_length' };
        
        $tot_hit_n += $hit_n;
        $tot_hit_length += $hit_length;
        
        my $avg_hit_length = $hit_length / $hit_n;
        
        my $avg_hits_per_repeat = $hit_n / $n;
        my $avg_prop_length_per_repeat = $avg_hit_length / $avg_length;
        
        my $amplification = $hit_length / $length;
        
        say $out
            "$repeat_group\t$n\t$length\t$avg_length\t" ,
            "$hit_n\t$hit_length\t$avg_hit_length\t" ,
            "$avg_hits_per_repeat\t$avg_prop_length_per_repeat\t" ,
            $amplification;
    
    }
    
    my $avg_length = $tot_length / $tot_n;
    my $avg_hit_length = $tot_hit_length / $tot_hit_n;
    
    my $avg_hits_per_repeat = $tot_hit_n / $tot_n;
    my $avg_prop_length_per_repeat = $avg_hit_length / $avg_length;
    
    my $genome_length = $tot_reference_table->{ 'length' };
    my $masked = 0;
    
    foreach my $ref_seq ( keys %{ $reference_table } ) {
    
        next unless defined $reference_table->{ $ref_seq }->{ 'seq' };
    
        my $seq = $reference_table->{ $ref_seq }->{ 'seq' };
        
        $masked += $seq =~ tr/1//;
    
    }
    
    my $amplification = $tot_hit_length / $tot_length;
    
    my $masked_prop = $masked / $genome_length;
    
    my $masked_amplification = $masked / $tot_length;
    
    my $overlap = $masked_amplification / $amplification;

    $out_file = $opt->{ 'output' } . ".reference_summary.csv";
    
    open ( $out , ">" , $out_file ) or die "$!";
    
    say $out "N\tLENGTH\tAVG_LENGTH\t" ,
        "HIT_N\tHIT_LENGTH\tAVG_HIT_LENGTH\t" ,
        "AVG_HITS_PER_REPEAT\tAVG_PROP_LENGTH_PER_REPEAT\t" ,
        "REPEAT_AMPLIFICATION\tGENOME_LENGTH\tMASKED\tMASKED_PROP\tMASKED_AMPLIFICATION\tUNIQUENESS";
    
    say $out "$tot_n\t$tot_length\t$avg_length\t" ,
            "$tot_hit_n\t$tot_hit_length\t$avg_hit_length\t" ,
            "$avg_hits_per_repeat\t$avg_prop_length_per_repeat\t" ,
            "$amplification\t$genome_length\t$masked\t$masked_prop\t$masked_amplification\t" ,
            $overlap;

}

sub print_repeat_lists {

    my ( $opt ) = @_;
    
    my $tot_repeat_table = $opt->{ 'tot_repeat_table' };
    my $repeat_table = $opt->{ 'repeat_table' };
    
    my @groups = sort keys %{ $tot_repeat_table };

    foreach my $repeat_group ( @groups ) {
    
        my $out_file = $opt->{ 'output' } . ".valid_repeats.${repeat_group}.csv";

        open ( my $out , ">" , $out_file ) or die "$!";

        say $out "REPEAT\tHITS\tPROP_HITS\tLENGTH\tPROP_LENGTH\tCUMULATIVE_LENGTH";
        
        my $tot_n =  $tot_repeat_table->{ $repeat_group }->{ 'tot_hit_n' };
        my $tot_length = $tot_repeat_table->{ $repeat_group }->{ 'tot_hit_length' };
        
        my @seen_repeats = sort {
            $repeat_table->{ $repeat_group }->{ $b }->{ 'hits_length' } <=>
            $repeat_table->{ $repeat_group }->{ $a }->{ 'hits_length' }
        } keys %{ $tot_repeat_table->{ $repeat_group }->{ 'seen' } };
        
        my $tmp_prop;
        
        foreach my $repeat ( @seen_repeats ) {
        
            my $hits_n = $repeat_table->{ $repeat_group }->{ $repeat }->{ 'hits_n' };
            my $hits_length = $repeat_table->{ $repeat_group }->{ $repeat }->{ 'hits_length' };
        
            $tmp_prop += $hits_length / $tot_length;
        
            say $out $repeat , "\t" ,
                $hits_n , "\t" , $hits_n / $tot_n , "\t" ,
                $hits_length , "\t" , $hits_length / $tot_length , "\t" ,
                $tmp_prop;
        
        }

    }

}

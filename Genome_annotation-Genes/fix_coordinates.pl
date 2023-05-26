#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

while (<>) {

    chomp;
    
    # seq_s_24711	ALN	cds	158010	158008	16	-	2	ID=cds17333;Parent=mRNA02015;Name=seq_s_24711_184;Target=ESS050247 210 209 +
    
    if ( $_ =~ m/^(\S+\s\S+\s\S+)\s(\d+)\s(\d+)\s(\S.+)\s(\d+)\s(\d+)\s(\S+)/ ) {
    
        my ( $pre , $start , $stop , $post , $cds_start , $cds_stop , $final ) =
        ( $1 , $2 , $3 , $4 , $5 , $6 , $7 );
        
        if ( $start > $stop ) {
        
            ( $start , $stop ) = ( $stop , $start );
        
        }
        
        if ( $cds_start > $cds_stop ) {
        
            ( $cds_start , $cds_stop ) = ( $cds_stop , $cds_start );
        
        }
        
        say "$pre\t$start\t$stop\t$post $cds_start $cds_stop $final";
    
    }
    
    else {
    
        say $_;
    
    }

}

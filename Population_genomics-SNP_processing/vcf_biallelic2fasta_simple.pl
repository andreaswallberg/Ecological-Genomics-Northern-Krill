#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

use Getopt::Long;

# Option table and defaults

my $opt = {
    
    'input' => [] ,

};

GetOptions(
    $opt ,
    
    # Input
    # -----
    
    # VCF input (VCF file)
    'input|in|i=s{,}' ,
    
);

my $seq_table;

foreach my $vcf_file ( @{ $opt->{ 'input' } } ) {

    parse_vcf( $vcf_file );
    
}

my @samples;

foreach my $sample ( @samples ) {
	
	say ">${sample}\n" , $seq_table->{ $sample };
    
}

sub parse_vcf {

    my ( $vcf_file ) = @_;

    my $vcf_in;

    if ( $vcf_file =~ m/\.gz$/ ) {
    
        open ( $vcf_in , "zcat $vcf_file |" ) or die "$!";
    
    }
    
    else {

        open ( $vcf_in , "<" , $vcf_file ) or die "$!";

    }
    
    say STDERR "Reading VCF SNPs from $vcf_file ...";

	my $j = 0;

    while (<$vcf_in>) {

        chomp;
        
        my $line = $_;
        
        if ( $line =~ m/^\#\#/ ) {
        
			next;
        
        }
        
        if ( $line =~ m/^#CHROM.+\s+FORMAT\s+(.+)/ ) {
        
            @samples = split ( /\s+/ , $1 );
        
        }

		elsif ( $line =~ m/(\S+)\s+(\S+)\s+\S+\s+(\w)\s+(\w)\s+\S+.+\s+\S+\:(GP|PL|GL|PS)\s+(.+)/ ) {
			
				my ( $chrom , $pos , $ref , $alt , $gt_strings ) = ( $1 , $2 , $3 , $4 , $6 );
				
				$j++;

				say STDERR "$chrom\t$pos\t$j";
				
				my $i = 0;
			
				foreach my $genotype ( split ( /\s+/ , $gt_strings ) ) {
						
					my $sample = $samples[ $i ];
						
					if ( $genotype =~ m/^0(\/|\|)0\:(\d+)/ ) {
					
						$seq_table->{ $sample } .= "${ref}${ref}";
					
					}
					
					elsif (
					
						$genotype =~ m/^0(\/|\|)1\:(\d+)/ or
						$genotype =~ m/^1(\/|\|)0\:(\d+)/
					
					) {
					
						$seq_table->{ $sample } .= "${ref}${alt}";
					
					}
					
					elsif ( $genotype =~ m/^1(\/|\|)1\:(\d+)/ ) {
					
						$seq_table->{ $sample } .= "${alt}${alt}";
					
					}
					
					else {
					
						$seq_table->{ $sample } .= "--";
					
					}
					
					$i++;
				
				}
			
			}

		else {

			die "$_";

		}

    }

}

#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

use Getopt::Long;
use Data::Dumper;

# Option table and defaults

my $opt = {
    
    'input' => [] ,
    'output' => "filtered.fasta" ,
    
    'min_fill_sample' => 0 ,
    'min_fill_position' => 0 ,

    'min_depth' => 0 ,
    'max_depth' => undef ,
    'max_chrom' => undef ,
    
    'haploid' => 0 ,
    
    'skip' => [] ,
    
    'spacing' => 0 ,

};

GetOptions(
    $opt ,
    
    # Input
    # -----
    
    # VCF input (VCF file)
    'input|in|i=s{,}' ,
    
    'depth_file=s' ,
    
    'output=s' ,
    
    'min_fill_sample=s' ,
    'min_fill_position=s' ,
    'min_depth=s' ,
    'max_depth=s' ,
    'max_chrom=s' ,

    'spacing=i' ,
    
    'haploid' ,

    'skip=s{,}' , 
    
);

my $ctg_length_table = {};

if ( defined $opt->{ 'depth_file' } ) {

	open ( my $dp_in , "<" , $opt->{ 'depth_file' } ) or die "$!";
	
	my $min_dp = $opt->{ 'min_depth' };
	
	while (<$dp_in>) {
	
		chomp;
		
		# ctg1    40      59      56857   471265
		
		if ( $_ =~ m/^(\S+)\s+\d+\s+(\d+)\s+\d+\s+(\d+)/ ) {
		
			my ( $ctg , $dp , $len_ht ) = ( $1 , $2 , $3 );
			
			$ctg_length_table->{ $ctg }->{ 'seen' } = 1;
			
			unless ( defined $ctg_length_table->{ $ctg }->{ 'len_ht' } ) {
			
				if ( $dp >= $min_dp ) {
				
					$ctg_length_table->{ $ctg }->{ 'len_ht' } = $len_ht;
				
					# say $_;
				
				}
			
			}
			
		}
	
	}

}

my $skip_table = {};

foreach my $sample ( @{ $opt->{ 'skip' } } ) {

	$skip_table->{ $sample } = 1;

}

my @samples;
my @sorted_samples;

my $seq_table = {};
my $fill_table = {};
my $dp_table = {};

my $geno_table = {};

open (
	my $pos_all_out ,
	">" ,
	$opt->{ 'output' } . ".fasta.positions.all.csv"
) or die "$!";

open (
	my $pos_pas_out ,
	">" ,
	$opt->{ 'output' } . ".fasta.positions.passed.csv"
) or die "$!";

say $pos_all_out "CHROM\tPOS\tREF\tALT\tDP\tGTS\tFILL_RATE\t00\t01\t11\tMISSING";
say $pos_pas_out "CHROM\tPOS\tREF\tALT\tDP\tGTS\tFILL_RATE\t00\t01\t11\tMISSING";

my $sites = 0;

open ( my $vcf_out , ">" , $opt->{ 'output' } . ".fasta.vcf" ) or die "$!";

foreach my $vcf_file ( @{ $opt->{ 'input' } } ) {

    parse_vcf( $vcf_file );
    
}

my $min_fill_sample = $opt->{ 'min_fill_sample' };

open ( my $fasta_out , ">" , $opt->{ 'output' } . ".fasta" ) or die "$!";
open ( my $samples_out , ">" , $opt->{ 'output' } . ".fasta.samples.csv" ) or die "$!";

open ( my $geno_out , ">" , $opt->{ 'output' } . ".geno" ) or die "$!";
open ( my $geno_samples_out , ">" , $opt->{ 'output' } . ".geno.samples.csv" ) or die "$!";

my @sorted_filled_samples;

foreach my $sample ( @sorted_samples ) {

	next if defined $skip_table->{ $sample };

	$fill_table->{ $sample } = 0 unless defined $fill_table->{ $sample };
	
	my $fill_prop = $fill_table->{ $sample } / $sites;
	
	if ( $fill_prop > $min_fill_sample ) {
	
		push @sorted_filled_samples , $sample;
	
	}
	
}

say $geno_samples_out $_ foreach @sorted_filled_samples;

for ( my $i = 0; $i < $sites; $i++ ) {

	foreach my $sample ( @sorted_filled_samples ) {

		print $geno_out $geno_table->{ $sample }->[ $i ];

	}
	
	say $geno_out "";
	
}

say $samples_out "SAMPLE\tSITES\tDP\tFILL_RATE";

foreach my $sample ( @sorted_samples ) {

	unless ( defined $skip_table->{ $sample } ) {

		say $samples_out
			$sample , "\t" ,
			$sites , "\t" ,
			$dp_table->{ $sample } /
			$fill_table->{ $sample } , "\t" ,
			$fill_table->{ $sample } / $sites;
			
	}
    
}

foreach my $sample ( @sorted_filled_samples ) {
	
	say $fasta_out ">${sample}\n" , $seq_table->{ $sample };
    
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
    
    say "Reading VCF SNPs from $vcf_file ...";

    my $min_fill_position = $opt->{ 'min_fill_position' };
    my $min_depth = $opt->{ 'min_depth' };
    my $max_depth = $opt->{ 'max_depth' };
    my $max_chrom = $opt->{ 'max_chrom' };
    
    my $last_chrom = "";
    my $last_pos = 0;
    
    while (<$vcf_in>) {

        chomp;
        
        my $line = $_;
        
        if ( $line =~ m/^\#\#/ ) {
        
			say $vcf_out $line;
			next;
        
        }
        
        if ( $line =~ m/^#CHROM.+\s+FORMAT\s+(.+)/ ) {
        
			say $vcf_out $line;
        
            @samples = split ( /\s+/ , $1 );
            @sorted_samples = sort @samples;
        
# 			foreach my $sample ( @samples ) {
#         
# 				$geno_table->{ $sample } = [];
#         
# 			}
        
        }
        
        # elsif ( $_ =~ m/(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+\S+.+DP\:RO\:QR\:AO\:QA\:GL\s+(.+)/ ) {

		elsif ( $line =~ m/(\S+)\s+(\S+)\s+\S+\s+(\w)\s+(\w)\s+\S+.+\;DP=(\d+).+\:(PL|GL)\s+(.+)/ ) {
		
            my ( $chrom , $pos , $ref , $alt , $dp , $gt_strings ) = ( $1 , $2 , $3 , $4 , $5 , $7 );           
			
			if ( $chrom eq "seq_s_3" ) {
			
				last;
			
			}
			
			if ( defined $opt->{ 'depth_file' } and not defined $ctg_length_table->{ $chrom }->{ 'seen' } ) {
			
				next;
			
			}
			
			# say "$chrom\t$pos";
			
			next if length $ref != 1;
			next if length $alt != 1;
			
			next if $dp < $min_depth;
			
			if ( defined $max_depth ) {
			
				next if $dp > $max_depth;
			
			}
			
			if ( defined $max_chrom ) {
			
				my $chrom_n = $chrom;
				$chrom_n =~ s/^ctg//;
				$chrom_n =~ s/(\-|\.).+$//;
				
				if ( $chrom_n > $max_chrom ) {
				
					next;
				
				}
			
			}
			
			if ( $chrom eq $last_chrom ) {
			
				if ( $pos - $last_pos < $opt->{ 'spacing' } ) {
				
					next;
				
				}
			
			}
			
			else {
			
				say $chrom;
			
				$last_chrom = "";
				$last_pos = 0;
			
			}
			
			my $data_table = {};
			my $tmp_dp_table = {};
			my $tmp_geno_table = {};

			my $gts = 0;
			my $gts_00 = 0;
			my $gts_01 = 0;
			my $gts_11 = 0;
			
			my $missing = 0;
			
            my $i = 0;
            
            foreach my $genotype ( split ( /\s+/ , $gt_strings ) ) {
                    
                my $sample = $samples[ $i ];
                
                unless ( defined $skip_table->{ $sample } ) {
                
					if ( $opt->{ 'haploid' } ) {
					
						if ( $genotype =~ m/^0\:(\d+)/ ) {
						
							$data_table->{ $sample } = $ref;
							$tmp_dp_table->{ $sample } += $1;
							$tmp_geno_table->{ $sample } = 1;
							
							$gts++;
							$gts_00++;
						
						}
						
						elsif ( $genotype =~ m/^1\:(\d+)/ ) {
						
							$data_table->{ $sample } = $alt;
							$tmp_dp_table->{ $sample } += $1;
							$tmp_geno_table->{ $sample } = 0;
							
							$gts++;
							$gts_11++;
						
						}
						
						else {
						
							$data_table->{ $sample } = "?";
							$tmp_geno_table->{ $sample } = 9;
							
							$missing++;
						
						}
					
					}
					
					else {
					
						if ( $genotype =~ m/^0(\/|\|)0\:(\d+)/ ) {
						
							$data_table->{ $sample } = "${ref}${ref}";
							$tmp_dp_table->{ $sample } = $2;
							$tmp_geno_table->{ $sample } = 2;
							
							$gts++;
							$gts_00++;
						
						}
						
						elsif (
						
							$genotype =~ m/^0(\/|\|)1\:(\d+)/ or
							$genotype =~ m/^1(\/|\|)0\:(\d+)/
						
						) {
						
							$data_table->{ $sample } = "${ref}${alt}";
							$tmp_dp_table->{ $sample } = $2;
							$tmp_geno_table->{ $sample } = 1;
							
							$gts++;
							$gts_01++;
						
						}
						
						elsif ( $genotype =~ m/^1(\/|\|)1\:(\d+)/ ) {
						
							$data_table->{ $sample } = "${alt}${alt}";
							$tmp_dp_table->{ $sample } = $2;
							$tmp_geno_table->{ $sample } = 0;
							
							$gts++;
							$gts_11++;
						
						}
						
						else {
						
							$data_table->{ $sample } = "??";
							$tmp_geno_table->{ $sample } = 9;
						
							$missing++;
						
						}
					
					}
                
                }
                
                $i++;
            
            }
            
            if ( ( $gts_00 > 0 and $gts_11 > 0 ) or $gts_01 > 0 ) {
            
				say $pos_all_out
					"$chrom\t$pos\t$ref\t$alt\t$dp\t" ,
					"$gts\t" , $gts / ( $gts + $missing ) , "\t" ,
					"$gts_00\t$gts_01\t$gts_11\t$missing";
					
				if ( ( $gts / ( $gts + $missing ) ) > $min_fill_position ) {
				
					say $vcf_out $line;
				
					say $pos_pas_out
						"$chrom\t$pos\t$ref\t$alt\t$dp\t" ,
						"$gts\t" , $gts / ( $gts + $missing ) , "\t" ,
						"$gts_00\t$gts_01\t$gts_11\t$missing";
				
					foreach my $sample ( @sorted_samples ) {
					
						unless ( defined $skip_table->{ $sample } ) {
					
							unless ( defined $data_table->{ $sample } ) {
							
								say "$chrom\t$pos\t$sample\t$_";
							
							}
						
							if ( $data_table->{ $sample } =~ m/\w/ ) {
							
								$fill_table->{ $sample }++;
							
								$dp_table->{ $sample } +=
								$tmp_dp_table->{ $sample };
							
							}
						
							$seq_table->{ $sample } .=
							$data_table->{ $sample };
							
							$geno_table->{ $sample }->[ $sites ] =
							$tmp_geno_table->{ $sample };
					
						}
					
					}

					$sites++;
					
					# say $sites;
					
					# last if $sites == 100000;
					
					$last_pos = $pos;
				
				}
			
			}
			
			$last_chrom = $chrom;
        
        }

    }

}

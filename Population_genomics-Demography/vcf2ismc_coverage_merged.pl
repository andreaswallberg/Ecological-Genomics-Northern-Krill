#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

use Getopt::Long;
use Data::Dumper;

# Option table and defaults

my $opt = {
    
    # Input files
    
    'vcf' => [] ,
	'files' => [] ,
    
    'group' => {} ,
    
};

GetOptions(
    $opt ,

    'files=s{,}' ,
    
    'sample=s' ,
    
    'coverage|c=s' ,
    
    'seqs|s=s' ,
    
    'fasta=s' ,
    
    'output|out|o=s' ,
   
    # Booleans
    
    'verbose|v' ,
    'help|h' ,

);

# Read coverage data
# ------------------

init_coverage( $opt );

# Parse phased data in vcf format
# ----------------------------------

my $seen_table = {};

$opt->{ 'seen_table' } = $seen_table;
	
my $seen_header = 0;

my $vcf_table = {};

my @headers;
my @chromosomes;

init_files( $opt );
init_vcf( $opt );

sub init_coverage {
    
    my ( $opt ) = @_;
    
    my $coverage_table = {};
    my $seq_table = {};
    
    $opt->{ 'coverage_table' } = $coverage_table;
    $opt->{ 'seq_table' } = $seq_table;
    
    if ( defined $opt->{ 'seqs' } and -e $opt->{ 'seqs' } ) {
    
        open ( my $in , "<" , $opt->{ 'seqs' } ) or die "$!";
    
        say "Reading sequences from $opt->{ 'seqs' } ...";
    
        while (<$in>) {
        
            chomp;
            
            if ( $_ =~ m/^(\S+)\s\d+\s(\d+)/ ) {
            
                unless ( defined $opt->{ 'get_chrom' }->{ $1 } ) {
            
                    $opt->{ 'get_chrom' }->{ $1 } = $2 + 1;
                    
                    push @{ $opt->{ 'get_chroms' } } , $1;
            
					say $1;
            
                }
            
            }
        
        }
    
    }
    
    if ( defined $opt->{ 'coverage' } and -e $opt->{ 'coverage' } ) {
        
        my $coverage_file = $opt->{ 'coverage' };
        
        say "Reading coverage stats from $coverage_file ...";

        my $coverage_in;
        
        if ( $coverage_file =~ m/\.gz/ ) {
        
            open ( $coverage_in , "zcat $coverage_file |" ) or die "$!";
        
        }
        
        else {
        
            open ( $coverage_in , "<" , $coverage_file ) or die "$!";
        
        }
        
        my $tmp_header = undef;
        
        while (<$coverage_in>) {

			if ( $_ =~ m/^\>(\S+)/ ) {
			
				$tmp_header = undef;
			
				if ( defined $opt->{ 'get_chrom' }->{ $1 } ) {
			
					$tmp_header = $1;
			
				}
			
			}
			
			elsif ( defined $tmp_header and defined $opt->{ 'get_chrom' }->{ $tmp_header } ) {
				
				$coverage_table->{ $tmp_header } .= $_;
			
			}
			
        }
        
    }
    
    if ( defined $opt->{ 'fasta' } and -e $opt->{ 'fasta' } ) {
        
        my $seq_file = $opt->{ 'fasta' };
        
        say "Reading FASTA seq from $seq_file ...";

        my $seq_in;
        
        if ( $seq_file =~ m/\.gz/ ) {
        
            open ( $seq_in , "zcat $seq_file |" ) or die "$!";
        
        }
        
        else {
        
            open ( $seq_in , "<" , $seq_file ) or die "$!";
        
        }
        
        my $tmp_header = undef;
        
        while (<$seq_in>) {

			chomp;
        
			if ( $_ =~ m/^\>(\S+)/ ) {
			
				$tmp_header = undef;
			
				if ( defined $opt->{ 'get_chrom' }->{ $1 } ) {
			
					$tmp_header = $1;
			
				}
			
			}
			
			elsif ( defined $tmp_header and defined $opt->{ 'get_chrom' }->{ $tmp_header } ) {
				
				$seq_table->{ $tmp_header } .= $_;
			
			}
			
        }
        
    }
    
}

sub init_files {

	my ( $opt ) = @_;

	my $file_table = {};
	
	my @files;
	
	foreach my $file ( @{ $opt->{ 'files' } } ) {
	
		open ( my $in , "<" , $file ) or die "$!";
	
		while (<$in>) {
		
			chomp;
			
			my ( $tmp_header , $file ) = split ( /\t/ , $_ );
		
			if ( defined $opt->{ 'get_chrom' }->{ $tmp_header } ) {
			
				unless ( defined $file_table->{ $file } ) {
			
					$file_table->{ $file } = 1;
			
					push @files , $file;
			
				}
			
			}
		
		}
	
	}
	
	foreach my $file ( @files ) {
	
		say "Will read SNPs from $file ...";
	
		push @{ $opt->{ 'vcf' } } , $file;
	
	}
	
}

sub init_vcf {
    
    my ( $opt ) = @_;
       
    # Read VCF
    # --------
    
    foreach my $vcf_file ( @{ $opt->{ 'vcf' } } ) {
        
        parse_vcf( $opt , $vcf_file );
        
    }
    
    my $seq_table = $opt->{ 'seq_table' };
    my $sample = $opt->{ 'sample' };
    
    my $coverage_table = $opt->{ 'coverage_table' };
    
    if ( @chromosomes ) {
    
		open ( my $out , ">" , $opt->{ 'output' } . ".vcf" ) or die "$!";
		
		say $out $_ foreach @headers;

		open ( my $fasta_out , ">" , $opt->{ 'output' } . ".fasta" ) or die "$!";
		open ( my $tab_out , ">" , $opt->{ 'output' } . ".tab" ) or die "$!";
    
		foreach my $chrom ( @chromosomes ) {
		
			say "Generating output for $chrom ...";
			
			my $coverage_seq = $coverage_table->{ $chrom };
			my $seq = $seq_table->{ $chrom };
			
			my $length_coverage_seq = length $coverage_seq;
			
			for ( my $i = 0; $i < $length_coverage_seq; $i++ ) {
			
				my $cov = substr ( $coverage_seq , $i , 1 );
				
				unless ( $cov ) {
				
					substr ( $seq , $i , 1 , "N" );
				
				}
			
			}
			
			my $max_length = ( int ( $length_coverage_seq / 100_000 ) ) * 100_000;
			
			my @positions = sort { $a <=> $b } keys %{ $vcf_table->{ $chrom } };
			
			foreach my $pos ( @positions ) {
			
				say $out $vcf_table->{ $chrom }->{ $pos };
			
			}
			
			say $tab_out "$chrom\t1\t$max_length\t0\t" , $max_length - 1 , "\t1\t$max_length";
			
			say $fasta_out ">$chrom\n" , substr ( $seq , 0 , $max_length );
			
		}
    
    }
    
}

sub parse_vcf {
    
    my ( $opt , $vcf_file ) = @_;
    
    my $coverage_table = $opt->{ 'coverage_table' };
    
    my $window = $opt->{ 'window' };
    
    my $seq_table = $opt->{ 'seq_table' };
    
    if ( -e $vcf_file ) {
        
        my $vcf_in;
        
        if ( $vcf_file =~ m/\.gz/ ) {
        
            open ( $vcf_in , "zcat $vcf_file |" ) or die "$!";
        
        }
        
        else {
        
            open ( $vcf_in , "<" , $vcf_file ) or die "$!";
        
        }
        
        if ( $opt->{ 'verbose' } ) {
            
            say "==> Reading from $vcf_file";
            
        }
        
		my $sample_id;
        
        while (<$vcf_in>) {
            
            chomp;
            
            if ( $_ =~ m/^\#\#/ ) {
            
                unless ( $seen_header ) {
                
					push @headers , $_;
                
                }
            
            }
            
            elsif ( $_ =~ m/^(#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT)\s+(.+)/ ) {
            
				my $i = 0;
							
				foreach my $sample ( split ( /\t/ , $2 ) ) {
				
					if ( $sample eq $opt->{ 'sample' } ) {
						
						$sample_id = $i;
						
						unless ( $seen_header ) {
						
							push @headers , "$1\t$sample";
							
							$seen_header = 1;
						
						}
						
						last;
				
					}
				
					$i++;
				
				}
				
				last;
                
            }
            
        }
        
        my $last_chrom = "";
        
        my $tmp_coverage_seq;
        
        while (<$vcf_in>) {
            
            chomp;
            
            if ( $_ =~ m/^(\S+)/ ) {
            
				unless ( defined $coverage_table->{ $1 } ) {
				
					next;
				
				}
            
            }
            
            if ( $_ =~m/^(\S+)\s+(\d+)\s+(\S+\s+\S\s+\S\s+\S+\s+.+\S+\s+GT\S*)\s+(.+)/ ) {
                
                my ( $chrom , $pos , $info , $data ) = ( $1  , $2 , $3 , $4 );
                
                if ( $chrom ne $last_chrom ) {
                
					say STDERR "Checking chromosome $chrom in $vcf_file ...";
                
                }
                
				if ( defined $coverage_table->{ $chrom } ) {
				
					if ( $chrom ne $last_chrom ) {
				
						$tmp_coverage_seq = $coverage_table->{ $chrom };
				
					}
					
					my $coverage_base = substr ( $tmp_coverage_seq , $pos - 1 , 1 );
					
					if ( $coverage_base != 0 ) {
				
						my @gts = split ( /\t/ , $data );
						
						my $gt = $gts[ $sample_id ];
						
						$info =~ s/\;ANN\=[^\|]+\|.+\s/\;\t/;
						
						$gt =~ s/\|/\//;
						
						if ( $gt =~ m/^0\W0/ ) {
						
							$vcf_table->{ $chrom }->{ $pos } = "$chrom\t$pos\t$info\t$gt";
						
						}
						
						else {
						
							if ( $gt =~ m/^0\W1/ ) {
							
								push @chromosomes , $chrom unless defined $seen_table->{ $chrom };
							
								$vcf_table->{ $chrom }->{ $pos } = "$chrom\t$pos\t$info\t$gt";
							
								$seen_table->{ $chrom } = 1;
							
							}
							
							elsif ( $gt =~ m/^(1\W0)(.+)/ ) {
							
								push @chromosomes , $chrom unless defined $seen_table->{ $chrom };
							
								$vcf_table->{ $chrom }->{ $pos } = "$chrom\t$pos\t$info\t0/1${2}";
							
								$seen_table->{ $chrom } = 1;
							
							}
							
							elsif ( $gt =~ m/^1\W1/ ) {
							
								$vcf_table->{ $chrom }->{ $pos } = "$chrom\t$pos\t$info\t$gt";
							
							}
						
						}
					
					}
					
					$last_chrom = $chrom;
                
                }
        
            }

        }
            
        close $vcf_in;
        
    }
    
    else {
        
        die "Unable to locate VCF formatted input file $vcf_file: $!";
        
    }
    
}

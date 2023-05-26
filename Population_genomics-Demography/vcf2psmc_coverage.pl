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

    'group' => {} ,
    
    'window' => 0 ,
    
};

GetOptions(
    $opt ,

    'sample=s' ,
    
    'coverage|c=s' ,
    
    'seqs|s=s' ,
    
    'vcf|input=s{,}' ,
    
    'window|windows=i' ,
    
    'output|out|o=s' ,
   
    # Booleans
    
    'verbose|v' ,
    'help|h' ,

);

# Possibly print help
# -------------------

#pod2usage( { -exitval => 1, -verbose => 2, -output => \*STDOUT } ) if not
#scalar @ARGV;

=head1 SYNOPSIS

 watterson_pi.pl [options]

 Options and switches:

=cut

# Possibly print help
# -------------------

pod2usage( { -exitval => 1, -verbose => 2, -output => \*STDOUT } ) if
$opt->{ 'help' };

# Set defaults
# ------------

init_defaults( $opt );

# Read coverage data
# ------------------

init_coverage( $opt );

# Parse phased data in vcf format
# ----------------------------------

init_vcf( $opt );

print_seq( $opt );

sub init_defaults {
    
    my ( $opt ) = @_;
    
    if ( defined $opt->{ 'vcf' } and not defined $opt->{ 'output' } ) {
        
        $opt->{ 'output' } = $opt->{ 'vcf' }->[ 0 ];
        
    }
    
}

sub init_coverage {
    
    my ( $opt ) = @_;
    
    my $coverage_table = {};
    
    $opt->{ 'coverage_table' } = $coverage_table;
    
    if ( defined $opt->{ 'seqs' } and -e $opt->{ 'seqs' } ) {
    
        open ( my $in , "<" , $opt->{ 'seqs' } ) or die "$!";
    
        say "Reading sequences from $opt->{ 'seqs' } ...";
    
        while (<$in>) {
        
            chomp;
            
            if ( $_ =~ m/^(\S+)/ ) {
            
                unless ( defined $opt->{ 'get_chrom' }->{ $1 } ) {
            
                    $opt->{ 'get_chrom' }->{ $1 } = 1;
                    
                    push @{ $opt->{ 'get_chroms' } } , $1;
            
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
        
			chomp;

			if ( $_ =~ m/^\>(.+)/ ) {
			
                $tmp_header = undef;
			
                if ( defined $opt->{ 'get_chrom' }->{ $1 } ) {
                
                    $tmp_header = $1;
                
                }
			
			}
			
			elsif ( defined $tmp_header ) {
			
                $coverage_table->{ $tmp_header } .= $_;
			
			}
			
        }
        
    }
    
}

sub init_vcf {
    
    my ( $opt ) = @_;
    
    my $length_table = {};
    
    $opt->{ 'length_all_table' } = $length_table;
    
    my $sample_seq_table = {};
    
    $opt->{ 'sample_seq_table' } = $sample_seq_table;
    
    my $sample_win_table = {};
    
    $opt->{ 'sample_win_table' } = $sample_win_table;
       
    # Read VCF
    # --------
    
    foreach my $vcf_file ( @{ $opt->{ 'vcf' } } ) {
        
        parse_vcf( $opt , $vcf_file );
        
    }
    
}

sub parse_vcf {
    
    my ( $opt , $vcf_file ) = @_;

    my $length_table = $opt->{ 'length_all_table' };  
    my $sample_seq_table = $opt->{ 'sample_seq_table' };
    my $sample_win_table = $opt->{ 'sample_win_table' };
    
    my $coverage_table = $opt->{ 'coverage_table' };
    
    my $window = $opt->{ 'window' };
    
    my $seen_table = {};
    
    $opt->{ 'seen_table' } = $seen_table;
        
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
            
                next;
            
            }
            
            elsif ( $_ =~ m/^#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT\s+(.+)/ ) {

                my $i = 0;
                               
                foreach my $sample ( split ( /\t/ , $1 ) ) {
                
                    if ( $sample eq $opt->{ 'sample' } ) {
                        
                        $sample_id = $i;
                        
                        last;
                
                    }
                
					$i++;
                
                }
                
                last;
                
            }
            
        }
        
        my $last_chrom = "";
        
        while (<$vcf_in>) {
            
            chomp;
            
            if ( $_ =~ m/^(\S+)/ ) {
            
				unless ( defined $coverage_table->{ $1 } ) {
				
					next;
				
				}
            
            }
            
            if ( $_ =~m/^(\S+)\s+(\d+)\s+\S+\s+(\S)\s+(\S)\s+\S+\s+.+\S+\s+GT\S*\s+(.+)/ ) {
                
                my ( $chrom , $pos , $ref , $alt , $data ) = ( $1  , $2 , $3 , $4 , $5 );
                
				my $window_pos = int ( ( $pos - 1 ) / $opt->{ 'window' } );
                
                if ( $chrom ne $last_chrom ) {
                
					say STDERR "Checking chromosome $chrom in $vcf_file ...";
                
                }
                
				if ( defined $coverage_table->{ $chrom } ) {

					my @gts = split ( /\t/ , $data );
					
					my $gt = $gts[ $sample_id ];
					
					if ( $gt =~ m/^0\W1/ ) {
					
						substr ( $coverage_table->{ $chrom } , $window_pos , 1 , "K" );
					
						$seen_table->{ $chrom } = 1;
					
					}
					
					elsif ( $gt =~ m/^1\W0/ ) {

						substr ( $coverage_table->{ $chrom } , $window_pos , 1 , "K" );
					
						$seen_table->{ $chrom } = 1;
					
					}
                
                }
                
                $last_chrom = $chrom;
        
            }

        }
            
        close $vcf_in;
        
    }
    
    else {
        
        die "Unable to locate VCF formatted input file $vcf_file: $!";
        
    }
    
}

sub print_seq {
    
    my ( $opt ) = @_;
    
    my $seen_table = $opt->{ 'seen_table' };
	my $coverage_table = $opt->{ 'coverage_table' };
	
    open ( my $out , ">" , $opt->{ 'output' } ) or die "$!";
    
	foreach my $chrom ( @{ $opt->{ 'get_chroms' } } ) {
    
		next unless defined $seen_table->{ $chrom };
        next unless defined $coverage_table->{ $chrom };
        
        my $seq = $coverage_table->{ $chrom };
        
        say $out ">$chrom";
        
        while ( $seq ) {
        
			my $seq_length = length $seq;
			
			if ( $seq_length < 60 ) {
			
				my $subseq = substr ( $seq , 0 , $seq_length , "" );
			
				say $out $subseq;
			
			}
			
			else {
        
				my $subseq = substr ( $seq , 0 , 60 , "" );
			
				say $out $subseq;
        
			}
        
        }
        
    }
    
}

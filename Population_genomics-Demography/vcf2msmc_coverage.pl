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

    # Groups of samples

    'groups' => [] ,

    'group' => {} ,
    
    'window' => 0 ,
    
};

GetOptions(
    $opt ,

    'group|groups|g=s{,}' ,
    
    'coverage|c=s' ,
    
    'seqs|s=s' ,
    
    'vcf|input=s{,}' ,
    
    'window|windows=i' ,
    
    'output|out|o=s' ,
   
    # Booleans
    
    'verbose|v' ,
    'help|h' ,

);

# Set defaults
# ------------

init_defaults( $opt );

# Create groups
# -------------

init_groups( $opt );

# Read coverage data
# ------------------

init_coverage( $opt );

# Parse phased data in vcf format
# ----------------------------------

init_vcf( $opt );

sub init_defaults {
    
    my ( $opt ) = @_;
    
    if ( defined $opt->{ 'vcf' } and not defined $opt->{ 'output' } ) {
        
        $opt->{ 'output' } = $opt->{ 'vcf' }->[ 0 ];
        
    }
    
}

sub init_groups {
    
    my ( $opt ) = @_;

    foreach my $group ( keys %{ $opt->{ 'group' } } ) {
    
        # Parse and generate groups from the CSV file
        
        my $group_file = $opt->{ 'group' }->{ $group };

        open ( my $in , "<" , $group_file ) or die "$!";
        
        say "Reading group file $group_file ..." if $opt->{ 'verbose' };
        
        while (<$in>) {
        
            chomp;
            
            my $sample = $_;
            
            $opt->{ 'app' }->{ 'sample_table' }->{ $sample } = $group;

            push @{
                $opt->{ 'app' }->{ 'group_table' }->{ $group }->{ 'samples' }
            } , $sample;
        
        }
	
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
        
        my @sample_names;
        
        my $group_table = {};
        
        $opt->{ 'group_table' } = $group_table;
        
        my @groups;
        
        while (<$vcf_in>) {
            
            chomp;
            
            if ( $_ =~ m/^\#\#/ ) {
            
                next;
            
            }
            
            elsif ( $_ =~ m/^#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT\s+(.+)/ ) {

                @sample_names = ( split ( /\s+/ , $1 ) );

                my $i = 0;
                
                my @these_samples;
                
                foreach my $sample ( @sample_names ) {
                
                    if ( defined $opt->{ 'app' }->{ 'sample_table' }->{ $sample } ) {
                
                        my $group =
                        $opt->{ 'app' }->{ 'sample_table' }->{ $sample };
                        
                        push @{ $group_table->{ $group }->{ 'ids' } } , $i;
                
						push @these_samples , ( $sample , $i );
                
                    }
                
					$i++;
                
                }
                
                say STDERR "Samples in this order: @these_samples";
                # exit;
                
                @groups = sort keys %{ $group_table };
                
                last;
                
            }
            
        }
        
        my $last_chrom = "";
        
        my $coverage_seq;
        
        while (<$vcf_in>) {
            
            chomp;
            
            if ( $_ =~ m/^(\S+)/ ) {
            
				unless ( defined $coverage_table->{ $1 } ) {
				
					next;
				
				}
            
            }
            
            if ( $_ =~m/^(\S+)\s+(\d+)\s+\S+\s+(\S)\s+(\S)\s+\S+\s+.+\S+\s+GT\S*\s+(.+)/ ) {
                
                my ( $chrom , $pos , $ref , $alt , $data ) = ( $1  , $2 , $3 , $4 , $5 );
                
                # say "$chrom\t$pos";
                
                # last if $pos > 100_000;
                
                if ( $chrom ne $last_chrom ) {
                
					say STDERR "Checking chromosome $chrom in $vcf_file ...";
                
                }
                
				if ( defined $coverage_table->{ $chrom } ) {
				
					unless ( $chrom eq $last_chrom ) {
				
						$coverage_seq = $coverage_table->{ $chrom };
						
						foreach my $group ( @groups ) {
						
							$group_table->{ $group }->{ 'last_pos' } = 0;
						
							open ( my $out , ">" , $opt->{ 'output' } . ".msmc.${chrom}.${group}" ) or die "$!";
						
							$group_table->{ $group }->{ 'out' } = $out;
						
						}
						
					}
					
					# last if $chrom eq "seq_s_2";
					
					# say "$chrom\t$pos";
				
					my $this_pos = $pos;
					$this_pos--;
					
					my $covered_base =
					substr ( $coverage_seq , $this_pos , 1 );
					
					next unless $covered_base;
					
					my @gts = split ( /\t/ , $data );
					
					my $seen = 0;
					
					foreach my $group ( @groups ) {
						
						my @group_gts = @gts[
							@{ $group_table->{ $group }->{ 'ids' } }
						];

						my $group_allele_table = {};

						my $i = 0;
						
						my $dna_string;
						
						foreach my $gt ( @group_gts ) {
						
							if ( $gt =~ m/^0\W0/ ) {
							
								$group_allele_table->{ $ref } += 2;
								
								$dna_string .= "${ref}${ref}";
							
							}
							
							elsif ( $gt =~ m/^1\W1/ ) {
							
								$group_allele_table->{ $alt } += 2;
								
								$dna_string .= "${alt}${alt}";
							
							}
							
							elsif ( $gt =~ m/^0\W1/ ) {
							
								$group_allele_table->{ $ref }++;
								$group_allele_table->{ $alt }++;
								
								$dna_string .= "${ref}${alt}";
							
							}
							
							elsif ( $gt =~ m/^1\W0/ ) {
							
								$group_allele_table->{ $ref }++;
								$group_allele_table->{ $alt }++;
								
								$dna_string .= "${alt}${ref}";
							
							}
						
						}
						
						# Keep track on variable SNPs in the group
						
						if ( keys %{ $group_allele_table } == 2 ) {
						
							my $last_pos = $group_table->{ $group }->{ 'last_pos' };
							
							my $length = $this_pos - $last_pos;
							
							my $subseq =
							substr ( $coverage_seq , $last_pos + 1 , $length );
							
							my $homozygous_sites = $subseq =~ tr/1//;
							
							my $out = $group_table->{ $group }->{ 'out' };
							
							say $out "$chrom\t$pos\t$homozygous_sites\t$dna_string";
							
							$group_table->{ $group }->{ 'last_pos' } = $this_pos;
						
						}
						
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

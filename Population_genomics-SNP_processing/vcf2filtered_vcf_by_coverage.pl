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
    
};

GetOptions(
    $opt ,
    
    'vcf=s{,}' ,
    
    'coverage|c=s' ,
    
    'seqs|s=s' ,
    
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
           
    # Read VCF
    # --------
    
    foreach my $vcf_file ( @{ $opt->{ 'vcf' } } ) {
        
        parse_vcf( $opt , $vcf_file );
        
    }
    
}

sub parse_vcf {
    
    my ( $opt , $vcf_file ) = @_;
    
    my $coverage_table = $opt->{ 'coverage_table' };
    
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
        
        open ( my $vcf_keep_out , ">" , $opt->{ 'output' } . ".keep.vcf" ) or die "$!";
        
        open ( my $vcf_remove_out , ">" , $opt->{ 'output' } . ".removed.vcf" ) or die "$!";
        
        open ( my $vcf_discarded_out , ">" , $opt->{ 'output' } . ".discarded.vcf" ) or die "$!";
        
        my $coverage_seq;
        my $last_chrom = "";
        
        while (<$vcf_in>) {
            
            chomp;
            
            if ( $_ =~ m/^\#/ ) {
            
                say $vcf_keep_out $_;
                say $vcf_remove_out $_;
                say $vcf_discarded_out $_;
                
            }
            
            elsif ( $_ =~m/^(\S+)\s+(\d+)\s\S+\s\w\s\w\s/ ) {
                
                my ( $chrom , $pos ) = ( $1 , $2 );
                
                if ( $chrom ne $last_chrom ) {
                
                    next unless defined $coverage_table->{ $chrom };
                
                    $coverage_seq = $coverage_table->{ $chrom };
                    
                    say "$chrom\t$pos";
                
                }
                
                my $covered_base =
                substr ( $coverage_seq , $pos - 1 , 1 );
                
                if ( $covered_base ) {
                
                    say $vcf_keep_out $_;
                
                }
                
                else {
                
                    say $vcf_remove_out $_;
                
                }
                
                $last_chrom = $chrom;

            }
            
            else {
            
                say $vcf_discarded_out $_;
            
            }
            
        }
        
    }
    
}

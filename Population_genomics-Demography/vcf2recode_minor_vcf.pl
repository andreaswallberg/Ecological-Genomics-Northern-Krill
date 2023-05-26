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
    
	'files' => [] ,
        
};

GetOptions(
    $opt ,

	'minor=s' ,
    'tag=s' ,
    'files=s{,}' ,
    
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

# Read minor allele data
# ----------------------

init_minor( $opt );
init_files( $opt );

sub init_minor {
    
    my ( $opt ) = @_;
    
    my $i = 1;
    
    if ( defined $opt->{ 'minor' } and -e $opt->{ 'minor' } ) {
    
        open ( my $in , "<" , $opt->{ 'minor' } ) or die "$!";
    
        say "Reading sequences from $opt->{ 'minor' } ...";
    
		my $header = <$in>;
    
        while (<$in>) {
        
            chomp;
            
            if ( $_ =~ m/^(\S+)\s(\d+)\s\S\s(\S)/ ) {
            
				unless ( defined $opt->{ 'minor_table' }->{ $1 } ) {
				
					say $i;
				
				}
            
                $opt->{ 'minor_table' }->{ $1 }->{ $2 } = $3;
            
				$i++;
            
            }
        
        }
    
    }
    
}

sub init_files {

	my ( $opt ) = @_;
	
	foreach my $file ( @{ $opt->{ 'files' } } ) {
	
		open ( my $in , "<" , $file ) or die "$!";
	
		while (<$in>) {
		
			chomp;
			
			my ( $file ) = split ( /\t/ , $_ );
			
			if ( -e $file ) {
			
				say "Reading $file ...";
		
				parse_vcf( $opt , $file );
		
			}
			
			else {
			
				say "No such file $file. Skipping!";
			
			}
		
		}
	
	}
	
}

sub parse_vcf {
    
    my ( $opt , $vcf_file ) = @_;
    
    my $tag = $opt->{ 'tag' };
    
    my $out_file = $vcf_file . ".${tag}.vcf.gz";
	
	if ( -e $out_file ) {
	
		say "$out_file already exists. Skipping!";
		return;
	
	}
	
    open ( my $out , "| gzip -c - > $out_file" ) or die "$!";
    
    say "Making $out_file ...";
    
    my $minor_table = $opt->{ 'minor_table' };
    
    if ( -e $vcf_file ) {
        
        my $vcf_in;
        
        if ( $vcf_file =~ m/\.gz/ ) {
        
            open ( $vcf_in , "zcat $vcf_file |" ) or die "$!";
        
        }
        
        else {
        
            open ( $vcf_in , "<" , $vcf_file ) or die "$!";
        
        }

        say "Reading from VCF file $vcf_file";
               
        while (<$vcf_in>) {
            
            chomp;
            
            if ( $_ =~ m/^\#/ ) {
            
                say $out $_;
            
            }
            
            elsif ( $_ =~ m/^(\S+)\s(\d+)\s(\S+)\s(\S)\s(\S)\s(.+)\s(GT\:\S+)\s(.+)/ ) {
            
				my ( $chrom , $pos , $id ) = ( $1 , $2 , $3 );
				
				my $minor;
            
				if ( defined $minor_table->{ "seq_s_" . $chrom }->{ $pos } ) {
				
					$minor = $minor_table->{ "seq_s_" . $chrom }->{ $pos };
				
				}
				
				elsif ( defined $minor_table->{ "seq_c_" . $chrom }->{ $pos } ) {
				
					$minor = $minor_table->{ "seq_c_" . $chrom }->{ $pos };
				
				}
				
				if ( defined $minor ) {
				
					my ( $ref , $alt ) = ( $4 , $5 );
					
					if ( $alt eq $minor ) {
					
						say $out $_;
					
					}
					
					elsif ( $ref eq $minor ) {
					
						my ( $info , $gt_format , $gt_string ) = ( $6 , $7 , $8 );
					
						my $new_gt_string;
					
						say "Recoding $chrom $pos ...";
					
						foreach my $gt ( split ( /\t/ , $gt_string ) ) {
					
							if ( $gt =~ m/^0\|0\:(.+)/ ) {
							
								$new_gt_string .= "\t" . "1|1\:" . $1;
							
							}
						
							elsif ( $gt =~ m/^0\|1\:(.+)/ ) {
							
								$new_gt_string .= "\t" . "1|0\:" . $1;
							
							}
							
							elsif ( $gt =~ m/^1\|0\:(.+)/ ) {
							
								$new_gt_string .= "\t" . "0|1\:" . $1;
							
							}
							
							elsif ( $gt =~ m/^1\|1\:(.+)/ ) {
							
								$new_gt_string .= "\t" . "0|0\:" . $1;
							
							}
						
						}
						
						say $out "$chrom\t$pos\t$id\t$alt\t$ref\t$info\t${gt_format}${new_gt_string}";
					
					}
				
				}
        
            }
            
        }
            
        close $vcf_in;
        
    }
    
    else {
        
        die "Unable to locate VCF formatted input file $vcf_file: $!";
        
    }
    
    close $out;
    
}

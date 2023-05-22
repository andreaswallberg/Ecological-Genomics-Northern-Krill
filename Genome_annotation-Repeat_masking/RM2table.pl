#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

use Getopt::Long;
use Data::Dumper;

# Option table and defaults

my $opt = {

	'psi' => [] ,
	'hmm' => [] ,
	
};

GetOptions(
    $opt ,

    # Input
    # -----

    # Allele file
    'tag=s' ,
    
    'rm=s' ,
    'rc=s' ,
    'te=s' ,
    'psi=s{,}' ,
    'hmm=s{,}' ,

);

my $tag = $opt->{ 'tag' };

my $fasta_file = $opt->{ 'rm' };
my $repeat_classifier_file = $opt->{ 'rc' };
my $te_file = $opt->{ 'te' };

my $ltr_table = {};

my @repeats;

my @origins = ( 'RepeatModeler' , 'RepeatClassifier' , 'TEsorter' , 'TransposonPSI' , 'HMMsearch' );

read_fasta( 'RepeatModeler' , $fasta_file );
read_fasta( 'RepeatClassifier' , $repeat_classifier_file );
read_fasta( 'TEsorter' , $te_file );

read_transposonpsi( 'TransposonPSI' );
read_hmmsearch( 'HMMsearch' );

print_table();

sub read_fasta {

	my ( $origin , $file ) = @_;
	
	open ( my $fasta_in , "<" , $file ) or die "$!";

	my $tmp_header;
	
	while (<$fasta_in>) {
	
		chomp;
		
		# >seq_c_23468_111710_115277_LTR#Unknown 
		
		if ( $_ =~ m/^\>/ ) {
		
			$tmp_header = undef;
		
		}
		
		if ( $_ =~ m/^\>(\S+)\#(\S+)/ ) {
		
			# say "$origin\t$1\t$2";
			
			unless ( defined $ltr_table->{ $1 } ) {
			
				push @repeats , $1;
			
			}
			
			$tmp_header = $1;
			
			$ltr_table->{ $1 }->{ $origin }->{ 'type' } = $2;
			$ltr_table->{ $1 }->{ $origin }->{ 'header' } = $_;
		
		}
		
		elsif ( not $_ =~ m/^\>/ and defined $tmp_header ) {
		
			$ltr_table->{ $tmp_header }->{ $origin }->{ 'seq' } .= $_;
		
		}
	
	}
	
}

sub read_transposonpsi {

	my ( $origin ) = @_;
	
	foreach my $transposonpsi_gff_file ( @{ $opt->{ 'psi' } } ) {
	
		open ( my $gff_in , "<" , $transposonpsi_gff_file ) or die "$!";

		# seq_c_100191_32326_40855_LTR#LTR/Gypsy	TransposonPSI	translated_nucleotide_match	1955	2701	130	+	.	ID=chain00001; Target=ltr_Roo; E=2e-10 233 451
		
		my $tmp_ltr_table = {};
		
		while (<$gff_in>) {
		
			chomp;
			
			if ( $_ =~ m/^(\S+)\#\S+\s+\S+\s+\S+\s+(\d+)\s(\d+)\s.+Target\=([^\;]+)/ ) {
			
				say "$origin\t$1\t$2\t$3\t$4";
				
# 				unless ( defined $ltr_table->{ $1 } ) {
# 				
# 					push @repeats , $1;
# 				
# 				}

				if ( defined $tmp_ltr_table->{ $1 } ) {
				
					$tmp_ltr_table->{ $1 }->{ $4 } += $3 - $2 + 1;
			
				}
			
			}
		
		}
		
		foreach my $repeat ( keys %{ $tmp_ltr_table } ) {
		
			my @sorted_types = sort {
				$tmp_ltr_table->{ $repeat }->{ $b } <=>
				$tmp_ltr_table->{ $repeat }->{ $a }
			} keys %{ $tmp_ltr_table->{ $repeat } };
		
			my $type = $sorted_types[ 0 ];
			
			$ltr_table->{ $repeat }->{ $origin }->{ 'type' } = $type;
		
		}
	
	}
	
}

sub read_hmmsearch {

	my ( $origin ) = @_;
	
	foreach my $hmm_file ( @{ $opt->{ 'hmm' } } ) {
	
		open ( my $hmm_in , "<" , $hmm_file ) or die "$!";
	
		my $tmp_ltr_table = {};
		
		# rnd-6_family-480_RM_UN#Unknown -          7SLRNA               DF0000016.4     256     312     474     414     494     403     640    -         3.2    5.7   3.7  ( Recon Family Size = 42, Final Multiple Alignment Size = 39 )
		
		while (<$hmm_in>) {
		
			chomp;
			
			if ( $_ =~ m/^(\w+\-\d+_\w+\-\d+).*#/ ) {

                my $repeat = $1;
			
                my @data = split ( /\s+/ , $_ );
                
                my ( $type , $e ) = @data[ 2 , 12 ];
                
                if ( $e <= 1e-3 and not $type =~ m/^DR/ ) {
                
                    unless ( defined $ltr_table->{ $repeat } ) {
                    
                        push @repeats , $repeat;
                    
                    }
                
                    if (
                    
                        defined
                        $tmp_ltr_table->{ $repeat }->{ $type } and $e <
                        $tmp_ltr_table->{ $repeat }->{ $type }
                        
                    ) {
                
                        $tmp_ltr_table->{ $repeat }->{ $type } = $e;
                
                    }
                
                    else {
                
                        $tmp_ltr_table->{ $repeat }->{ $type } = $e;
        
                    }
        
                }
			
			}
		
		}
		
		foreach my $repeat ( keys %{ $tmp_ltr_table } ) {
		
			my @sorted_types = sort {
				$tmp_ltr_table->{ $repeat }->{ $a } <=>
				$tmp_ltr_table->{ $repeat }->{ $b }
			} keys %{ $tmp_ltr_table->{ $repeat } };
		
			my $type = $sorted_types[ 0 ];
			
			# say $repeat;
			
			$ltr_table->{ $repeat }->{ $origin }->{ 'type' } = $type;
		
		}
	
	}
	
}

sub print_table {

	say "NR\tREPEAT\t" , join ( "\t" , @origins ) , "\tFINAL_TYPE";

	# say "NR\tREPEAT\t" , join ( "\t" , map { $_ . "\t" . $_ } @origins );

	my $i = 1;

	foreach my $repeat ( @repeats ) {
	
		my $string = "$i\t$repeat";

		foreach my $origin ( @origins ) {

			my $type = "";
			
			if ( defined $ltr_table->{ $repeat }->{ $origin }->{ 'type' } ) {
			
				$type = $ltr_table->{ $repeat }->{ $origin }->{ 'type' };
			
			}
			
			$string .= "\t$type";

		}
		
		my $new_repeat = $repeat . "_" . $tag;
		my $type = "";
		
		if (
		
			defined $ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'type' } and not
			$ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'type' } =~ m/unknown/i
			
		) {
		
			$new_repeat .= "_ReC";
			$type = $ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'type' };
		
		}
		
		elsif (
		
			defined $ltr_table->{ $repeat }->{ 'TEsorter' }->{ 'type' }
			and not
			$ltr_table->{ $repeat }->{ 'TEsorter' }->{ 'type' } =~ m/unknown/i
			
		) {
		
			$new_repeat .= "_TEs";
			$type = $ltr_table->{ $repeat }->{ 'TEsorter' }->{ 'type' };
		
		}
		
		elsif (
		
			defined $ltr_table->{ $repeat }->{ 'TransposonPSI' }->{ 'type' }
			and not
			$ltr_table->{ $repeat }->{ 'TransposonPSI' }->{ 'type' } =~ m/unknown/i
			
		) {
		
			$new_repeat .= "_PSI";
			$type = $ltr_table->{ $repeat }->{ 'TransposonPSI' }->{ 'type' };
		
		}
		
		elsif (
		
			defined $ltr_table->{ $repeat }->{ 'RepeatModeler' }->{ 'type' }
			and not
			$ltr_table->{ $repeat }->{ 'RepeatModeler' }->{ 'type' } =~ m/unknown/i
			
		) {
		
			$new_repeat .= "_RM";
			$type = $ltr_table->{ $repeat }->{ 'RepeatModeler' }->{ 'type' };
		
		}
		
		elsif (
		
			defined $ltr_table->{ $repeat }->{ 'RepeatModeler' }->{ 'type' }
			and not
			$ltr_table->{ $repeat }->{ 'RepeatModeler' }->{ 'type' } =~ m/unknown/i
			
		) {
		
			$new_repeat .= "_RM";
			$type = $ltr_table->{ $repeat }->{ 'RepeatModeler' }->{ 'type' };
		
		}

		elsif (
		
			defined $ltr_table->{ $repeat }->{ 'HMMsearch' }->{ 'type' }
			
		) {
		
			$new_repeat .= "_DFA";
			
			$type = "SINE/" . $ltr_table->{ $repeat }->{ 'HMMsearch' }->{ 'type' };
		
            # say "$repeat\t$type";
		
		}
		
		my $repeat_group = {

			# TEsorter
			
			"LTR/Bel-Pao" => "LTR/Pao" ,
		
			# TransposonPSI
		
			"LINE" => "LINE" ,
			"ltr_Roo" => "LTR/Pao" ,
			"TY1_Copia" => "LTR/Copia" ,
			"gypsy" => "LTR/Gypsy" ,
			"hAT" => "DNA/hAT" ,
			"DDE_1" => "DNA/DDE" ,
			"piggybac" => "DNA/PiggyBac" ,
			"mariner" => "DNA/TcMar-Mariner" ,
			"MuDR_A_B" => "DNA/MULE-MuDR" ,
			"helitronORF" => "RC/Helitron" ,
			"ISC1316" => "DNA/ISC1316" ,
			"cacta" => "DNA/CMC-EnSpm" ,
			"P_element" => "DNA/P_element" ,
			"Crypton" => "DNA/Crypton" ,
			"mariner_ant1" => "DNA/TcMar-Ant1" ,

		};
		
		unless ( $type ) {
		
			$new_repeat .= "_UN";
			$type = "Unknown";
		
		}
		
		if ( defined $repeat_group->{ $type } ) {
		
			$type = $repeat_group->{ $type };
		}
		
		$string .= "\t$type";
		
		say $string;
	
		my $header = $ltr_table->{ $repeat }->{ 'RepeatModeler' }->{ 'header' };
	
		$header =~ s/^\S+\#\S+//;
	
		say STDERR ">$new_repeat#$type $header\n" , $ltr_table->{ $repeat }->{ 'RepeatModeler' }->{ 'seq' };
	
		$i++;

	}

}

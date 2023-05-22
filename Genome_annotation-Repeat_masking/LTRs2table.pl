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
    
    'lr=s' ,
    'rc=s' ,
    'te=s' ,
    'psi=s{,}' ,
    'hmm=s{,}' ,

);

my $tag = $opt->{ 'tag' };

my $fasta_file = $opt->{ 'lr' };
my $repeat_classifier_file = $opt->{ 'rc' };
my $te_file = $opt->{ 'te' };

my $ltr_table = {};

my @repeats;

my @origins = ( 'LTR_Retriever' , 'RepeatClassifier' , 'TEsorter' , 'TransposonPSI' , 'HMMsearch' );

read_fasta( 'LTR_Retriever' , $fasta_file );
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
		
		if ( $_ =~ m/^\>(\S+\_\w+\_\d+\_\d+\_\d+)\_\w+\#(\S+)/ ) {
		
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
			
			if ( $_ =~ m/^(\S+\_\w+\_\d+\_\d+\_\d+)\_\S+\s+\S+\s+\S+\s+(\d+)\s(\d+)\s.+Target\=([^\;]+)/ ) {
			
				# say "$origin\t$1\t$2\t$3\t$4";
				
				unless ( defined $ltr_table->{ $1 } ) {
				
					push @repeats , $1;
				
				}
				
				
				$tmp_ltr_table->{ $1 }->{ $4 } += $3 - $2 + 1;
			
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
		
		# seq_s_73_823510_826043_LTR#Unknown         -           2534 Gypsy-31-LTR_DR      DF0002487.1   954   7.2e-11   34.9  27.5   1   2   1.2e-08   5.3e-07   22.1  13.4    41   331   152   472   120   483 0.38 -
		
		while (<$hmm_in>) {
		
			chomp;
			
			if ( $_ =~ m/^(\S+)\_LTR\#/ ) {

                my $repeat = $1;
			
                my @data = split ( /\s+/ , $_ );
                
                my ( $type , $start , $stop ) = @data[ 3 , 15 , 16 ];
			
				if ( not $type =~ m/^RLT/i ) {
				
                    unless ( defined $ltr_table->{ $repeat } ) {
                    
                        push @repeats , $repeat;
                    
                    }
                    
                    my $length = $stop - $start + 1;
                    
                    if ( $length >= 100 ) {
                    
                        $tmp_ltr_table->{ $repeat }->{ $type } += $length;
			
                    }
			
                }
			
			}
		
		}
		
		foreach my $repeat ( keys %{ $tmp_ltr_table } ) {
		
			my @sorted_types = sort {
				$tmp_ltr_table->{ $repeat }->{ $b } <=>
				$tmp_ltr_table->{ $repeat }->{ $a }
			} keys %{ $tmp_ltr_table->{ $repeat } };
		
			my $type = $sorted_types[ 0 ];
			
			say $repeat;
			
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
		
			defined $ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'type' }
			and
			$ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'type' } =~ m/LTR/i
			
		) {
		
			$new_repeat .= "_ReC";
			$type = $ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'type' };
		
		}
		
		elsif (
		
			defined $ltr_table->{ $repeat }->{ 'TEsorter' }->{ 'type' }
			and
			$ltr_table->{ $repeat }->{ 'TEsorter' }->{ 'type' } =~ m/LTR/i
			
		) {
		
			$new_repeat .= "_TEs";
			$type = $ltr_table->{ $repeat }->{ 'TEsorter' }->{ 'type' };
		
		}
		
		elsif (
		
			defined $ltr_table->{ $repeat }->{ 'TransposonPSI' }->{ 'type' }
			and
			$ltr_table->{ $repeat }->{ 'TransposonPSI' }->{ 'type' } =~ m/^(LTR|TY1_Copia|gypsy)/i
			
		) {
		
			$new_repeat .= "_PSI";
			$type = $ltr_table->{ $repeat }->{ 'TransposonPSI' }->{ 'type' };
		
		}
		
		elsif (
		
			defined $ltr_table->{ $repeat }->{ 'HMMsearch' }->{ 'type' }
			
		) {
		
			$new_repeat .= "_DFA";
			
			$type = $ltr_table->{ $repeat }->{ 'HMMsearch' }->{ 'type' };
		
            if ( $type =~ m/Copia/i ) {
            
                $type = "LTR/Copia";
            
            }
            
            elsif ( $type =~ m/Bel/i ) {
            
                $type = "LTR/Pao";
            
            }
            
            elsif ( $type =~ m/Gypsy/i ) {
            
                $type = "LTR/Gypsy";
            
            }
            
            elsif ( $type =~ m/DIRS/i ) {
            
                $type = "LTR/DIRS";
            
            }
		
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
		
		if ( defined $repeat_group->{ $type } ) {
		
			$type = $repeat_group->{ $type };
		}
		
		unless ( $type ) {
		
            if (
            
                defined $ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'type' }
                
            ) {
            
                $new_repeat .= "_ReC";
                $type = $ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'type' };
            
            }
		
		}
		
		$string .= "\t$type";
		
		if ( $type =~ m/LTR/i or $type =~ m/Unknown/i ) {
		
			say $string;
		
			say STDERR ">$new_repeat#$type\n" , $ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'seq' };
		
			$i++;
		
		}

	}

}

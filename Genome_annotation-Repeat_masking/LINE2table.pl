#!/usr/bin/perl -w

# Copyright (c) 2023 Andreas Wallberg (Uppsala University)

use strict;
use warnings;
use 5.010;

my $tag = shift @ARGV;
my $fasta_file = shift @ARGV;
my $repeat_classifier_file = shift @ARGV;
my $te_file = shift @ARGV;

my $repeat_group = {

	# TEsorter
	
	"LTR/Bel-Pao" => "LTR/Pao" ,

	# TransposonPSI

	"LTR/Roo" => "LTR/Pao" ,
	
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

my $ltr_table = {};

my @origins = ( 'TransposonPSI' , 'RepeatClassifier' , 'TEsorter' );

my @repeats;

read_fasta( 'TransposonPSI' , $fasta_file );
read_fasta( 'RepeatClassifier' , $repeat_classifier_file );
read_fasta( 'TEsorter' , $te_file );

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
		
			my $type = $2;
			
			if ( defined $repeat_group->{ $type } ) {
			
				$type = $repeat_group->{ $type };
			}
			
			unless ( defined $ltr_table->{ $1 } ) {
			
				push @repeats , $1;
			
			}
			
			$tmp_header = $1;
			
			$ltr_table->{ $1 }->{ $origin }->{ 'type' } = $type;
			$ltr_table->{ $1 }->{ $origin }->{ 'header' } = $_;
		
		}
		
		elsif ( not $_ =~ m/^\>/ and defined $tmp_header ) {
		
			$ltr_table->{ $tmp_header }->{ $origin }->{ 'seq' } .= $_;
		
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
		
		# my @origins = ( 'LTR_Retriever' , 'RepeatClassifier' , 'TEsorter' , 'TransposonPSI' );
		
		my $new_repeat = $repeat . "_" . $tag;
		my $type = "";
		
		if (
		
			defined $ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'type' }
			
		) {
		
			unless (
			
				$ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'type' } =~ m/^LINE/ or
				$ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'type' } =~ m/^Unknown/
				
			) {
			
				next;
			
			}
		
		}
		
		if (
		
			defined $ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'type' }
			and
			$ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'type' } =~ m/^LINE/i
			
		) {
		
			$new_repeat .= "_ReC";
			$type = $ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'type' };
		
		}
		
		elsif (
		
			defined $ltr_table->{ $repeat }->{ 'TEsorter' }->{ 'type' }
			and
			$ltr_table->{ $repeat }->{ 'TEsorter' }->{ 'type' } =~ m/^LINE/i
			
		) {
		
			$new_repeat .= "_TEs";
			$type = $ltr_table->{ $repeat }->{ 'TEsorter' }->{ 'type' };
		
		}
		
		elsif (
		
			defined $ltr_table->{ $repeat }->{ 'TransposonPSI' }->{ 'type' }
			and
			$ltr_table->{ $repeat }->{ 'TransposonPSI' }->{ 'type' } =~ m/^LINE/i
			
		) {
		
			$new_repeat .= "_PSI";
			$type = $ltr_table->{ $repeat }->{ 'TransposonPSI' }->{ 'type' };
		
		}
		
		$string .= "\t$type";
		
		if ( $type =~ m/LINE/i ) {
		
			say $string;
		
			say STDERR ">$new_repeat#$type\n" , $ltr_table->{ $repeat }->{ 'RepeatClassifier' }->{ 'seq' };
		
			$i++;
		
		}

	}

}

#!/usr/bin/perl

# ==============================================================================

# The script provided below is a modification of the purge.pl script originally
# written by Michael Roach (Australian Wine Research Institute) and part of the
# open source Purge Haplotigs tool v1.1.0. That tool can be found here:
# https://bitbucket.org/mroachawri/purge_haplotigs/wiki/Updates

# The script was modified by Andreas Wallberg (Uppsala University)

# It retains the original copyright and license notice below.

# ==============================================================================

# Copyright (c) 2017 Michael Roach (Australian Wine Research Institute)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


use strict;
use warnings;
use Getopt::Long;
use threads;

# NOTE: load shared
# --------------------

use threads::shared;
use Thread::Semaphore;
use Thread::Queue;
use FindBin qw($RealBin);
use List::Util qw/sum/;
use lib "$RealBin/../lib";
use PipeUtils;

# NOTE: for debugging
# ----------------------

use Data::Dumper;

#---INPUT PARAMETERS---

my $genome;
my $coverage_stats;
my $bam_file;
my $threads = 30; # NOTE:  = increase number of threads
my $limit_io;
my $maxmatch_cutoff = 250;
my $bestmatch_cutoff = 70;
my $low_cutoff = 5;
my $hit_cutoff = 1;
my $out_prefix = 'curated';
my $dotplots;
my $falconNaming;
my $repeats;
my $min_win_size = 5000;
my $max_windows = 200;
# my $minimap2_parameters = '-p 1e-5 -f 0.001 -N 100000';

# NOTE: change minimap2 parameters to fit large genome
# -------------------------------------------------------

my $minimap2_parameters = "-p 1e-5 -f 0.001 -N 10000";

my $minimiser_drop = '4G';
my $xvfb=" ";
our $verbose;

#---HELP MESSAGE---

my $usage = "
USAGE:
purge_haplotigs  purge  -g genome.fasta  -c coverage_stats.csv

REQUIRED:
-g / -genome        Genome assembly in fasta format. Needs to be indexed with samtools faidx.
-c / -coverage      Contig by contig coverage stats csv file from the previous step.

OPTIONAL:
-t / -threads       Number of worker threads to use. DEFAULT = $threads
-o / -outprefix     Prefix for the curated assembly. DEFAULT = \"$out_prefix\"
-r / -repeats       BED-format file of repeats to ignore during analysis.
-d / -dotplots      Generate dotplots for manual inspection.
-b / -bam           Samtools-indexed bam file of aligned and sorted reads/subreads to the
                    reference, required for generating dotplots.
-f / -falconNaming  Rename contigs in the style used by FALCON/FALCON-unzip

ADVANCED:
-a / -align_cov     Percent cutoff for identifying a contig as a haplotig. DEFAULT = $bestmatch_cutoff
-m / -max_match     Percent cutoff for identifying repetitive contigs. Ignored when 
                    using repeat annotations (-repeats). DEFAULT = $maxmatch_cutoff
-I                  Minimap2 indexing, drop minimisers every N bases, DEFAULT = $minimiser_drop
-v / -verbose       Print EVERYTHING.
-limit_io           Limit for I/O intensive jobs. DEFAULT = -threads
-wind_min           Min window size for BED coverage plots (for dotplots). DEFAULT = $min_win_size
-wind_nmax          Max windows per contig for BED coverage plots (for dotplots). DEFAULT = $max_windows

";


#---CHECK PROGRAMS---

if (!(check_programs qw/bedtools minimap2 samtools/)){
    err('ONE OR MORE REQUIRED PROGRAMS IS MISSING');
}


#---PARSE ARGUMENTS---

my $args = "@ARGV";
GetOptions (
    "genome=s" => \$genome,
    "coverage_stats=s" => \$coverage_stats,
    "bam=s" => \$bam_file,
    "threads=i" => \$threads,
    "limit_io=i" => \$limit_io,
    "outprefix=s" => \$out_prefix,
    "align_cov=i" => \$bestmatch_cutoff,
    "max_match=i" => \$maxmatch_cutoff,
    "I=s" => \$minimiser_drop,
    "dotplots" => \$dotplots,
    "falconNaming" => \$falconNaming,
    "repeats=s" => \$repeats,
    "wind_min=i" => \$min_win_size,
    "wind_nmax=i" => \$max_windows,
    "verbose" => \$verbose
) or die $usage;

# required parameters
(($genome) && ($coverage_stats)) or die $usage;
($limit_io) or $limit_io = $threads;


# require bam file if generating dotplots and check X11 or Xvfb
if ($dotplots){
    ($bam_file) && (-s $bam_file) or err("BAM file required when producing dotplots\n$usage");
    if (!(check_programs qw/Rscript/)){
        err('ONE OR MORE REQUIRED PROGRAMS IS MISSING');
    }
    if (system("Rscript -e 'png()' 2> /dev/null")!=0){
        if (check_programs("xvfb-run")){
            $xvfb = "xvfb-run -a";
        } else {
            msg('

WARNING:
Neither X11 nor xvfb is available; unable to produce dotplots.
');
            undef $dotplots;
        }
    }
}


# check the files
if (!(check_files($genome, $coverage_stats))){
    msg('One or more files missing, exiting');
    die $usage;
}

# check optional files
if (($repeats) && !(check_files($repeats))){
    msg('One or more files missing, exiting');
    die $usage;
}
if (($dotplots) && !(check_files($bam_file))){
    msg('One or more files missing, exiting');
    die $usage;
}


#---SET UP LOGGING---
our $LOG;
my $TMP_DIR = 'tmp_purge_haplotigs';
(-d $TMP_DIR) or mkdir $TMP_DIR;
open $LOG, '>', "$TMP_DIR/purge_haplotigs_purge.log" or die "failed to open log file for writing";



#---GLOBAL VARIABLES---

# directories
my $MINCE_DIR = "$TMP_DIR/CONTIGS";
my $TMP_ALN = "$TMP_DIR/TMP_ALN";
my $COV_DIR = "$TMP_DIR/COV";
my $ERR_DIR = "$TMP_DIR/STDERR";
my $ASSIGNED = 'dotplots_reassigned_contigs';
my $UNASSIGNED = 'dotplots_unassigned_contigs';
my $PTH;

# files
my $assembly_fasta = "$TMP_DIR/assembly.fasta";
my $assembly_index = "$TMP_DIR/assembly.$minimiser_drop.mmi";
my $suspects_fasta = "$TMP_DIR/suspects.fasta";
my $minimap2PAF = "$TMP_DIR/minimap2.$minimiser_drop.paf";
my $minimap2IDX = "$minimap2PAF.index";
my $hit_summary = "$TMP_DIR/hit_summary.tsv";
my $assembly_repeats = "$TMP_DIR/assembly.repeats.bed";
my $assembly_cov = "$TMP_DIR/assembly.coverage.bed";
my $assembly_logcov = "$TMP_DIR/assembly.logcov.bed";

# NOTE: files
# -----------

my $minimap2PAF_simple = "$minimap2PAF.simplified.paf";
my $minimap2IDX_simple = "$minimap2PAF_simple.index";

# file indexes
my %mmIndex;                    # $mmIndex{contig} = file read start for minimap2 alignments
my %faIndex;                    #                    file read start for genome fasta file
my %repIndex;                   #                    file read start for repeats bed file
my %covIndex;                   #                    file read start for coverage bed file

# some clean file names
my $genome_file_name = $genome;
$genome_file_name =~ s/.*\///;

# reassignment step
my %suspects;   # suspects flagged from coverage analysis
my %junk;       # junk flagged from coverage analysis

my %contigLEN :shared;          # length of contig
my %contigMLEN :shared;         # length minus repetitive regions
my %contigHIT1 :shared;         # first best remaining hit
my %contigHIT2 :shared;         # second best remaining hit
my %contigBM :shared;           # best match score (percent contig aligned to either hit)
my %contigMM :shared;           # max match score (sum of all alignments as perc of contig)
my %contigASSIGN :shared;       # contig reassignment, r=repeat, h=haplotig, n=no reassign, u=unknown
my %contigREASSIGN :shared;     # 1/0, flag if it's been reassinged
my %contigRENAME;               # for renaming haplotigs
my %contigOVPURG :shared;       # flag, contig is over-purged
my %hits :shared;               # hash of annonymous array refs of contig hits
my $over_purge_mode;            # for over-purge checking, mode will either be 'detect' or 'fix'

# NOTE: Add some variables
# ------------------------

my %suspects_simple;            # An updated and simplified table of suspects will go here
my %length_table_after_repeats; # A table with contig lengths after correcting for repeats
my %repeat_table;               # A table with repeat coordinates
my $n_contigs_done :shared;     # A counter

# ------------------------

#---OUTPUT FILES---

my $contig_paths = "$out_prefix.contig_associations.log";
my $out_fasta = "$out_prefix.fasta";
my $out_haplotigs = "$out_prefix.haplotigs.fasta";
my $out_artefacts = "$out_prefix.artefacts.fasta";
my $out_reassignments = "$out_prefix.reassignments.tsv";



#---THREADS---

my $available_threads = Thread::Semaphore->new($threads);
my $max_jobs = Thread::Semaphore->new($limit_io);
my $writing_to_out = Thread::Semaphore->new(1);
my $queueBedtools = Thread::Queue->new();
my $queueHitSummary = Thread::Queue->new();
my $queuePairwise;
my $queueOverPurge;

#---OPEN DIRS ETC---

pre_flight();



#---PRINT PARAMETERS---
my $param_message = "
Beginning Pipeline

PARAMETERS:
Genome fasta:           $genome
Coverage csv:           $coverage_stats";
if ($bam_file){
    $param_message .= "
Produce dotplots:       TRUE
Bam file:               $bam_file
Min cov window len:     $min_win_size bp
max ctg cov windows:    $max_windows";
} else {
    $param_message .= "
Produce dotplots:       FALSE";
}
$param_message .= ($falconNaming) ?
"\nFalcon-style naming:    TRUE" :
"\nFalcon-style naming:    FALSE";
if ($repeats){
    $param_message .= "
Repeat annotations:     $repeats";
}
$param_message .= "
Threads:                $threads
I/O intense jobs:       $limit_io
Cutoff, alignment:      $bestmatch_cutoff %
Cutoff, repeat:         $maxmatch_cutoff %
Cutoff, suspect:        $low_cutoff %
Out prefix:             $out_prefix
minimap2 parameters:    '$minimap2_parameters'

Running using command:
purge_haplotigs purge $args\n\n";

msg($param_message);


#---PIPELINE BEGIN---

msg("\n\nPREPARATION\n");

# read in fasta.fai
read_fasta_fai();

# read in coverage stats
read_cov_stats();

# mince genome, this will make later steps run much faster
index_genome();

# make bed windows, reads per window, for each contig. Used later for juxtaposing read-depth coverage with dotplots
get_window_coverage() if ($dotplots);

# make minimap2 index, run alignments, sumamrise hits
run_minimap2_alignments();

# NOTE: parse the repeats before we generate the hit summary
# ----------------------------------------------------------

# parse repeats if provided
if ($repeats){

    parse_repeats();
    
    # NOTE: simplify the minimap2 PAF and generate the hit summary
    # ------------------------------------------------------------

    if (-s $minimap2PAF){
        msg('Simplifying the minimap2 PAF');
        
        simplify_paf();
        
    }
    
    # Set simplified PAF as main PAF

	$minimap2PAF = $minimap2PAF_simple;
	$minimap2IDX = $minimap2IDX_simple;

	# Re-read the index
	
	%mmIndex = ();
	
	# -------------------------------
	
}

# read alignment index
read_minimap2_index();

# NOTE: No need to re-generate the hit summary here
# -------------------------------------------------

# generate the hit summary
# hit_summary();



#---ITERATIVE STEP---

my $convergence = 0;

while(!($convergence)){
	
	# NOTE: print a message here
	# --------------------------
	
	msg("\n\nChecking for convergence ...\n");
	
    # read through blastn hit summary and get top 2 matches for each suspect contig
    get_contig_hits();

    # run mummer steps
    pairwise_alignments();
    
    # remove conflict reassignments
    check_assignments();
    
    # add to reassignments list
    add_reassignments();
}


#---GENERATE OUTPUT---


# check overpurging
over_purge_check();

msg("\n\nGENERATING OUTPUT\n");

# get the reassignment paths
write_contig_associations();

# write the table and new assembly
write_assembly();


msg("\n\nPURGE HAPLOTIGS HAS COMPLETED SUCCESSFULLY!\n");

exit(0);


#---SUBROUTINES---

sub pre_flight {
    
    # directories
    for my $dir ($TMP_DIR, $MINCE_DIR, $TMP_ALN, $COV_DIR, $ERR_DIR){
        if (!(-d $dir)){
            mkdir $dir or err("failed to create directory $dir");
        }
    }

    if ($dotplots){
        for my $dir ($ASSIGNED, $UNASSIGNED){
            if (!(-d $dir)){
                mkdir $dir or err("failed to create directory $dir");
            }
        }
    }
    
    # cleanup
    for my $file ($out_artefacts, $out_fasta, $out_haplotigs, $out_reassignments, $contig_paths){
        if (-s $file){
            unlink $file or err("failed to clean up previous run output file: $file");
        }
    }
    
    return;
}



sub read_fasta_fai {
    # index if needed
    if (!(-s "$genome.fai")){
        msg("Indexing $genome");
        runcmd({ command => "samtools faidx $genome 2> $ERR_DIR/samtools.faidx.stderr",
                 logfile => "$ERR_DIR/samtools.faidx.stderr",
                 silent => 1 });
    }
    
    # read contig lengths
    msg("Reading $genome.fai");
    
    open my $FAI, '<',  "$genome.fai" or err ("failed to open  $genome.fai for reading");

	# NOTE: read FAI index faster
	# ---------------------------

    while(my $l = <$FAI>){
    
		chomp $l;
		
		if ( $l =~ m/^(\S+)\s(\d+)/ ) {
            $contigLEN{$1} = $2;
            $contigMLEN{$1} = $2;
		}
    }
    
    # ---------------------------

    # ORIG:
    
#     while(my $l = <$FAI>){
#         my @line = split /\s+/, $l;
#         if (!defined($line[0]) || !defined($line[1])){
#             err("bad entry in $genome.fai index file? line:\n$l");
#         } else {
#             $contigLEN{$line[0]} = $line[1];
#             $contigMLEN{$line[0]} = $line[1];
#         }
#     }

    close $FAI;

    return;
}



sub read_cov_stats {
    msg("Reading $coverage_stats");
    
    open my $CSV, '<', $coverage_stats or err("failed to open $coverage_stats for reading");
    
    while(my $l = <$CSV>){
        next if ($l =~ /^#/);
        
        my @line = split /,/, $l;
        if (($line[0]) && !($contigLEN{$line[0]})){
            err("no contig \"$line[0]\" in genome index file, csv line:\n$l");
        } elsif ($line[1] eq 's'){
            $suspects{$line[0]} = 1;
        } elsif ($line[1] eq 'j'){
            $junk{$line[0]} = 1;
        }
    }
    
    if (!(keys %suspects)){
        if (!(keys %junk)){
            err('No contigs flagged as either suspects or artefacts, nothing to do, exiting...');
        }
        msg('WARNING: no contigs flagged as suspects, dropping flagged artefacts and exiting...');
        write_assembly();
        exit(0);
    }
    
    close $CSV;
    
    return;
}



sub index_genome {
    
    msg("Scanning $genome");
    
    open my $GEN, '<', $genome or err("failed to open $genome for reading");
    my $offset = tell($GEN);
    while(<$GEN>){
        if($_=~/>(\S+)/){
            $faIndex{$1}=$offset;
        }
        $offset = tell($GEN);
    }
    close $GEN;
    
    return;
}



sub get_seq {
    my $ctg = $_[0];
    my $seq;
    
    open my $GEN, '<', $genome or err("failed to open $genome for reading");
    
    # check the seq was indexed
    defined($faIndex{$ctg}) or err("$ctg not in $genome");
    
    # fast forward
    seek($GEN, $faIndex{$ctg}, 0);
    
    # check the header
    $seq = <$GEN>;
    if($seq=~/>(\S+)/){
        ($ctg eq $1) or err("wrong seq returned from $genome, expected $ctg, returned $1");
    } else {
        err("error reading $ctg from $genome");
    }
    
    # slurp the seq
    while(<$GEN>){
        if($_!~/^>/){
            $seq .= $_;
        } else {
            last;
        }
    }
    
    return $seq;
}



sub get_window_coverage {
    # check the bam is indexed
    if (!(-s "$bam_file.bai")){
        msg("Indexing $bam_file");
        runcmd({ command => "samtools index $bam_file 2> $ERR_DIR/samtools.index.stderr",
                 logfile => "$ERR_DIR/samtools.index.stderr",
                 silent => 1 });
    }

    # check for temp files from previous runs
    if (-s "$assembly_cov.tmp"){
        unlink "$assembly_cov.tmp" or err("Failed to clean up temp file $assembly_cov.tmp");
    }
    
    # make the bed windows for each contig
    if (!(-s $assembly_cov)){
        msg('Getting windowed read-depth for each contig');
        bedtools_multicov();
    } else {
        msg('Reusing windowed read-depth from previous run');
    }
    
    # convert to log2 read-depth
    if (!(-s $assembly_logcov)){
        msg('Generating log2(read-depth / average read-depth) coverages for plotting later');
        
        cov_2_logcov(average_read_depth());
    }
    
    index_coverage();

    return;
}



sub bedtools_multicov {
    for my $ctg ( sort { $contigLEN{$b} <=> $contigLEN{$a} } keys %suspects){
        # add the job to the queue
        $queueBedtools->enqueue($ctg);
    }
    
    # finalise the queue
    $queueBedtools->end();
    
    # spawn the worker threads
    for (1..$limit_io){
        $available_threads->down(1);
        threads->create(\&multicov_job);
    }
    
    # wait for workers to finish
    $available_threads->down($threads);
    $available_threads->up($threads);
    
    # join workers
    for my $thr (threads->list()){
        $thr->join();
    }
    
    # finish multicov file
    rename "$assembly_cov.tmp", $assembly_cov;
    
    return;
}



sub multicov_job {
    # iterate over the queue until all jobs are done
    while (defined(my $ctg = $queueBedtools->dequeue())) {
        
        # calculate window len
        my $win = int($contigLEN{$ctg} / $max_windows);
        $win = $min_win_size if ($win < $min_win_size);
        my $step = int($win / 2);

        # print genome file
        open my $GFAI, '>', "$COV_DIR/$ctg.g" or err("failed to open $COV_DIR/$ctg.g for writing");
        print $GFAI "$ctg\t$contigLEN{$ctg}\n";
        close $GFAI;

        # make windows
        runcmd({ command => "bedtools makewindows -g '$COV_DIR/$ctg.g' -w $win -s $step > '$COV_DIR/$ctg.wind' 2> $ERR_DIR/bedtools.mkwind.stderr",
                 logfile => "$ERR_DIR/bedtools.mkwind.stderr",
                 silent => 1 });
        
        my @mcov;
        
        # use a pipe to run bedtools multicov
        my $pipeCmd = "bedtools multicov -bams $bam_file -bed '$COV_DIR/$ctg.wind'";
        open my $INTMP, '-|', $pipeCmd or err("failed to open pipe: $pipeCmd");
        while(<$INTMP>){
            push @mcov, $_;
        }
        close $INTMP or err("failed to close Bedtools pipe: $pipeCmd | (this script)");
        
        # append to the assembly coverage
        $writing_to_out->down(1);
        open my $TMP, '>>', "$assembly_cov.tmp" or err("failed to open $assembly_cov.tmp for appending");
        print $TMP @mcov;
        close $TMP;
        $writing_to_out->up(1);
        
        # cleanup
        for my $file ("$COV_DIR/$ctg.g", "$COV_DIR/$ctg.wind"){
            unlink $file if (-e $file);
        }
    }
    
    # exit thread
    $available_threads->up(1);

    return;
}



sub average_read_depth {
    my @average_depth;

    # generate histogram table files for each contig, dump them in the minced directory
    open my $MCV, '<', "$assembly_cov" or err("Failed to open $assembly_cov for reading");
    
    while(<$MCV>){
        chomp;
        if ($_){
            my @l = split /\s+/;
            my $d = $l[3] / ($l[2] - $l[1]);
            push @average_depth, $d;
        }
    }
    close $MCV;

    # calculate the average
    my $avg = sum(@average_depth) / @average_depth;
    
    return $avg;
}



sub cov_2_logcov {
    my $average = $_[0];
    
    # open cov file for reading and tmp file for writing
    open my $COV, '<', "$assembly_cov" or err("Failed to open $assembly_cov for reading");
    open my $TMP, '>', "$assembly_logcov.tmp" or err("Failed to open $assembly_logcov.tmp for writing");

    # convert to log2 read-depth
    while(<$COV>){
        my@l = split /\s+/;
        $l[3] = $l[3] / ($l[2] - $l[1]);
        if ($l[3]==0){
            $l[3] = -1.5;
        } else {
            $l[3] = log($l[3] / $average) / log(2);
        }
        
        # we'll cap the extremes to -1.5 and 1.5
        $l[3] = -1.5 if ($l[3] < -1.5);
        $l[3] = 1.5 if ($l[3] > 1.5);

        print $TMP "$l[0]\t$l[1]\t$l[2]\t$l[3]\n";
    }
    
    close $COV;
    close $TMP;

    rename "$assembly_logcov.tmp", $assembly_logcov;
    
    return;
}


sub index_coverage {
    
    msg("Scanning $assembly_logcov");
    
    open my $LCV, '<', $assembly_logcov or err("failed to open $assembly_logcov for reading");
    
    my $curCtg = 'init';
    my $offset = tell($LCV);
    
    # generate the index file
    while(<$LCV>){
        my@l=split/\s+/;
        if ($curCtg ne $l[0]){
            $curCtg = $l[0];
            $covIndex{$l[0]} = $offset;
        }
        $offset = tell($LCV);
    }
    close $LCV;
    
    return;
}


sub get_cov {
    my $ctg = $_[0];
    my $out;
    
    open my $COV, '<', $assembly_logcov or err("Failed to open $assembly_logcov for reading");
    
    seek($COV, $covIndex{$ctg}, 0);
    
    while(<$COV>){
        my@l=split/\s+/;
        if ($l[0] eq $ctg){
            $out .= $_;
        } else {
            last;
        }
    }
    
    return $out;
}


sub run_minimap2_alignments {

	# NOTE: skip checking for the PAF index
	# -------------------------------------

#     # build index if needed
#     if (!(-s $assembly_index)){
#         msg('Building assembly index for minimap2');
#         
#         # open minimap2 pipe job
#         my $errLog = "$ERR_DIR/minimap2.idx.stderr";
#         my $pipeCmd = "minimap2 -I $minimiser_drop -t $threads -d $assembly_index.tmp - 2> $errLog";
#         open my $MMI, '|-', $pipeCmd or err("Failed to open minimap2 pipe command $pipeCmd");
#         
#         # sig handler
#         local $SIG{PIPE} = sub { 
#             err("Minimap2 pipe has died:\n(this script) | $pipeCmd\n\nCheck $errLog for possible cause.") 
#         };
#         
#         # pass seqs to pipe
#         for my $contig (sort keys %contigLEN){
#             if (!($junk{$contig})){
#                 print $MMI get_seq($contig);
#             }
#         }
#         
#         # done
#         rename "$assembly_index.tmp", $assembly_index;
#         close $MMI or err("Failed to close Minimap2 pipe: $pipeCmd");
#     }
#     

	# ----------------------------------

    # run the alingments if needed
    if ((-s $minimap2PAF)){
        msg('Reusing minimap2 alignments from previous run');
    } else {
        
        # start the minimap2 alignments
        run_minimap2();
        
        msg('Finished minimap2 alignments');
    }
    
    return;
}


sub run_minimap2 {

    msg('Performing minimap2 alignments');

    # cleanup old files
    for my $file("$minimap2PAF.tmp", "$minimap2IDX.tmp", $minimap2IDX){
        if (-s $file){
            unlink $file;
        }
    }
    
    my $errLog = "$ERR_DIR/minimap2.stderr";
    
    # NOTE: generate tmp suspects seq and "not_junk" sequence files
    # -------------------------------------------------------------
    
    open ( my $tmp_s_out , ">" , "minimap2.tmp.suspects.fa" ) or die "$!";
    open ( my $tmp_j_out , ">" , "minimap2.tmp.not_junk.fa" ) or die "$!";

    for my $contig (sort keys %contigLEN){
        if ($suspects{$contig}){
            print $tmp_s_out get_seq($contig);
        }
		if (!($junk{$contig})){
			print $tmp_j_out get_seq($contig);
		}
    }
    close $tmp_s_out;
    close $tmp_j_out;
    
	print "minimap2 -I $minimiser_drop $minimap2_parameters -t $threads -o $minimap2PAF.tmp minimap2.tmp.not_junk.fa minimap2.tmp.suspects.fa 2> $errLog\n";
	system ( "minimap2 -I $minimiser_drop $minimap2_parameters -t $threads -o $minimap2PAF.tmp minimap2.tmp.not_junk.fa minimap2.tmp.suspects.fa 2> $errLog && touch $minimap2PAF.tmp.done" );

    print "Finished mapping with minimap2!\n";
    
    # finalize the alignment output
    rename "$minimap2PAF.tmp", $minimap2PAF;

	# NOTE: No need to index here as we will make a simplified PAF
	# ------------------------------------------------------------

	# ORIG:
	
#     ## build the index of the alignments ###
#     msg('Indexing minimap2 alignments');
#     
#     # open the alignments for reading
#     open my $FWR, '<', $minimap2PAF or err("failed to open $minimap2PAF for reading");
#     
#     # open filehandle for alignment index
#     open my $IDX, '>', "$minimap2IDX.tmp" or err("failed to open $minimap2IDX.tmp for writing");
#     
#     # tracking for index file
#     my $curCtg = 'init';
#     my $offset = tell($FWR);
#     
#     # generate the index file
#     while(<$FWR>){
#         my ($ctg) = $_ =~ m/^(.+?)\s/;
#         if ($curCtg ne $ctg){
#             print $IDX $ctg, "\t", $offset, "\n";
#             $curCtg = $ctg;
#         }
#         $offset = tell($FWR);
#     }
#     close $FWR;
#     close $IDX;
#     
#     # finalize the index file
#     rename "$minimap2IDX.tmp", $minimap2IDX;
# 
#     
    return;
}


sub run_minimap2_orig {
    msg('Performing minimap2 alignments');
    
    # cleanup old files
    for my $file("$minimap2PAF.tmp", "$minimap2IDX.tmp", $minimap2IDX){
        if (-s $file){
            unlink $file;
        }
    }
    
    ### run minimap2 ###
    # open a pipe to the minimap2 job
    my $errLog = "$ERR_DIR/minimap2.stderr";
    my $pipeCmd = "minimap2 $minimap2_parameters -t $threads $assembly_index - 2> $errLog  > $minimap2PAF.tmp";
    open my $MMA, '|-', $pipeCmd or err("failed to open minimap2 pipe command $pipeCmd");
    
    # sig handler
    local $SIG{PIPE} = sub { 
        err("Minimap2 pipe has died:\n(this script) | $pipeCmd\n\nCheck $errLog for possible cause.") 
    };
    
    # pass seqs to pipe job
    for my $contig (sort keys %contigLEN){
        if ($suspects{$contig}){
            print $MMA get_seq($contig);
        }
    }
    close $MMA;
    
    # finalize the alignment output
    rename "$minimap2PAF.tmp", $minimap2PAF;
    
    ### build the index of the alignments ###
    msg('Indexing minimap2 alignments');
    
    # open the alignments for reading
    open my $FWR, '<', $minimap2PAF or err("failed to open $minimap2PAF for reading");
    
    # open filehandle for alignment index
    open my $IDX, '>', "$minimap2IDX.tmp" or err("failed to open $minimap2IDX.tmp for writing");
    
    # tracking for index file
    my $curCtg = 'init';
    my $offset = tell($FWR);
    
    # generate the index file
    while(<$FWR>){
        my ($ctg) = $_ =~ m/^(.+?)\s/;
        if ($curCtg ne $ctg){
            print $IDX $ctg, "\t", $offset, "\n";
            $curCtg = $ctg;
        }
        $offset = tell($FWR);
    }
    close $FWR;
    close $IDX;
    
    # finalize the index file
    rename "$minimap2IDX.tmp", $minimap2IDX;
    
    return;
}



sub read_minimap2_index {
    msg('Reading index of minimap2 alignments');
    
    open my $IDX, '<', $minimap2IDX or err("failed to open $minimap2IDX for reading");
    
    # NOTE: read the index faster
    # ---------------------------
    
	while(<$IDX>){
		chomp;
		if ( $_ =~ m/^(\S+)\s+(\S+)/ ) {
        
			$mmIndex{$1} = $2;
        
        }
    }
    
    # ORIG:

#     while(<$IDX>){
#         my@l=split/\s+/;
#         $mmIndex{$l[0]} = $l[1];
#     }

    close $IDX;
    
    return;
}


sub hit_summary {

    # NOTE: We never enter this method in this modified version of the script
    # -----------------------------------------------------------------------

    # check if needed
    if (-s $hit_summary){
        msg('Reusing hit summary file from previous run');
    } else {
    
        msg('Preparing contig hit summary');
        
        # cleanup
        if (-s "$hit_summary.tmp"){
            unlink "$hit_summary.tmp";
        }
        
        # queue up suspect contigs
        for my $contig (sort { $contigLEN{$b} <=> $contigLEN{$a} } keys %suspects){
            $queueHitSummary->enqueue($contig);
        }
        
        # finalise queue
        $queueHitSummary->end();
        
        # spawn worker threads for hit summary
        for (1..$threads){
            $available_threads->down(1);
            threads->create(\&hit_summary_worker);
        }
        
        # wait on remaining jobs
        $available_threads->down($threads);
        $available_threads->up($threads);
        
        # join workers
        for my $thr (threads->list()){
            $thr->join();
        }
        
        # finalise the hit summary
        if (-s "$hit_summary.tmp"){
            rename "$hit_summary.tmp", $hit_summary;
        } else {
            if (!(keys %junk)){
                msg('WARNING: There were no hits for suspect contigs and no flagged artefact contigs');
                msg('WARNING: Nothing left to do, exiting...');
                exit(0);
            } else {
                msg('WARNING: There were no hits for suspect contigs');
                msg('WARNING: Dropping flagged artefact contigs and exiting...');
                write_assembly();
                exit(0);
            }
        }
    }
    
    # slurp the hits
    open my $TSV, '<', $hit_summary or err("failed to open $hit_summary for reading");
    
	# NOTE: read the hit summary into a more memory efficient string
    # --------------------------------------------------------------
    
    # NOTE: this changes the data structure from an array to a concatenated string
    # which requires splitting later on.
    
    while(my $l = <$TSV>){
        
        if ( $l =~ m/^(\S+)\s+(\S+)/ ) {
        
			my ( $ctg , $other_ctg ) = ( $1 , $2 );
			
			$hits{ $ctg } .= "|$other_ctg";
        
        }
        
    }
    close $TSV;
    
    foreach my $ctg ( keys %hits ) {
    
		substr ( $hits{ $ctg } , 0 , 1 , "" );
    
    }
    
    my $n_contigs = keys %hits;
    
	msg("Read hit summaries for $n_contigs contigs");
    
    # --------------------------------------------------------------
    
    # ORIG:
    
#     while(my $l = <$TSV>){
#         my @line = split /\s+/, $l;
#         if (!($hits{$line[0]})){
#             my @h :shared;
#             $hits{$line[0]} = \@h;
#         }
#         push @{$hits{$line[0]}}, $line[1];
#     }
#     close $TSV;

    return;
}



sub hit_summary_worker {
    # iterate over the queue until no more jobs
    while (defined(my $ctg = $queueHitSummary->dequeue())) {
        
        # skip if not alignments for contig
        if (!(defined($mmIndex{$ctg}))){
            next;
        }
        
        # collect hit contig total matching residues and print PAF
        my %hits;
        
        # open the alignment file for reading
        open my $FWR, '<', $minimap2PAF or err("failed to open $minimap2PAF for reading");
        
        # fast-forward to to the correct location in the alignment file
        seek($FWR, $mmIndex{$ctg}, 0);

        # NOTE: grab the hits with a fast rexexp without splitting into an array
        # ----------------------------------------------------------------------
        
        # grab the hits
        while(<$FWR>){
            
			if ( $_ =~ m/^(\S+)\s+\d+\s+\d+\s+\d+\s+\S+\s+(\S+)\s+\d+\s+\d+\s+\d+\s+(\d+)/ ) {
			
				# q contig = $1
				# t contig = $2
				# matches = $3
				
				last if ($1 ne $ctg);
				if ($1 ne $2){
					$hits{$2}+=$3;
				}
				
			}
        
        }
        
        # -----------------------------------------------------
        
        # ORIG:
        
#         # grab the hits
#         while(<$FWR>){
#             my@l=split/\s+/;
#             last if ($l[0] ne $ctg);
#             if ($l[0] ne $l[5]){
#                 $hits{$l[5]}+=$l[9];
#             }
#         }
        
        # filter the hits
        for my $hit (keys %hits){
			my $val = ($hits{$hit} / $contigLEN{$ctg}) * 100;
			msg("$hit for $ctg: $val");
            if ($val * 100 < $hit_cutoff){ 
                delete $hits{$hit};
            }
#             else {
# 				msg("Keeping $hit for $ctg: $val");
#             }
        }
        
        # prepare the output
        my @out;
        for my $outHit (sort { $hits{$b} <=> $hits{$a} } keys %hits){
                push @out, "$ctg\t$outHit\n";
        }
        
        # write the hits to the summary file
        if (@out){
            $writing_to_out->down(1);
            open my $HIT, '>>', "$hit_summary.tmp" or err("failed to open $hit_summary.tmp for appending");
            print $HIT @out;
            close $HIT;
            $writing_to_out->up(1);
        }
    }
    
    $available_threads->up(1);
    
    return;
}



sub parse_repeats {
    
    if (-s $assembly_repeats){
        msg('Reusing repeats from tmp directory');
    } else {
        msg("Parsing repeat annotations from $repeats");
        
        # we'll pass this through bedtools as a format test and to ensure overlaps are merged
        runcmd({ command => "cat $repeats | sort -k1,1 -k2,2n -k3,3n | bedtools merge > $assembly_repeats.tmp 2> $ERR_DIR/bedtools_merge.stderr",
                 logfile => "$ERR_DIR/bedtools_merge.stderr",
                 silent => 1 });
        
        # finalize
        rename "$assembly_repeats.tmp", $assembly_repeats;
    }
    
    # open the newely generated repeats file for indexing
    open my $BED, '<', $assembly_repeats or err("Failed to open $assembly_repeats for reading");
    
    # iterate repeats
    my $ctg = 'init';
    my $offset = tell($BED);
    
    while(<$BED>){
        next if (/^#/); 
        next if (/^\s+$/);
        my @l = split /\s+/;
        
        # check if new contig
        if ($ctg ne $l[0]){
            if ($contigLEN{$l[0]}){
                $repIndex{$l[0]} = $offset;
            }
            $ctg = $l[0];
        }
        
        $offset = tell($BED);
        
        # subtract mask interval from the contig len
        $contigMLEN{$l[0]} -= ($l[2] - $l[1]);
        
    }
    close $BED;
    
    # NOTE: read the repeats file once and populate the repeats table in memory
    # -------------------------------------------------------------------------
    
    unless ( keys %repeat_table ) {
    
		msg("Reading repeats in $assembly_repeats ...");
    
		# open the repeats file
		open my $BED, '<', $assembly_repeats or err("Failed to open $assembly_repeats for reading");
		
		while (<$BED>) {
			
			if ( $_ =~ m/^(\S+)\s/ ) {
			
				$repeat_table{ $1 } .= $_;
			
			}
		
		}
    
    }
    
    return;
}


sub simplify_paf {

	# NOTE: simplify PAF by keeping hits above the treshold
	# -----------------------------------------------------

	my @ctgs_w_repeats = keys %repeat_table;
	my $n_ctgs_w_repeats = @ctgs_w_repeats;
	
	msg("We have repeat intervals for $n_ctgs_w_repeats contigs...");
    
	msg("Simplifying PAF $minimap2PAF into $minimap2PAF_simple");
		
	# Make a log
	
	open ( my $paf_out_info , ">" , "$minimap2PAF_simple.info.tsv" ) or die "$!";
	
	# Read the PAF one time to find the candidate contigs
	
	open( my $paf_in , "$minimap2PAF" ) or die "$!";
	
	my $last_ctg;
	my $last_ctg_length;
	my $last_ctg_length_after_repeats;
	my $t_ctg_length;

	my %tmp_hits;
	
	my $valid_hits = {};
	my %suspects_simple;
	
	my $i = 1;
	
	while (<$paf_in>) {

		chomp;
		
		# ctg1	638930	4	638927	+	ctg1	638930	4	638927	549643	638923	60	tp:A:P	cm:i:94108	s1:i:549643	s2:i:5133	dv:f:0	rl:i:188611
		
		if ( $_ =~ m/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)/ ) {
		
			my ( $q_ctg , $q_ctg_length , $q_start , $q_stop , $t_ctg , $t_ctg_length , $t_start , $t_stop ) =
			( $1 , $2 , $3 , $4 , $5 , $6 , $7 , $8 );
			
			if ( defined $last_ctg and ( $q_ctg ne $last_ctg ) ) {
				
				if ( keys %tmp_hits ) {
				
					simplify_paf_finish_ctg(
						$i ,
						$last_ctg ,
						$last_ctg_length ,
						\%tmp_hits ,
						$last_ctg_length_after_repeats ,
						$paf_out_info ,
						$valid_hits ,
						\%suspects_simple ,
					);
				
				}
				
				$i++;
				
				%tmp_hits = ();
				
				$last_ctg = undef;
				$last_ctg_length = undef;
				$last_ctg_length_after_repeats = undef;
				$t_ctg_length = undef;
				
			}
			
			# Check query contig
			# ------------------
			
			if ( defined $suspects{ $q_ctg } ) {
			
				# Query
			
				unless ( defined $last_ctg_length_after_repeats ) {

					# Check repeats
					
					if ( defined $length_table_after_repeats{ $q_ctg } ) {
					
						$last_ctg_length_after_repeats = $length_table_after_repeats{ $q_ctg };
					
					}
					
					else {
					
						if ( defined $repeat_table{ $q_ctg } ) {
						
							my $tmp_seq = 1 x $q_ctg_length;
						
							foreach my $repeat ( split ( /\n/ , $repeat_table{ $q_ctg } ) ) {
							
								if ( $repeat =~ m/^\S+\s(\d+)\s(\d+)/ ) {
								
									my $length = $2 - $1 + 1;
									
									substr ( $tmp_seq , $1 , $length , 2 x $length );
								
								}
							
							}
							
							my $repeats = $tmp_seq =~ tr/2//;
					
							$last_ctg_length_after_repeats = $q_ctg_length - $repeats;
						
						}
						
						else {
						
							$last_ctg_length_after_repeats = $q_ctg_length;
						
						}
						
						$length_table_after_repeats{ $q_ctg } = $last_ctg_length_after_repeats;
				
					}
				
				}

				$last_ctg = $q_ctg;
				$last_ctg_length = $q_ctg_length;
				
				if ( $q_ctg eq $t_ctg ) {
				
					next;
				
				}
			
				my $t_prop = ( $t_ctg_length / $last_ctg_length_after_repeats ) * 100;
					
				my $match_length = $q_stop - $q_start + 1;
				
				$tmp_hits{ $t_ctg }->{ 'match_length' } += $match_length;
				
				push @{ $tmp_hits{ $t_ctg }->{ 'matches' } } , [ $q_start , $q_stop , $match_length ];
				
			}
			
		}
		
	}

	if ( keys %tmp_hits ) {

		simplify_paf_finish_ctg(
			$i ,
			$last_ctg ,
			$last_ctg_length ,
			\%tmp_hits ,
			$last_ctg_length_after_repeats ,
			$paf_out_info ,
			$valid_hits ,
			\%suspects_simple ,
		);

	}
	
	close $paf_in;
	
	# Output files
	
    # check if needed
    
    my $hit_summary_out;
    
    if (-s $hit_summary){
        msg('Reusing hit summary file from previous run');
    } else {
    
        msg('Preparing contig hit summary');
        
        # cleanup
        if (-s "$hit_summary.tmp"){
            unlink "$hit_summary.tmp";
            
        }

        open ( $hit_summary_out , ">" , "$hit_summary.tmp" ) or die "$!";
        
	}

	my $check_hits_table = {};
	
	foreach my $q_ctg ( sort keys %suspects_simple ) {
		
		my @t_ctgs = sort { $valid_hits->{ $q_ctg }->{ $b } <=> $valid_hits->{ $q_ctg }->{ $a } } keys %{ $valid_hits->{ $q_ctg } };
		
		foreach my $t_ctg ( @t_ctgs ) {
		
			if ( $hit_summary_out ) {
		
				say $hit_summary_out "$q_ctg\t$t_ctg";
			
			}
			
			$check_hits_table->{ $q_ctg }->{ $t_ctg } = 1;
			$check_hits_table->{ $t_ctg }->{ $q_ctg } = 1;
			$check_hits_table->{ $t_ctg }->{ $t_ctg } = 1;
			
		}
		
		$check_hits_table->{ $q_ctg }->{ $q_ctg } = 1;
	
	}
	
	rename "$hit_summary.tmp", $hit_summary;
	
	%hits = ();

	open ( my $in , "<" , $hit_summary ) or die "$!";

    while(my $l = <$in>){
        
        if ( $l =~ m/^(\S+)\s+(\S+)/ ) {
        
			my ( $ctg , $other_ctg ) = ( $1 , $2 );
			
			$hits{ $ctg } .= "|$other_ctg";
        
        }
        
    }
    close $in;

    foreach my $ctg ( keys %hits ) {

		substr ( $hits{ $ctg } , 0 , 1 , "" );

    }

    my $n_contigs = keys %hits;

	msg("Read hit summaries for $n_contigs contigs");
	
	# Read the PAF a second time to filter it
	
	open( $paf_in , "$minimap2PAF") or die "$!";
	
	# Make the simplified PAF
	
	open ( my $paf_out , ">" , "$minimap2PAF_simple.tmp" ) or die "$!";
	
	$i = 1;
	
	while (<$paf_in>) {

		chomp;
		
		# ctg1	638930	4	638927	+	ctg1	638930	4	638927	549643	638923	60	tp:A:P	cm:i:94108	s1:i:549643	s2:i:5133	dv:f:0	rl:i:188611
		
		if ( $_ =~ m/^(\S+)\s+\d+\s+\d+\s+\d+\s+\S+\s+(\S+)\s+\d+\s+\d+\s+\d+/ ) {
		
			my ( $q_ctg ,  $t_ctg ) = ( $1 , $2 );
			
			if ( defined $check_hits_table->{ $q_ctg }->{ $t_ctg } ) {
			
				msg("SECOND PASS\t$q_ctg\t$t_ctg");
			
				say $paf_out "$_";
			
			}
			
		}
		
	}
	
	%tmp_hits = ();
	
	# finalize the alignment output
	rename "$minimap2PAF_simple.tmp", $minimap2PAF_simple;
	
	### build the index of the alignments ###
	msg('Indexing simplified minimap2 alignments');
	
	# open the alignments for reading
	open my $FWR, '<', $minimap2PAF_simple or err("failed to open $minimap2PAF_simple for reading");
	
	# open filehandle for alignment index
	open my $IDX, '>', "$minimap2IDX_simple.tmp" or err("failed to open $minimap2IDX_simple.tmp for writing");
	
	# tracking for index file
	my $curCtg = 'init';
	my $offset = tell($FWR);
	
	# generate the index file
	while(<$FWR>){
		my ($ctg) = $_ =~ m/^(.+?)\s/;
		if ($curCtg ne $ctg){
			print $IDX $ctg, "\t", $offset, "\n";
			$curCtg = $ctg;
		}
		$offset = tell($FWR);
	}
	close $FWR;
	close $IDX;
	
# 	# finalize the index file
	rename "$minimap2IDX_simple.tmp", $minimap2IDX_simple;
	
	return;
	
}

sub simplify_paf_finish_ctg {

	my ( $i , $last_ctg , $last_ctg_length , $tmp_hits_ref , $last_ctg_length_after_repeats , $paf_out_info , $valid_hits , $suspects_simple ) = @_;
	
	msg("Checking contig $i $last_ctg ...");
	
	# my $local_bestmatch_cutoff = $bestmatch_cutoff - 10;
	
	foreach my $t_ctg ( keys %{ $tmp_hits_ref } ) {
		
		# Check query contig
		# ------------------
		
		# Only bother to check if total sum of alignments are longer than the treshold
		# Otherwise there is no point in continuing
	
		if ( ( $tmp_hits_ref->{ $t_ctg }->{ 'match_length' } / $last_ctg_length_after_repeats ) * 100 >= $bestmatch_cutoff ) {

			# Count aligned bases more carefully
		
			my $seq = 0 x $last_ctg_length;
		
			foreach my $match_ref ( @{ $tmp_hits_ref->{ $t_ctg }->{ 'matches' } } ) {
				
				substr ( $seq , $match_ref->[ 0 ] , $match_ref->[ 2 ] , 1 x $match_ref->[ 2 ] );
					
			}
	
			my $matches = $seq =~ tr/1//;
			
			# Only bother to check this alignment if it has the potential to be long enough
			# Otherwise there is no point in continuing
			
			if ( ( $matches / $last_ctg_length_after_repeats ) * 100 >= $bestmatch_cutoff ) {
			
				# Remove matches that overlap with repeats
				
				if ( defined $repeat_table{ $last_ctg } ) {
				
					foreach my $repeat ( split ( /\n/ , $repeat_table{ $last_ctg } ) ) {
					
						if ( $repeat =~ m/^\S+\s(\d+)\s(\d+)/ ) {
						
							my $length = $2 - $1 + 1;
							
							substr ( $seq , $1 , $length , 2 x $length );
						
						}
					
					}
				
				}
				
				my $matches_after_repeats = $seq =~ tr/1//;
				
				# Keep target contig if proportion of matched sites is above
				# treshold after removing repeat regions
				
				my $aligned = ( $matches_after_repeats / $last_ctg_length_after_repeats ) * 100;

				if ( $aligned >= $bestmatch_cutoff ) {
					
					print $paf_out_info "HIT\t$last_ctg\t$last_ctg_length\t$last_ctg_length_after_repeats\t$t_ctg\t$matches\t$matches_after_repeats\t$aligned\n";
					
					msg("HIT\t$last_ctg\t$last_ctg_length\t$last_ctg_length_after_repeats\t$t_ctg\t$matches\t$matches_after_repeats\t$aligned");
					
					$valid_hits->{ $last_ctg }->{ $t_ctg } = $aligned;
					
					$suspects_simple->{ $last_ctg } = 1;
					
				}
			
			}
		
		}
		
	}
	
}

sub get_reps {
    my $ctg = $_[0];
    my $bed;

	# NOTE: just return the pre-read repeat intervals
	# -----------------------------------------------
	
    if ( defined $repeat_table{ $ctg } ) {
    
		$bed = $repeat_table{ $ctg };
    
    }
    
    # ORIG:
    
#     # open the repeats file
#     open my $REP, '<', $assembly_repeats or err("Failed to open $assembly_repeats for reading");
#     
#     # fast forward to contig
#     if (defined($repIndex{$ctg})){
#         seek($REP, $repIndex{$ctg}, 0);
#         
#         # read the repeat annotations
#         while(<$REP>){
#             my@l=split/\s+/;
#             if($l[0] eq $ctg){
#                 $bed .= $_;
#             } else {
#                 last;
#             }
#         }
#     }

    # return the bed file as a string
    return $bed;
}


#---ITERATIVE STEP---



sub get_contig_hits {
    msg('Reading contig hits from hit summary');
    
    # remove reference contigs if they themselves have been reassigned 
    for my $ctg (keys %suspects){
        next if ($contigREASSIGN{$ctg});
        if ($contigHIT1{$ctg}){
            if ($contigREASSIGN{$contigHIT1{$ctg}}){
                undef $contigHIT1{$ctg};
                undef $contigHIT2{$ctg};
                $contigASSIGN{$ctg} = 0;
            }
        }
        if ($contigHIT2{$ctg}){
            if ($contigREASSIGN{$contigHIT2{$ctg}}){
                undef $contigHIT2{$ctg};
                $contigASSIGN{$ctg} = 0;
            }
        }
    }
    
    # NOTE: add counter here
    # ----------------------
    
    my $i = 0;
    # add best non-reassigned reference contigs for each contig 
    for my $ctg (keys %suspects){
        if ($hits{$ctg}){
			# NOTE: 
			msg("\tContig $i") unless $i % 100_000;
			for my $hit (split ( /\|/ , $hits{$ctg} ) ) {
            # for my $hit (@{$hits{$ctg}}){
                if (!($contigREASSIGN{$ctg}) && !($contigREASSIGN{$hit})){
                    if (!($contigHIT1{$ctg})){
                        $contigHIT1{$ctg} = $hit;
                        $contigASSIGN{$ctg} = 0;
                    } elsif (!($contigHIT2{$ctg}) && ($hit ne $contigHIT1{$ctg})){
                        $contigHIT2{$ctg} = $hit;
                        $contigASSIGN{$ctg} = 0;
                    }
                }
            }
            $i++;
        }
    }
    
    return;
}



sub pairwise_alignments {
    msg('Performing pairwise comparisons on contig hits');
    
    $queuePairwise = Thread::Queue->new();
    
    # NOTE: reset threat settings
    # ---------------------------
    
	$available_threads = Thread::Semaphore->new($threads);
	$max_jobs = Thread::Semaphore->new($limit_io);
    
    CTG: for my $contig (sort { $contigLEN{$b} <=> $contigLEN{$a} } keys %suspects){
        next CTG if ($contigREASSIGN{$contig});
        next CTG if ($contigASSIGN{$contig});
        
        # skip if no hits
        if (!($contigHIT1{$contig})){
            $contigASSIGN{$contig} = 'n';
            $contigBM{$contig} = '-';
            $contigMM{$contig} = '-';
            next CTG;
        }
        
        # queue the pairwise comparison
        $queuePairwise->enqueue(['purge',$contig]);
        
    }
    
    # finalise queue
    $queuePairwise->end();
    
    # spawn worker threads
    for (1..$limit_io){
        $available_threads->down(1);
        threads->create(\&pairwise_worker);
    }
    
    # wait on remaining jobs
    $available_threads->down($threads);
    $available_threads->up($threads);

    # join workers
    for my $thr (threads->list()){
        $thr->join();
    }
    
    return;
}



sub pairwise_worker {
    
    while (defined(my $args = $queuePairwise->dequeue())) {
        my $job = @{$args}[0];
        my $contig = @{$args}[1];
        
        my $alignments;
        my $bestmatch=0;
        my $maxmatch=0;
        my $assignment;
        
        my $ref1;
        my $ref2;
        
        $ref1 = $contigHIT1{$contig};
        $ref2 = $contigHIT2{$contig} if ($contigHIT2{$contig});
        
        # check if the contig is 100% repeat
        if ($contigMLEN{$contig} == 0){
            $assignment = "r";
        } else {

            # get minimap2 alignments for ref1 and ref2
            ($ref1, $ref2, $alignments) = minimap2_alignments($ref1, $ref2, $contig, $job);
            
            # get 'maxmatch' and 'bestmatch' coverages from the sorted general output of minimap2, subtract masking regions on the fly
            ($bestmatch, $maxmatch) = get_bm_mm(\@$alignments, $contig);

            msg("$contig vs $ref1: $bestmatch, $maxmatch\n"); # NOTE: debug
            
            # guess the assignment
            $assignment = guess_assignment($bestmatch, $maxmatch);
        }
        
        # print the reassignments
        $writing_to_out->down(1);
        $contigASSIGN{$contig} = $assignment;
        $contigBM{$contig} = $bestmatch;
        $contigMM{$contig} = $maxmatch;
        $writing_to_out->up(1);
        
        # NOTE: print a message here
        # -----------------------------
        
		$n_contigs_done++;
		
		print "$n_contigs_done contigs done!\n";
        
    }

    $available_threads->up(1);
    
    return;
}



sub minimap2_alignments {
    my $ref1 = $_[0];
    my $ref2 = $_[1];
    my $contig = $_[2];
    my $job = $_[3];
    
    my @r1paf;
    my @r2paf;
    my @rbed;
    
    # NOTE: skip this variable
    # ---------------------------

    # my @rsbed;
    
    # NOTE: read PAF without splitting line into array
    # ------------------------------------------------
    
    # alignments
    open my $FWR, '<', $minimap2PAF or err("failed to open $minimap2PAF for reading");
    seek($FWR, $mmIndex{$contig}, 0);
    while(<$FWR>){
		if ( $_ =~ m/^(\S+)\s+\d+\s+\d+\s+\d+\s+\S+\s+(\S+)\s+\d+/ ) {
			last if ($1 ne $contig);
			if ($2 eq $ref1){
				push @r1paf, $_;
			} elsif (($ref2) && ($2 eq $ref2)){
				push @r2paf, $_;
			}
		}
    }
    close $FWR;
    
    # ------------------------------------------------
    
    # ORIG:
    
#     # alignments
#     open my $FWR, '<', $minimap2PAF or err("failed to open $minimap2PAF for reading");
#     seek($FWR, $mmIndex{$contig}, 0);
#     while(<$FWR>){
#         my@l=split/\s+/;
#         last if ($l[0] ne $contig);
#         if ($l[5] eq $ref1){
#             push @r1paf, $_;
#         } elsif (($ref2) && ($l[5] eq $ref2)){
#             push @r2paf, $_;
#         }
#     }
#     close $FWR;

    if (!(@r1paf)){
        $ref1 = 0;
    }
    if($ref2){
        if (!(@r2paf)){
            $ref2 = 0;
        }
    }
    
    # make rdotplot files for plotting dotplots
    if (($dotplots)&&($job eq 'recheck')){
        my $TDP;
        if ($ref1){
            open $TDP, '>', "$TMP_ALN/$contig.1.rdotplot" or err("failed to open $TMP_ALN/$contig.1.rdotplot for writing");
            for my $line (@r1paf){
                my @l=split/\s+/,$line;
                if ($l[4] eq '+'){
                    print $TDP "$l[2]\t$l[7]\n$l[3]\t$l[8]\nNA\tNA\n";
                } else {
                    print $TDP "$l[3]\t$l[7]\n$l[2]\t$l[8]\nNA\tNA\n";
                }
            }
            close $TDP;
        }
        if ($ref2){
            open $TDP, '>', "$TMP_ALN/$contig.2.rdotplot" or err("failed to open $TMP_ALN/$contig.2.rdotplot for writing");
            for my $line (@r2paf){
                my @l=split/\s+/,$line;
                if ($l[4] eq '+'){
                    print $TDP "$l[2]\t$l[7]\n$l[3]\t$l[8]\nNA\tNA\n";
                } else {
                    print $TDP "$l[3]\t$l[7]\n$l[2]\t$l[8]\nNA\tNA\n";
                }
            }
            close $TDP;
        }
    }

    # NOTE: read PAF without splitting line into array
    # ---------------------------------------------------
    
	# only need the contig coordinates
    if ($ref1){
        for my $line (@r1paf){
			if ( $line =~ m/^\S+\s+\d+\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s+\d+/ ) {
				push @rbed, [$1,$2];
            }
        }
    }
    
    if ($ref2){
        for my $line (@r2paf){
			if ( $line =~ m/^\S+\s+\d+\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s+\d+/ ) {
				push @rbed, [$1,$2];
            }
        }
    }
    
    # ------------------------------------------------
    
    # ORIG:
    
#     # only need the contig coordinates
#     if ($ref1){
#         for my $line (@r1paf){
#             my @l=split/\s+/,$line;
#             push @rbed, [$l[2],$l[3]];
#         }
#     }
#     
#     if ($ref2){
#         for my $line (@r2paf){
#             my @l=split/\s+/,$line;
#             push @rbed, [$l[2],$l[3]];
#         }
#     }
    
    # sort
    @rbed = sort { $a->[0] <=> $b->[0] } @rbed;
    
    # NOTE: faster return from repeat handling
    # ----------------------------------------
    
    # subtract repeats if needed
    if ($repeats){
        return ( $ref1, $ref2, sub_repts(\@rbed, $contig) );
    } else {
        return ( $ref1, $ref2, \@rbed );
    }
    
    # ORIG:
    
#     # subtract repeats if needed
#     if ($repeats){
#         @rsbed = sub_repts(\@rbed, $contig);
#     } else {
#         @rsbed = @rbed;
#     } 
#    
#     return ($ref1, $ref2, \@rsbed);

}


# NOTE: this rewritten sub_repts method avoids the external call to bedtools to get repeat info
# ---------------------------------------------------------------------------------------------

sub sub_repts {

    my $bed = $_[0];
    my $contig = $_[1];
    
	# get sequence length
    
    my $tmp_seq_length = $contigLEN{ $contig };
    
	# make a sequence mask
    
    my $tmp_seq = 0 x $tmp_seq_length;
    
    foreach my $al ( @{ $bed } ) {
    
		my ( $start , $stop ) = ( $al->[ 0 ] , $al->[ 1 ] );
    
		my $al_length = $stop - $start + 1;
		
		# Mark alignments
		
		substr ( $tmp_seq , $start , $al_length , 1 x $al_length );
    
    }
    
    # get the repeat annotations
    my $repAnnotations = get_reps($contig);
    # my $repAnnotations;
    
    if ( $repAnnotations ) {
        
        foreach my $rep_line ( split ( /\n/ , $repAnnotations ) ) {
        
			if ( $rep_line =~ m/^\S+\s(\d+)\s(\d+)/ ) {
			
				my $rep_length = $2 - $1 + 1;
				
				# Reset repeat regions
				
				substr ( $tmp_seq , $1 , $rep_length , 0 x $rep_length );
			
			}
        
        }
        
	}
	
	my @out;
	
	while ( $tmp_seq =~ m/1+/g ) {
	
		my ( $start , $stop ) = ( @- , @+ );
		$stop--;
	
		push @out , [ $start , $stop ];
	
	}

    # NOTE: return array ref instead of array
    # ---------------------------------------
    
    return \@out;

}


sub sub_repts_bedtools_original {
    my $bed = $_[0];
    my $contig = $_[1];
    
    my $rep = "$TMP_ALN/$contig.bed";
    my $ali = "$TMP_ALN/$contig.ali";
    
    # get the repeat annotations
    my $repAnnotations = get_reps($contig);
    
    my @out;
    
    if ($repAnnotations){
        
        # dump repeats to tmp file
        open my $REP, '>', $rep or err("Failed to open $rep for writing");
        print $REP $repAnnotations;
        close $REP;
        
        # dump alignment coords as tmp bed file
        open my $TMP, '>', $ali or err("Failed to open $ali for writing");
        for my $al (@{$bed}){
            print $TMP "$contig\t@{$al}[0]\t@{$al}[1]\n";
        }
        close $TMP;
        
        # subtract alignments
        my $pipeCmd = "bedtools subtract -a '$ali' -b '$rep'";
        open $TMP, '-|', $pipeCmd or err("Failed to open pipe: $pipeCmd");
        
        # NOTE: faster reading without array splitting
        # --------------------------------------------
        
        while(<$TMP>){
			chomp;
			if ( $_ =~ m/^\S+\s(\S+)\s(\S+)/ ) {
				# print "Remaining: $1\t$2\n";
				push @out, [ $1 , $2 ];
            }
        }
        close $TMP or err("Failed to close pipe: $pipeCmd | (this script)");
        
        # --------------------------------------------
        
        # ORIG:
        
#         while(<$TMP>){
#             my@l=split/\s+/;
#             push @out, [ $l[1],$l[2] ];
#         }

        # cleanup
        for my $file ($rep, $ali){
            if (-e $file){
                unlink $file;
            }
        }
        
    } else {
        @out = @{$bed};
    }
    
    # return @out;
    

    # NOTE: return array ref instead of array
    # ---------------------------------------
    
    return \@out;
    
}



sub get_bm_mm {
    my @bed = @{$_[0]};
    my $contig = $_[1];
    
    my $bm=0;
    my $mm=0;
    
    
    my @p;
    
	# NOTE: less array dereferencing and accessing
	# --------------------------------------------

    for my $l (@bed){
		my ( $start , $stop ) = ( $l->[0] , $l->[1] );
        $mm+=($stop - $start);
        if (@p){
            next if (($start > $p[0]) && ($stop < $p[1]));
            if ($start > $p[1]){
                $bm+=($p[1]-$p[0]);
                @p=($start, $stop);
            } elsif ($p[1] < $stop) {
                $p[1] = $stop;
            } 
        } else {
            @p=($start, $stop);
        }
    }
    
    # --------------------------------------------
    
    # ORIG:
    # --------------------------------------------

#     for my $l (@bed){
#         $mm+=($l->[1] - $l->[0]);
#         if (@p){
#             next if (($l->[0] > $p[0]) && ($l->[1] < $p[1]));
#             if ($l->[0] > $p[1]){
#                 $bm+=($p[1]-$p[0]);
#                 @p=($l->[0], $l->[1]);
#             } elsif ($p[1] < $l->[1]) {
#                 $p[1] = $l->[1];
#             } 
#         } else {
#             @p=($l->[0], $l->[1]);
#         }
#     }

	# --------------------------------------------

    $bm+=($p[1]-$p[0]) if (@p);

    if ($repeats){
        $mm = sprintf "%.2f", ($mm/$contigMLEN{$contig}) * 100;
        $bm = sprintf "%.2f", ($bm/$contigMLEN{$contig}) * 100;
    } else {
        $mm = sprintf "%.2f", ($mm/$contigLEN{$contig}) * 100;
        $bm = sprintf "%.2f", ($bm/$contigLEN{$contig}) * 100;
    }
    
    return ($bm, $mm);
}



sub guess_assignment {
    my $assignment;
    if ($_[0] >= $bestmatch_cutoff){
        if ($_[1] >= $maxmatch_cutoff){
            $assignment = "r";
        } else {
            $assignment = "h";
        }
    } elsif ($_[0] < $low_cutoff){
        $assignment = "n";
    } else {
        $assignment = "u";
    }
    return $assignment;
}



sub check_assignments {
    msg('Checking contig assignments for conflicts');
    
    # check all assignments for conflicts
    for my $ctg (sort keys %suspects){
        next if ($contigASSIGN{$ctg} !~ /[rh]/i);
        next if ($contigREASSIGN{$ctg});
        if ($contigREASSIGN{$contigHIT1{$ctg}}){
            err('ref seq was already reassigned, this should not have happened');
        }
        
        my $r_ctg = $contigHIT1{$ctg};
        if (($contigASSIGN{$r_ctg})&&($contigASSIGN{$r_ctg} =~ /[rh]/i)){
            msg("CONFLICT: $ctg and it's match $r_ctg both flagged for reassignment");
            
            # just keep the longer contig
            if ($contigLEN{$ctg} > $contigLEN{$r_ctg}){
                $contigASSIGN{$ctg} = 0;
                msg("\tKeeping longer contig $ctg");
            } else {
                $contigASSIGN{$r_ctg} = 0;
                msg("\tKeeping longer contig $r_ctg");
            }
        }
    }
    return;
}



sub add_reassignments {
    msg('Logging reassignments and checking for convergence');
    
    my $convergence_check = 1;
    
    for my $ctg(sort keys %suspects){
        next if ($contigREASSIGN{$ctg});
        
        if ($contigASSIGN{$ctg} =~ /[rh]/i){
            $contigREASSIGN{$ctg} = 1;
            $convergence_check = 0;
        }
    }
    
    # convergence check
    if ($convergence_check){
        $convergence = 1;
        msg('Convergence reached!');
    } else {
        msg('Convergence not reached, more passes needed');
    }
    
    return;
}



sub over_purge_check {
    msg('Checking for over-purging');
    
    # first round is to detect over-purged contigs
    $over_purge_mode = 'recheck';
    
    $queueOverPurge = Thread::Queue->new();
    
    # NOTE: reset threat settings
    # ---------------------------
    
	$available_threads = Thread::Semaphore->new($threads);
	$max_jobs = Thread::Semaphore->new($limit_io);
    
    # queue up all reassigned contigs
    for my $ctg(sort keys %suspects){
        if ($contigREASSIGN{$ctg}){
            $queueOverPurge->enqueue($ctg);
        }
    }
    
    # finalise queue
    $queueOverPurge->end();
    
    # spawn workers to check if still satisfy conditions of haplotig
    for (1..$limit_io){
        $available_threads->down(1);
        threads->create(\&over_purge_worker);
    }
    
    # wait on remaining jobs
    $available_threads->down($threads);
    $available_threads->up($threads);

    # join workers
    for my $thr (threads->list()){
        $thr->join();
    }
    
    msg('Fixing over-purged contigs');
    
    # second round is to iteratively fix over-purged contigs
    $over_purge_mode = 'fix';
    
    $queueOverPurge = Thread::Queue->new();
    
    # queue up all over-purged contigs
    for my $ctg(sort { $contigLEN{$b} <=> $contigLEN{$a} } keys %contigOVPURG){
        $queueOverPurge->enqueue($ctg);
    }
    
    # finalise queue
    $queueOverPurge->end();
    
    # fix contigs
    $available_threads->down(1);
    over_purge_worker();
    
    return;
}



sub over_purge_worker {
    while (defined(my $ctg = $queueOverPurge->dequeue())) {
        
        # add hits for contig
        undef $contigHIT1{$ctg};
        undef $contigHIT2{$ctg};
        my $alignments;
        
        # NOTE: split string instead of using array as in the original code
        
        # current best hits for contig
        for my $hit (split ( /\|/ , $hits{$ctg} ) ) {
        # for my $hit (@{$hits{$ctg}}){
            if (!($contigREASSIGN{$hit})){
                if (!($contigHIT1{$ctg})){
                    $contigHIT1{$ctg} = $hit;
                } elsif (!($contigHIT2{$ctg}) && ($hit ne $contigHIT1{$ctg})){
                    $contigHIT2{$ctg} = $hit;
                }
            }
        }
        
        # not hits = over-purged
        if (!($contigHIT1{$ctg})){
            if ($over_purge_mode eq 'recheck'){
                $contigOVPURG{$ctg} = 1;
            } else {
                $contigREASSIGN{$ctg} = 0;
                msg("\tContig $ctg added back to primary assembly");
            }
        } else {
            my $assignment;
            
            # check if the contig is 100% repeat
            if ($contigMLEN{$ctg} == 0){
                $assignment = "r";
            } else {
                # get minimap2 alignments for ref1 and ref2
                ($contigHIT1{$ctg}, $contigHIT2{$ctg}, $alignments) = minimap2_alignments($contigHIT1{$ctg}, $contigHIT2{$ctg}, $ctg, 'recheck');
                
                # get 'maxmatch' and 'bestmatch' coverages from the sorted general output of minimap2, subtract masking regions on the fly
                my ($bestmatch, $maxmatch) = get_bm_mm(\@$alignments, $ctg);
                
                # guess the assignment
                $assignment = guess_assignment($bestmatch, $maxmatch);
            }
            
            # dotplots
            if ($dotplots){
                dotplots($contigHIT1{$ctg},$contigHIT2{$ctg},$ctg);
            }
            
            # check if over-purged
            if ($assignment !~ /[rh]/i){
                if ($over_purge_mode eq 'recheck'){
                    $contigOVPURG{$ctg} = 1;
                } else {
                    $contigREASSIGN{$ctg} = 0;
                    msg("\tContig $ctg added back to primary assembly");
                }
            } else {
                if ($over_purge_mode eq 'fix') {
                    # confirmed haplotig
                    $contigREASSIGN{$ctg} = 1;
                }
            }
        }
    }
    
    $available_threads->up(1);
    
    return;
}



sub dotplots {
    my $ref1 = $_[0];
    my $ref2 = $_[1];
    my $contig = $_[2];
    my $cmd;
    
    # dump the coverage
    my $tmpCov = "$COV_DIR/$contig.logcov";
    $tmpCov =~ s/\|.+//;
    open my $COV, '>', $tmpCov or err("Failed to open $tmpCov for writing");
    print $COV get_cov($contig);
    close $COV;
    
    # dotplot command 
    if (($ref2) && ($ref1)){
        $cmd = "$xvfb $RealBin/../scripts/dot_plot.Rscript '$UNASSIGNED/$contig.png' '$contig' $contigLEN{$contig} '$tmpCov' '$ref1' '$TMP_ALN/$contig.1.rdotplot' $contigLEN{$ref1} '$ref2' '$TMP_ALN/$contig.2.rdotplot' $contigLEN{$ref2} 1> '$TMP_ALN/$contig.Rscript.stderr' 2>&1\n";
    } elsif (($ref1) && !($ref2)) {
        $cmd = "$xvfb $RealBin/../scripts/dot_plot.Rscript '$UNASSIGNED/$contig.png' '$contig' $contigLEN{$contig} '$tmpCov' '$ref1' '$TMP_ALN/$contig.1.rdotplot' $contigLEN{$ref1} 1> '$TMP_ALN/$contig.Rscript.stderr' 2>&1\n";
    } elsif (($ref2) && !($ref1)){
        $cmd = "$xvfb $RealBin/../scripts/dot_plot.Rscript '$UNASSIGNED/$contig.png' '$contig' $contigLEN{$contig} '$tmpCov' '$ref2' '$TMP_ALN/$contig.2.rdotplot' $contigLEN{$ref2} 1> '$TMP_ALN/$contig.Rscript.stderr' 2>&1\n";
    }
    
    # make the dotplot
    if ($cmd){
        runcmd({ command => $cmd, logfile => "$TMP_ALN/$contig.Rscript.stderr", silent => 1 });
    } else {
        err("Contig $contig returned alignment score but no rdotplot files");
    }
    
    # cleanup
    for my $file ($tmpCov, "$TMP_ALN/$contig.1.rdotplot", "$TMP_ALN/$contig.2.rdotplot", "$TMP_ALN/$contig.Rscript.stderr"){
        if (-e $file){
            unlink $file;
        }
    }
    
    return;
}



#---END ITERATIVE STEP---



sub write_contig_associations {
    msg('Writing contig associations');
    
    my %primaries;
    my $p = 0;
    
    # iterate contigs, get all haplotigs for each primary, make new names for the primary contigs
    for my $ctg (sort { $contigLEN{$b} <=> $contigLEN{$a} } keys %contigLEN){
        if ($contigREASSIGN{$ctg}){
            push @{$primaries{$contigHIT1{$ctg}}}, $ctg;
            # rename assignments for the table output
            if ($contigASSIGN{$ctg} eq "h"){
                $contigASSIGN{$ctg} = "HAPLOTIG";
            } elsif ($contigASSIGN{$ctg} eq "r"){
                $contigASSIGN{$ctg} = "REPEAT";
            } else {
                err("unknown reassignment $contigASSIGN{$ctg}, this shouldn't have happened")
            }
            # move dotplot if exists
            if(-s "$UNASSIGNED/$ctg.png"){
                rename "$UNASSIGNED/$ctg.png", "$ASSIGNED/$ctg.png";
            }
        } else {
            $contigRENAME{$ctg} = ($falconNaming) ? sprintf("%06d", $p) . "F" : $ctg;
            $p++;
        }
    }
    
    # iterate primary contigs with haplotigs, make new names for the haplotigs, print associations to log file
    open $PTH, '>', $contig_paths or err("failed to open $contig_paths for writing");
    for my $ctg (sort { $contigLEN{$b} <=> $contigLEN{$a} } keys %primaries){
        print $PTH "$ctg,PRIMARY";
        my $h = 0;
        for my $htg (@{$primaries{$ctg}}){
            print $PTH " "x length("$ctg,PRIMARY") if ($h >= 1);
            print $PTH " -> $htg,", $contigASSIGN{$htg} , "\n";
            $contigRENAME{$htg} = ($falconNaming) ? $contigRENAME{$ctg} . "_" . sprintf("%03d", $h) : $htg;
            $h++;
        }
        print $PTH "\n";
    }
    close $PTH;
    return;
}



sub write_assembly {
    msg('Writing the reassignment table and new assembly files');
    
    # Init table of reassignments
    open my $CUT, '>', $out_reassignments or err("failed to open $out_reassignments for writing");
    print $CUT ($falconNaming) ? 
        "#reassigned_contig\ttop_hit_contig\tsecond_hit_contig\tbest_match_coverage\tmax_match_coverage\treassignment\tnew_name\n" : 
        "#reassigned_contig\ttop_hit_contig\tsecond_hit_contig\tbest_match_coverage\tmax_match_coverage\treassignment\n";
    
    # HAPLOTIGS
    open my $CUH, '>', $out_haplotigs or err("failed to open $out_haplotigs for writing");
    for my $ctg (sort { $contigRENAME{$a} cmp $contigRENAME{$b} } keys %contigLEN){
        if (($contigREASSIGN{$ctg}) && !($junk{$ctg})){
            my $c2 = $contigHIT2{$ctg} || "-";
            print $CUT ($falconNaming) ? 
                "$ctg\t$contigHIT1{$ctg}\t$c2\t$contigBM{$ctg}\t$contigMM{$ctg}\t$contigASSIGN{$ctg}\t$contigRENAME{$ctg}\n" :
                "$ctg\t$contigHIT1{$ctg}\t$c2\t$contigBM{$ctg}\t$contigMM{$ctg}\t$contigASSIGN{$ctg}\n";
            write_seq($CUH, $ctg);
        }
    }
    close $CUH;
    
    # PRIMARY CONTIGS
    print $CUT "#contigs_kept\n";
    open my $CUP, '>', $out_fasta or err("failed to open $out_fasta for writing");
    for my $ctg (sort { $contigLEN{$b} <=> $contigLEN{$a} } keys %contigLEN){
        if ( !($junk{$ctg}) && !($contigREASSIGN{$ctg}) ){
            my $c1 = $contigHIT1{$ctg} || "-";
            my $c2 = $contigHIT2{$ctg} || "-";
            my $bm = $contigBM{$ctg} || "-";
            my $mm = $contigMM{$ctg} || "-";
            print $CUT ($falconNaming) ?
                "$ctg\t$c1\t$c2\t$bm\t$mm\tKEEP\t$contigRENAME{$ctg}\n" :
                "$ctg\t$c1\t$c2\t$bm\t$mm\tKEEP\n";
            write_seq($CUP, $ctg);
        }
    }
    close $CUP;
    
    # ARTEFACTS
    open my $CUA, '>', $out_artefacts or err("failed to open $out_artefacts for writing");
    print $CUT "#junk_contigs\n";
    for my $ctg (sort { $contigLEN{$b} <=> $contigLEN{$a} } keys %contigLEN){
        if ($junk{$ctg}){
            print $CUT ($falconNaming) ?
                "$ctg\t-\t-\t-\t-\tJUNK\t$contigRENAME{$ctg}\n" :
                "$ctg\t-\t-\t-\t-\tJUNK\n";
            write_seq($CUA, $ctg);
        }
    }
    close $CUA;
    close $CUT;
    
    return;
}



sub write_seq {
    my $fh = $_[0];
    my $ctg = $_[1];
    
    # grab the contig seq
    my @seq = split/\n/, get_seq($ctg);
    
    # rename if renaming
    if($falconNaming){
        $seq[0] = ">$contigRENAME{$ctg}";
    }
    
    # print
    for my $l (@seq){
        print $fh "$l\n";
    }

    return;
}





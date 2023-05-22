# Scripts for processing repeat annotations

These scripts were developed for this study and have relatively narrow and study-specific scopes. Their overall purpose was to generate custom repeat library.

## LTRs2table.pl

This script processes multiple sources of annotations for LTR retrotransposons and reclassifies them according to the detected annotations. It is used to parse and annotate output from LTR_Harvest and LTR_Finder. It only outputs transposons with at least one hit against a LTR template or that remain "Unknown", i.e. repeats that get reclassified as something else than LTRs are discarded.

Usage example:

    ./LTRs2table.pl \
      -tag FIN \ # A tag, such as "FIN" as a short for "LTR_Finder"
      -lr LTRs.fa \ # The input LTR library
      -rc LTRs.fa.classified \ # Annotations detected using the RepeatModeler tool RepeatClassifier
      -te LTRs.fa.rexdb-metazoa.cls.lib \ # Annotations detected using a search with TEsorter against the RexDB/Metazoa database
      -psi LTRs.fa.allHits.chains.bestPerLocus.gff3 \ # Annotations detected using TransposonPSI
      -hmm LTRs.fa.dfam.hmm.domtbl \ # Annotations detected using HMMSEARCH against the Dfam database
          1> LTRs.fa.table.tsv \ # Tabular output with annotations and scores
          2> LTRs.fa.reclassified.fasta \ # Reclassified LTR repeats in FASTA format
          
It can read multiple files following "-psi" and "-hmm". It selects "winning" annotations according to following order:

1.  If a RepeatClassifier annotation is detected, use this annotation and assign the repeat the "ReC" label. Else, continue to the next step:
2.  If a TEsorter annotation is detected, use this annotation and assign the repeat the "TEs" label. Else, continue to the next step:
3.  If a TransposonPSI annotation is detected, use this annotation and assign the repeat the "PSI" label. Else, continue to the next step:
4.  If a HMMsearch annotation is detected, use this annotation and assign the repeat the "DFA" label.

The script outputs LTR classifications with headers that are compatible with RepeatMasker:

    ">repeat_LABEL#Order/Superfamily"
    
    Example:
    
    >seq_c_98391_5186_12351_FIN_ReC#LTR/Pao

Here "seq_c_98391" is the name of sequence in the reference genome, "5186_12351" specifies the start/stop of the transposable element, "FIN" indicates it was detected using "LTR_Finder", "ReC" that it was classified with RepeatClassifier and "LTR/Pao" is the biological classification of the repeat.

## subdivide_repeats.pl

This script subdivides LTRs according to the identity of the 5' and 3' LTR region, as specified in a GFF3 file produced by LTR_retriever. The user can specify one or more lower identify thresholds to partition the LTRs accordingly.

### Usage example:

    ./subdivide_repeats.pl \
	    LTRs_pass.gff3 \ # The GFF3 file
	    LTRs.fa.reclassified.fasta \ # The reclassified LTRs (see script above)
	    0.85 0.99 # One or more identity tresholds (0.85=85% or higher identity; 0.99=99% or higher identity)
  
## collect_hits_hmm.pl

This script parses the tabular output from a HMMSEARCH against the LTR GyDB database. It sub-classifies the LTRs into being "complete" or "incomplete" according to how many of the expected LTR GAG/Pol protein domains are detected.

In the preceeding search, the LTRs need to have been translated into amino acids in all six reading frames.

### Usage example:

		./collect_hits_hmm.pl \
			LTRs.6_frames.aa.fa.GyDB.hmm.tbl \ # Basename of the output file
			5 \ # The minimum number of domains to be considered complete
			LTRs.6_frames.aa.fa.GyDB.hmm.tbl # The HMMSEARCH table
      
The downstream script **subdivide_repeats_list.pl** can subdivide the repeats based a the list generated, such as that generated in the preceeding step.

### Usage example:

		./subdivide_repeats_list.pl \
			complete \ # Repeats in the list will go in a file tagged with this label
			incomplete \ # Repeats not in the list will go in a file tagged with this label
			LTRs.6_frames.aa.fa.GyDB.hmm.tbl.hits.complete.csv \ # A repeat list
			LTRs.fa.reclassified.fasta # The input repeat FASTA file
      
 The downstream script **tag_seq.pl** adds tags preceeding the "#" in the repeat header, such as "C" for complete or "I" for incomplete.
 
    ./tag_seq.pl "C" LTRs.fa.reclassified.fasta.id_0.85.fasta.complete.fasta > LTRs.fa.reclassified.fasta.id_0.85.fasta.complete.tagged.fasta

## check_masking_rate.pl

Taking the output from RepeatMasker, this script removes repeats from a FASTA file that are masked by >= 80% (hard-coded) by other repeats in a repeat library. A list of repeats to keep is generated as output.

### Usage example:

    ./check_masking_rate.pl \
        LTRs.fa.reclassified.fasta.id_0.85.fasta.complete.fasta.masked \
        > LTRs.fa.reclassified.fasta.id_0.85.fasta.complete.fasta.masked.keep.list.csv

## DNA2table.pl, LINE2table.pl, RC2table.pl, LTR2table.pl, RM2table.pl

These scripts annotate various groups of transposons. They assume that the primary output from the search tool (e.g. TransposonPSI) has already been subdivided according to the major classes/orders of repeats. The annotation priority follow the same order as above.

    DNA => DNA transposons
    LINE => LINEs
    RC => Rolling circles
    LTR => LTR retrotransposons

## DNA2table.pl

This script processes multiple sources of annotations for DNA transposons and reclassifies them according to the detected annotations. It is used to parse and annotate output from TransposonPSI and dnaPipeTE. It works similarly to **LTRs2table.pl** above but has a simpler input and output interface. It only outputs transposons with at least one hit against a DNA transposon template.

Usage example:

    ./DNA2table.pl \
        PSI \ # A tag, such as "PSI" as a short for "TransposonPSI"
        DNA.fa.non_copia_gypsy.fasta \ # The input transposon library
        DNA.fa.non_copia_gypsy.fasta.classified \ # Annotations detected using the RepeatModeler tool RepeatClassifier
        DNA.fa.non_copia_gypsy.fasta.renamed.rexdb-tir.cls.lib \ # Annotations detected using a search with TEsorter against the RexDB/TIR database
          1> DNA.fa.non_copia_gypsy.fasta.table.tsv # Tabular output with annotations and scores
          2> DNA.fa.non_copia_gypsy.fasta.reclassified # Reclassified transposons in FASTA format
  
## LINE2table.pl

This script processes multiple sources of annotations for LINEs and reclassifies them according to the detected annotations. It is used to parse and annotate output from TransposonPSI and dnaPipeTE. It works similarly to **LTRs2table.pl** above but has a simpler input and output interface. It only outputs transposons with at least one hit against a LINE template.

Usage example:
 
    ./LINE2table.pl \
        PSI \ # A tag, such as "PSI" as a short for "TransposonPSI"
        LINE.fa.non_copia_gypsy.fasta \ # The input transposon library
        LINE.fa.non_copia_gypsy.fasta.classified \ # Annotations detected using the RepeatModeler tool RepeatClassifier
        LINE.fa.non_copia_gypsy.fasta.clustered.renamed.rexdb-metazoa.cls.lib \ # Annotations detected using a search with TEsorter against the RexDB/Metazoa database
          1> LINE.fa.non_copia_gypsy.fasta.table.tsv \ # Tabular output with annotations and scores
          2> LINE.fa.non_copia_gypsy.fasta.reclassified # Reclassified transposons in FASTA format
          
## RC2table.pl

This script processes multiple sources of annotations for Rolling Circle helitrons and reclassifies them according to the detected annotations. It is used to parse and annotate output from TransposonPSI and dnaPipeTE. It works similarly to **LTRs2table.pl** above but has a simpler input and output interface. It only outputs transposons with at least one hit against a helitron template.

    ./RC2table.pl \
        PSI \ # A tag, such as "PSI" as a short for "TransposonPSI"
        RC.fa.non_copia_gypsy.fasta \ # The input transposon library
        RC.fa.non_copia_gypsy.fasta.classified \ # Annotations detected using the RepeatModeler tool RepeatClassifier
        RC.fa.non_copia_gypsy.fasta.clustered.renamed.rexdb-metazoa.cls.lib \ # Annotations detected using a search with TEsorter against the RexDB/Metazoa database
          1> RC.fa.non_copia_gypsy.fasta.table.tsv \ # Tabular output with annotations and scores
          2> RC.fa.non_copia_gypsy.fasta.reclassified # Reclassified transposons in FASTA format

## LTR2table.pl

This script processes multiple sources of annotations for LTRs and reclassifies them according to the detected annotations. It is used to parse and annotate output from TransposonPSI and dnaPipeTE. It works similarly to **LTRs2table.pl** above but has a simpler input and output interface. It only outputs transposons with at least one hit against a LTR template.

    ./LTRtable.pl \
        PSI \ # A tag, such as "PSI" as a short for "TransposonPSI"
        LTR.fa.non_copia_gypsy.fasta \ # The input transposon library
        LTR.fa.non_copia_gypsy.fasta.classified \ # Annotations detected using the RepeatModeler tool RepeatClassifier
        LTR.fa.non_copia_gypsy.fasta.clustered.renamed.rexdb-metazoa.cls.lib \ # Annotations detected using a search with TEsorter against the RexDB/Metazoa database
          1> LTR.fa.non_copia_gypsy.fasta.table.tsv \ # Tabular output with annotations and scores
          2> LTR.fa.non_copia_gypsy.fasta.reclassified # Reclassified transposons in FASTA format

## RM2table.pl

This script processes multiple sources of repeat annotations and reclassifies them according to the detected annotations. It is used to parse and annotate output from RepeatModeler. It works similarly to **LTRs2table.pl** above.

	./RM2table.pl \
      -tag RM \
      -rm repeatmodeler.fa.orig \
		  -rc repeatmodeler.fa.classified \ # Annotations detected using the RepeatModeler tool RepeatClassifier
		  -te repeatmodeler.fa.renamed.rexdb-metazoa.cls.lib \ # Annotations detected using a search with TEsorter against the RexDB/Metazoa database
		  -psi repeatmodeler.fa.classified.allHits.chains.bestPerLocus.gff3 \ # Annotations using TransposonPSI \
		  -hmm LTRs.fa.dfam.hmm.domtbl \ # Annotations detected using HMMSEARCH against SINEs in the Dfam database
		    1> repeatmodeler.fa.orig.table.csv \
		    2> repeatmodeler.fa.orig.reclassified

## dnapipete2repeatmasker_library.pl

This script filters and selects putative *de novo* assembled repeats from the dnaPipeTE output:

	dnapipete2repeatmasker_library.pl \
		0.1 \ minimum proportion of the contig annotated as a repeat in the heat
		0.1 \ minimum proportion of the repeat template in the hit
		one_RM_hit_per_Trinity_contigs \ # Tabular input based on BLAST hits
		annoted.fasta # basename of output file
  

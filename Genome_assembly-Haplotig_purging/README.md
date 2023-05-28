Scripts for processing assembly haplotigs

The script provided below modifies the **purge.pl** script originally written by Michael Roach (Australian Wine Research Institute) and part of the open source Purge Haplotigs tool v1.1.0. That tool can be found here:
https://bitbucket.org/mroachawri/purge_haplotigs/wiki/Updates

The original purge.pl script has a permissive license:

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

This script is executed when running the "purge_haplotigs purge" subcommand.

The modified version of the script is derived from Michael Roach's code and re-distributed under the same license. The modified version of the script was implemented to increase run-time performance in processing the very large and repeated genome by:

- using minimap2 without making an mmi index file on disk
- using faster parsing of tabular text files in many places: instead of splitting lines, it performs regular expression pattern matching to extract data points from specific fields
- using string-based instead of array based data structures, reducing memory usage. This pertains to storing associations between different contigs.
- using a faster internal handling of repeats instead of repeatedly calling and executing the external tool "bedtools".

To accomplish the latter, it first calculates the length of each query contig (that is also a suspect) after removing repeat sequences. Target contigs that are long enough to have the potential to produce an alignment of significant length against that corrected query length are then considered in-depth. The script creates a mask for the full length of the query contig that is specific for this target contig, setting all bases to the symbol "0". Aligned regions reported in the PAF file (alignment output from minimap2) are then rewritten using the state "1". Repeated regions are then written on top the mask using the state "2". The remaining "1"s are counted and if they surpass the treshhold, for example 70%, of the corrected query length the query and target are kept as potential haplotigs. In this version of the script, **purge.pl** parses the original PAF twice: in the first pass, the significant alignments are detected and in the second pass, they are printed to a simplified PAF that only contains these. This PAF is then indexed and used for the downstream analyses. A "hit summary" file is written to disk but not read back as the necessary information is already kept in memory. The script then proceeds with the original iterative algorithm.


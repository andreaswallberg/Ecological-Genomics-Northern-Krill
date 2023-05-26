# Ecological-Genomics-Northern-Krill

This repository contais code used to perform genome assembly, annotation and population genomics analyses of krill.

Some of these scripts and tools are study-specific (e.g. processing of repeats), whereas others have general-purpose use (e.g. processing of FASTA files, genome-masking and population genetics tools). Hints about this are given in the description for each respective script.

Unless stated otherwise, all scripts assume that the tab character separate fields in tabular text files (even when the file ending says "csv").

All code is open source. Unless stated otherwise for specific scripts, it is licensed under the same terms as Perl, which means it is dually-licensed under either the Artistic or GPL licenses.

One script incorporates a modified version of the BioPerl method used to compute Tajima's D. That original code was written by Jason Stajich and is available here:
https://metacpan.org/release/CJFIELDS/BioPerl-1.6.924/source/Bio/PopGen/Statistics.pm

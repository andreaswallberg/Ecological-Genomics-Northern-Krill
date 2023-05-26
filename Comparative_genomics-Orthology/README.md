# Scripts for processing orthology data

These scripts processes input or output from comparative analyses of orthology assessments and gene family evolution

## process_groups_all.pl

This scripts takes the tabular gene family output table from SwiftOrtho and reformats it into a format that is compatible with CAFE.

The output format from SwiftOrtho and input format for CAFE are described here, respectively:
https://github.com/Rinoahu/SwiftOrtho
https://github.com/hahnlab/CAFE5

It is specific to this study and contains hard-coded species labels in the script.

### Usage example:

  ./process_groups_all.pl all_matches.out.30_30.orth.apc > all_matches.out.30_30.orth.apc.table.ALL.csv

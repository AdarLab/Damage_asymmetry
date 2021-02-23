#!/bin/bash
dirName=$1
faName=$2
#################################################################################################
# This script gets exons/introns fasta file and merge the sequences that related to the same gene.
#################################################################################################
awk 'NR%2{printf "%s ",$0;next;}1' $dirName/sequence/$faName | sed 's/::.*)//' | awk 'NF>1{a[$1] = a[$1]""$2};END{for(i in a)print i"\n"a[i]}' > $dirName/sequence/continuous_$faName
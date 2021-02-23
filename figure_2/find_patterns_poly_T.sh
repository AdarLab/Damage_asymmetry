#!/bin/bash
#fuzznuc version EMBOSS:6.6.0.0

###########################################################################
# This script finds all occurances of poly-T sequences in the human genome.
###########################################################################


path=your_path
refGenome=hg38.fa
mkdir -p $path/patterns/wholeGenome

declare -a patterns=("T" "TT" "TTT" "TTTT" "TTTTT" "TTTTTT" "TTTTTTT" "TTTTTTTT" "TTTTTTTTT" "TTTTTTTTTT") # T(1) to T(10).
# Scan all the pattern and use fuzznuc tool for extracting all positions of the oligonucleotide in the human genome:
# {T} - avoid patterns with T upstream and downstream (e.g. TTT wonn't be considered as TT).
# -complement true - finds the occurences on the opposite strand either.
for pattern in "${patterns[@]}"
do
  gffName="$pattern"".gff"
  bedName="$pattern"".bed"
  fuzznuc -sequence $refGenome -pattern {T}"$pattern"{T} -complement true -outfile $path/patterns/wholeGenome/$gffName -rformat gff 
  wait
  #Arrange the files
  sed '/^#/ d' $path/patterns/wholeGenome/$gffName | awk 'BEGIN{OFS = "\t"} {if(length($1) < 6 ){print($1,$4,$5-1,0,$6,$7)}}' > $path/patterns/wholeGenome/$bedName
done
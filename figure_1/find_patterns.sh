#!/bin/bash
# fuzznuc version EMBOSS:6.6.0.0
# bedtools v2.27.1
######################################################################################################################
# This script gets bed file of genomic element and finds all occurances of a specific oligonucleotide in this element.
######################################################################################################################

path=your_data_path
refGenome=/hg38.fa
bedFiles=$path/annotation
mkdir -p $path/patterns/wholeGenome
pattern="TT"

gffName="$pattern"".gff"
bedName="$pattern"".bed"
# Use fuzznuc tool for extracting all positions of the oligonucleotide in the human genome, "-complement true" - finds the occurences on the opposite strand either.
fuzznuc -sequence $refGenome -pattern $pattern -complement true -outfile $path/patterns/wholeGenome/$gffName -rformat gff" 
wait
sed '/^#/ d' $path/patterns/wholeGenome/$gffName | awk 'BEGIN{OFS = "\t"} {if(length($1) < 6 ){print($1,$4-1,$5,0,$6,$7)}}' > $path/patterns/wholeGenome/$bedName
wait
mkdir -p $path/patterns/transcribedStrand 
mkdir -p $path/patterns/nonTranscribedStrand 

for elementFile in $bedFiles/*.bed
do
  for patternFile in $path/patterns/wholeGenome/*.bed
  do
    elementName=`basename $elementFile`
    patternName=`basename $patternFile`
    outputName="${elementName%.bed*}""_""$patternName"
    #1) "-s" force strandedness. That is, only report hits that overlap on the same strand.
    #2) "-S" require different strandedness. That is, only report hits that overlap on the opposite strand.
    bedtools intersect -a $patternFile -b $elementFile -wa -s > $path/patterns/nonTranscribedStrand/$outputName
    bedtools intersect -a $patternFile -b $elementFile -wa -S > $path/patterns/transcribedStrand/$outputName
  done
done
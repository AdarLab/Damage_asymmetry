#!/bin/bash
# fuzznuc version EMBOSS:6.6.0.0
# bedtools v2.27.1
#################################################################################################################
# This script gets bed file of genomic element and finds all occurances of poly T coordinates within this element
# using bedtools slop and intersect commands for plotting average profiles of poly Ts.
#################################################################################################################

genome=hg38.chrom.sizes
currentPath=hg38_poly_T
elementsPath=annotation
patternPath=$currentPath/patterns/wholeGenome

mkdir -p $currentPath/patterns/transcribedStrand 
mkdir -p $currentPath/patterns/nonTranscribedStrand    

for patternFile in $currentPath/patterns/wholeGenome/*.bed
do
  for elementFile in $elementsPath/*.bed
  do
    patternName=`basename $patternFile`
    elementName=`basename $elementFile`
    outputName="${elementName%.bed*}""_""$patternName"
    #1) "-s" force strandedness. That is, only report hits that overlap on the same strand.
    #2) "-S" require different strandedness. That is, only report hits that overlap on the opposite strand.
    bedtools slop -i $elementFile -g $genome -b 10 | bedtools intersect -a $patternFile -b - -wa -s > $currentPath/patterns/nonTranscribedStrand/$outputName
    bedtools slop -i $elementFile -g $genome -b 10 | bedtools intersect -a $patternFile -b - -wa -S > $currentPath/patterns/transcribedStrand/$outputName
  done
done   
    

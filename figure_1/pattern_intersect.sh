#!/bin/bash
# bedtools version v2.26.0
#######################################################################################################################
# This script gets bed files containing TT coordinates in the entire genome,
# bed files of humen genomic elements and extract the damages the interseceted the given elements for further analyses.
#######################################################################################################################
genome=hg38.chrom.sizes
path=your_path
profilesPath=$path/averageProfilesData
bedFiles=$profilesPath/annotation



declare -a patterns=("TT")

mkdir -p $profilesPath/patterns/transcribedStrand 
mkdir -p $profilesPath/patterns/nonTranscribedStrand 

# Scan the elements' bed files and expand the coordinates using bedtools slop in order to avoid partial overlapped TTs.
# Use bedtools intersect with -wa option to write the original Damage-seq reads for each overlap.
for elementFile in $bedFiles/*.bed
do
  for patternFile in $path/patterns/wholeGenome/*.bed
  do
    elementName=`basename $elementFile`
    patternName=`basename $patternFile`
    outputName="${elementName%.bed*}""_""$patternName"
    #1) "-s" force strandedness. That is, only report hits that overlap on the same strand.
    #2) "-S" require different strandedness. That is, only report hits that overlap on the opposite strand.
    bedtools slop -i $elementFile -g $genome -b 10 | bedtools intersect -a $patternFile -b - -wa -s > $profilesPath/patterns/nonTranscribedStrand/$outputName
    bedtools slop -i $elementFile -g $genome -b 10 | bedtools intersect -a $patternFile -b - -wa -S > $profilesPath/patterns/transcribedStrand/$outputName
  done
done
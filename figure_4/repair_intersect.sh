#!/bin/bash
# bedtools version v2.26.0
##############################################################################################
# This script gets XR-seq bed files of 4 nucleotides, bed files of humen genomic elements
# and extract the damages the interseceted the given elements for further analyses.
##############################################################################################
genome=/cs/icore/elisheva.h/hg38-ref_genome/hg38.chrom.sizes
dir=$1
readFile=$2
elementFile=$3
outputName=$4
intervals=$5
echo "$dir"
echo "$readFile"
echo "$elementFile"
#1) "-s" force strandedness. That is, only report hits that overlap on the same strand.
#2) "-S" require different strandedness. That is, only report hits that overlap on the opposite strand.

# If the intervals are for plotting averge profiles extend the intervals by 10 bases on each side.
if [ $intervals == "plot" ]
then
  bedtools slop -i $elementFile -g $genome -b 10 | bedtools intersect -a $readFile -b - -wa -s > $dir/nonTranscribedStrand/$outputName
  bedtools slop -i $elementFile -g $genome -b 10 | bedtools intersect -a $readFile -b - -wa -S > $dir/transcribedStrand/$outputName
else
  bedtools intersect -a $readFile -b $elementFile -wa -s > $dir/nonTranscribedStrand/$outputName
  bedtools intersect -a $readFile -b $elementFile -wa -S > $dir/transcribedStrand/$outputName
fi


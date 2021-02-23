#!/bin/bash
# Bedtools version: v2.26.0

############################################################################################
# This script gets C->T mutation coordinates in Melanoma and usees Bedtools coverage command
# for evaluting the mutation counts on each strand of human enhancers.
############################################################################################

dataDir=your_path
enhancersFile=$dataDir/stranded_enhancers/hg19_enhancers_500_sorted.bed
mutationsFile=$1
mkdir -p $dataDir/mutationsCoverage
coverageDir=$dataDir/mutationsCoverage

# Report how many C->T mutations are in each strand of the enhancers and normalize the counts to the element length.
mkdir -p $coverageDir/transcribedStrand
mkdir -p $coverageDir/nonTranscribedStrand 
mutationsFileName=`basename $mutationsFile`
echo "$mutationsFileName"
function mutationsCoverage
{
	enhancersFile=$1
	enhancersFileName=`basename $enhancersFile _500_sorted.bed`
	bedtools coverage -s -a $enhancersFile -b $mutationsFile | awk '{print($0=$0"\t"($7/($3-$2)))}' > $coverageDir/nonTranscribedStrand/"$enhancersFileName"_"$mutationsFileName" #-s: select only the reads on the same strand of the original element.
	bedtools coverage -S -a $enhancersFile -b $mutationsFile | awk '{print($0=$0"\t"($7/($3-$2)))}' > $coverageDir/transcribedStrand/"$enhancersFileName"_"$mutationsFileName" #-s: select only the reads on the opposite strand.
}
mutationsCoverage $dataDir/stranded_enhancers/hg19_enhancers_500_sorted.bed


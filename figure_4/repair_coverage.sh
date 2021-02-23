#!/bin/bash
#Bedtools version: v2.26.0
#################################################################################################
# This script gets XR-seq bed files of 4 nucleotides, bed files of humen genomic elements and
# count the number of appearances of the 4 nucleotides in every element.
#################################################################################################
dir=$1
readFile=$2
elementFile=$3
echo "$dir"
echo "$readFile"
echo "$elementFile"
readName=`basename $readFile`
outputName="${readName%.bed*}""_cov.bed"
bedtools coverage -a $elementFile -b $readFile > $dir/$outputName

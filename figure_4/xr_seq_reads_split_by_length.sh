#!/bin/bash
# bedtools v2.26.0

#####################################################################################################################################
# This script gets the XR-seq filtered reads, splits them into sapareated files acording to their length and extract their sequences.
#####################################################################################################################################

dataPath=XR_Seq_data
refGenome=hg38.fa
for bedFile in $dataPath/bedFiles/*.bed
do
  fileName=`basename $bedFile _filtered_unique_sorted.bed` # Extract the "pure" file name.
  mkdir -p $dataPath/length/bedFiles/"$fileName"_length
  # Enter into the output directory.
  cd $dataPath/length/bedFiles/"$fileName"_length
  awk 'BEGIN {FS = OFS = "\t"}; ($3-$2 > 14 && $3-$2 < 36 && length($1) < 6){print > (""$3-$2".bed")}' $bedFile
done
wait

# Scan all the directories of the lengths
for dir in $dataPath/length/bedFiles/*
do
  dirName=`basename $dir` # Extract the "pure" file name.
  mkdir -p $dataPath/length/fastaFiles/$dirName
  # Scan the bed files and extract the sequences of the intervals.
  for bedFile in $dir/*.bed
  do
    fileName=`basename $bedFile .bed`
    faName="$fileName".fa
    bedtools getfasta -fi $refGenome -bed $bedFile -s -name > $dataPath/length/fastaFiles/$dirName/$faName
  done
  
done
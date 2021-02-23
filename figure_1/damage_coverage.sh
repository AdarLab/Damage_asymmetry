#!/bin/bash
#Bedtools version: v2.26.0
#################################################################################################
# This script gets Damage-seq bed files of 4 nucleotides, bed files of humen genomic elements and
# count the number of appearances of the 4 nucleotides in every element.
#################################################################################################
dataDir=your_data_path
human_elements_path=$dataDir/annotation
damage_seq_bedFiles_path=$dataDir/filteredReads/bedFiles

mkdir -p $dataDir/damageCoverage
coverageDir=$dataDir/damageCoverage

declare -a elements=("uniq_coding_genes" "chopped_exons" "chopped_introns") # All the elements.
declare -a strands=("transcribed" "nonTranscribed") # None_transcribed and transcribed strand respectively.

# Scan all the genomic elements, count how many damages are in each element and normalize the counts to the element length.
for element in "${elements[@]}"
do
  elementFile=$human_elements_path/*"$element"*
  mkdir -p $coverageDir/transcribed
  mkdir -p $coverageDir/nonTranscribed 
  for damageFile in $damage_seq_bedFiles_path/*.bed
  do
    damageFileName=`basename $damageFile`
    echo "$damageFileName"
    bedtools coverage -s -a $elementFile -b $damageFile | awk '{print($0=$0"\t"($7/($3-$2)))}' > $dataDir/damageCoverage/nonTranscribed/$element"_"$damageFileName #-s: select only the reads on the same strand of the original element (transcribed strand).
    bedtools coverage -S -a $elementFile -b $damageFile | awk '{print($0=$0"\t"($7/($3-$2)))}' > $dataDir/damageCoverage/transcribed/$element"_"$damageFileName #-S: select only the reads on the opposite strand (coding strand).
  done    
done

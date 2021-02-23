#!/bin/bash
# bedtools version v2.26.0
##############################################################################################
# This script gets Damage-seq bed files of 4 nucleotides, bed files of humen genomic elements
# and extract the damages the interseceted the given elements for further analyses.
##############################################################################################

genome=hg38.chrom.sizes
dataDir=your_data_path
human_elements_path=$dataDir/annotation
damage_seq_bedFiles_path=$dataDir/filteredReads/bedFiles

mkdir -p $dataDir/damageIntersect
intersectDir=$dataDir/damageIntersect


###################################
# Intersection for average profiles
###################################
plotDir=$dataDir/averageProfilesData
mkdir -p $plotDir/damageIntersect
intersectDirPlot=$plotDir/damageIntersect
# Scan the elements' bed files and expand the coordinates using bedtools slop in order to avoid partial overlapped damages.
# Use bedtools intersect with -wa option to write the original Damage-seq reads for each overlap.
for elementFile in $human_elements_path/*.bed
do
  elementName=`basename $elementFile .bed`
  mkdir -p $intersectDirPlot/transcribed
  mkdir -p $intersectDirPlot/nonTranscribed
  for damageFile in $damage_seq_bedFiles_path/*.bed
  do
    damageFileName=`basename $damageFile`
    echo "$damageFileName"
    bedtools slop -i $elementFile -g $genome -b 10 | bedtools intersect -a $damageFile -b - -wa -S  > $intersectDirPlot/transcribed/$elementName"_"$damageFileName #-s: select only the reads on the opposite (transcribed) strand.
    bedtools slop -i $elementFile -g $genome -b 10 | bedtools intersect -a $damageFile -b - -wa -s > $intersectDirPlot/nonTranscribed/$elementName"_"$damageFileName #-s: select only the reads on the same (coding) strand of the original element.    
  done    
done
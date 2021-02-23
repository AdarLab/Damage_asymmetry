#!/bin/bash
# bedtools version v2.26.0
##################################################################################################
# This script gets Damage-seq bed files of 4 nucleotides, TT coordinates in the entire genome,
# XR-seq bed files of 4 nucleotides, bed files of humen genomic elements and extract the elements
# that interseceted with human enhancers for further analyses.
##################################################################################################
dataDir=your_data_path
damage_seq_bedFiles_path=Damage_seq_path/bedFiles
xr_seq_bedFiles_path=XR_seq_path/4_nt/bedFiles
pattern_path=patterns/wholeGenome
hg38_genome=hg38.chrom.sizes

mkdir -p $dataDir/intersectedData/damage
mkdir -p $dataDir/intersectedData/TT

intersectDir=$dataDir/intersectedData

# Scan the enhancers' bed files and send them to the intersect_enhancers script.

for bedFile in $dataDir/stranded_enhancers/*.bed
  do
  dir=`dirname $bedFile`
  dirName=`basename $dir`
  mkdir -p $intersectDir/damage/minus/$version
  mkdir -p $intersectDir/damage/plus/$version 
      
  mkdir -p $intersectDir/TT/minus/$version
  mkdir -p $intersectDir/TT/plus/$version 
      
  mkdir -p $intersectDir/repair/minus/$version
  mkdir -p $intersectDir/repair/plus/$version 
  
  # Damage.
  sbatch intersect_enhancers.sh $bedFile $damage_seq_bedFiles_path/NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed $intersectDir/damage/minus/"$version" $intersectDir/damage/plus/"$version" $hg38_genome
  
  # TT.
  sbatch intersect_enhancers.sh $bedFile $pattern_path/TT.bed $intersectDir/TT/minus/"$version" $intersectDir/TT/plus/"$version" $hg38_genome
  
  # Repair
  for repairFile in $xr_seq_bedFiles_path/*.bed
  do
    sbatch intersect_enhancers.sh $bedFile $repairFile $intersectDir/repair/minus/"$version" $intersectDir/repair/plus/"$version" $hg38_genome
  done
done 

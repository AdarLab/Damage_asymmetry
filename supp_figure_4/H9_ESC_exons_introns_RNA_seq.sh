#!/bin/bash

dataDir=RNA_seq
exonsFile=chopped_exons_seq.fa
intronsFile=chopped_introns_seq.fa
##########################################################################################################################################################
# This script gets human exon and introns FASTA files, merge the elements which are related to same transcript and genearates dinucleotide counting table.
##########################################################################################################################################################
mkdir -p $dataDir/sequence 
mkdir -p $dataDir/tables/counting
mkdir -p $dataDir/tables/ratio

function merge_elements
{  
  faFile=$1
  element=$2
  faName=`basename $faFile`

  awk -F'_'$element'' '{print $1}' $faFile > $dataDir/sequence/$faName
  wait
  merge_elements.sh $dataDir $faName
  wait
  python tables.py $dataDir/tables $dataDir/sequence/continuous_$faName continuous_$faName" 

}

merge_elements $exonsFile exon
merge_elements $intronsFile intron

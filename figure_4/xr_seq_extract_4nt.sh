#!/bin/bash

###############################################################################
# This script gets XR-seq read bed files and exract the positions of 
# the 4 nt sourranding the dimer based on the nucleotide frequencies bar plots. 
###############################################################################

dataPath=XR_Seq_data
lengthPath=length
mkdir -p $dataPath/4_nt/bedFiles
declare -a lengths=(23 24 25 26 27 28 29 30) #All the lengths.
declare -a NHF1_CPD_poses=(16 19 17 20 18 21 18 21 19 22 20 23 21 24 22 25) #The start and end of the 4 nucleotides. Each following values represent the start and end of the 4 nt.
declare -a CSB_CPD_poses=(17 20 17 20 18 21 18 21 19 22 20 23 21 24 21 24)
declare -a XPC_CPD_poses=(16 19 17 20 17 20 18 21 19 22 20 23 21 24 21 24)

function extract_coordinates()
{
  dirName=$1
  outputPath=$2
  outputName=$3
  posesName=$4[@] # Get the array
  poses=("${!posesName}") # Get the array itself.
  
  mkdir -p $outputPath/$outputName
  for i in "${!lengths[@]}"
  do 
  bedFile=$lengthPath/"$dirName"/"${lengths[$i]}".bed
  fileName=`basename $bedFile`
  pose_1=$((2*$i)) # The position of the first nuc out of the 4.
  pose_4=$((2*$i+1)) # The position of the last nuc.
  # Extract the coordinates of the 4 nucleotides from the bed file based on the given positions of the positions array.
  awk -v pos1=${poses[$pose_1]} -v pos4=${poses[$pose_4]} '{OFS="\t"}; {if($6 == "+"){print($1,$2+pos1-1,$2+pos4,$4,$1":"$2"-"$3,$6)} else {print($1,$3-pos4, $3-pos1+1,$4,$1":"$2"-"$3,$6)}}' $bedFile > $outputPath/$outputName/$fileName
  done
  wait
  cat $outputPath/$outputName/*.bed > $outputPath/"$outputName"_4nt.bed
}

extract_coordinates "NHF1_CPD_length" "$dataPath/4_nt/bedFiles" "NHF1_CPD" NHF1_CPD_poses
extract_coordinates "CSB_CPD_length" "$dataPath/4_nt/bedFiles" "CSB_CPD" CSB_CPD_poses
extract_coordinates "XPC_CPD_length" "$dataPath/4_nt/bedFiles" "XPC_CPD" XPC_CPD_poses
 

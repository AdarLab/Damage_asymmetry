#!/bin/bash
# bedtools version v2.26.0
#############################################################################################
# This script gets stranded enhancer coordinates, extract 150 bases relative to the center and
# determines the strand as - before the center and + after the center.
#############################################################################################

path=your_path
enhancersFile=$path/enhancers_sorted_uniq.bed 

function extract_coordinates
{
  bedFile=$1
  extension=$2
  outputDir=$3
  fileName=`basename $bedFile _sorted_uniq.bed`

  # Extend X nt from each side of the given midpoint coordinates (columns 7 & 8).
  awk -v var=$extension '{OFS="\t"}; {print($1, $7-var, $7, $4, $5, "'-'")}' $bedFile > $outputDir/"$fileName"_plot_"$extension".bed
  awk -v var=$extension '{OFS="\t"}; {print($1, $8, $8+var, $4, $5, "'+'")}' $bedFile >> $outputDir/"$fileName"_plot_"$extension".bed
  # Extract 140 nt on each side of the enhancer's center on both strands.
  awk -v var=$extension '{OFS="\t"}; {print($1, $7-var, $8+var, $4, $5, "'+'")}' $bedFile > $outputDir/"$fileName"_plot_plus_"$extension".bed
  awk -v var=$extension '{OFS="\t"}; {print($1, $7-var, $8+var, $4, $5, "'-'")}' $bedFile > $outputDir/"$fileName"_plot_minus_"$extension".bed
  
  wait
  bedtools sort -i $outputDir/"$fileName"_plot_"$extension".bed > $outputDir/"$fileName"_plot_"$extension"_sorted.bed
  bedtools sort -i $outputDir/"$fileName"_plot_plus_"$extension".bed > $outputDir/"$fileName"_plot_plus_"$extension"_sorted.bed
  bedtools sort -i $outputDir/"$fileName"_plot_minus_"$extension".bed > $outputDir/"$fileName"_plot_minus_"$extension"_sorted.bed
  wait 
  rm -f $outputDir/"$fileName"_plot_"$extension".bed $outputDir/"$fileName"_plot_plus_"$extension".bed $outputDir/"$fileName"_plot_minus_"$extension".bed # Remove unneccesary files  
}

mkdir -p $path/stranded_enhancers/given_midpoint

extract_coordinates $bedFile 500 $path/stranded_enhancers/given_midpoint # Expand 500 bp for the box plots analysis.

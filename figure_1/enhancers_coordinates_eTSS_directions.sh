#!/bin/bash
# bedtools version v2.26.0
#######################################################################################################################################
# This script gets stranded enhancer coordinates, extract 750 bases upstream and down stream the eTSS for the averge profiles analysis.
#######################################################################################################################################

path=your_path
genesFile_hg38=$path/hg38_RefSeq_genes_17_12_19.bed

function filter_enhancers
{
  elementFile=$1
  genesFile=$2
  outputDir=$3

  fileName=`basename $elementFile .bed`
  # Keep the longest one of the same 'TSS'(using the sort commands).
  cat $elementFile | sort -rn -k5,5 | sort -u -k1,1 -k2,2 |  sortBed -i - > $outputDir/sorted_uniqTSS.bed  
  
  # Remove enhancers that have other enhancers or coding genes at a distance of less than 2 Kb from each side.
  bedtools sort -i $elementFile | bedtools closest -a $outputDir/sorted_uniqTSS.bed -b - -D ref -N | awk '{OFS="\t"}; ($25 > 2000 || $25 < -2000) {print($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12)}' | sort -u | bedtools sort -i - > $outputDir/closest_enhancers.bed
  awk '{OFS="\t"}; {print($1,$2,$3,$4,$5,$6)}' $geneFile | bedtools sort -i - | bedtools closest -a $outputDir/closest_enhancers.bed -b - -D ref -N | awk '{OFS="\t"}; ($19 > 2000 || $19 < -2000) {print($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12)}' | sort -u > $path/closest_genes.bed
  #Removes duplicates - if there are exactly the same genes and random chromosomes.
  awk '(length($1) < 6)' $path/closest_genes.bed > $outputDir/"$fileName"_sorted_uniq.bed 
  wait
  rm -f $outputDir/*TSS* $outputDir/closest_enhancers.bed $outputDir/closest_genes.bed
  
}

function extract_coordinates
{
  # This function gets the FANTOM5 enhancers coordinates and extract X bases on each side of the eTSS on the *coding* strand.
  bedFile=$1
  extension=$2
  outputDir=$3
  fileName=`basename $bedFile _sorted_uniq.bed`

  # Extend X nt from each side of the eTSSs.
  awk -v var=$extension '{OFS="\t"}; {print($1, $2-var, $2, $4, $5, "'+'")}' $bedFile > $outputDir/"$fileName"_"$extension"_left.bed # 750bp upstream the eTSS.
  awk -v var=$extension '{OFS="\t"}; {print($1, $3, $3+var, $4, $5, "'+'")}' $bedFile > $outputDir/"$fileName"_"$extension"_right.bed # 750bp downstream the eTSS.
  wait
  bedtools sort -i $outputDir/"$fileName"_"$extension"_left.bed > $outputDir/"$fileName"_"$extension"_left_sorted.bed
  bedtools sort -i $outputDir/"$fileName"_"$extension"_right.bed > $outputDir/"$fileName"_"$extension"_right_sorted.bed
  wait 
  rm -f $outputDir/"$fileName"_"$extension"_left.bed $outputDir/"$fileName"_"$extension"_right.bed # Remove unneccesary files.  
}


mkdir -p $path/eTSS_analysis_data/filtered_enhancers
filter_enhancers $path/hg38_enhancers.bed $genesFile_hg38 $path/eTSS_analysis_data/filtered_enhancers
wait
mkdir -p $path/eTSS_analysis_data/stranded_enhancers
extract_coordinates $path/filtered_enhancers/hg38_enhancers_sorted_uniq.bed 750 $path/eTSS_analysis_data/stranded_enhancers # Expand 750 bp for the average profiles analysis.


#!/bin/bash
# bedtools version v2.26.0
###################################################################################################################################################
# This script gets bed files containing enhancer coordinates defined by FANTOM5 project and filter enhancers that overlap other enhancers or genes.
###################################################################################################################################################

path=your_data_path
genesFile=hg38_RefSeq_genes_17_12_19.bed
mkdir -p $path/filtered_enhancers
for bedFile in $path/*.bed
do
  fileName=`basename $bedFile .bed`
  # In case several enhancers have the same TSS, keep the longest one.(using the sort commands).
  cat $bedFile | sort -rn -k5,5 | sort -u -k1,1 -k2,2 | sortBed -i - > $path/sorted_uniqTSS.bed  

  # Find overlaps and report the number of merged intervals. 
  bedtools merge -i $path/sorted_uniqTSS.bed -c 4,4,6 -o collapse,count,distinct > $path/merged_uniqTSS.bed 
  
  #Keeps only those that don't overlap (have value of '1' in the fifth column obtained by bedtools merge).
  awk '($5 == '1')' $path/merged_uniqTSS.bed > $path/no_overlaps.bed 
  wait
  # Remove enhancers that have other enhancers or coding genes at a distance of less than 2 Kb from each side.
  bedtools sort -i $bedFile | bedtools closest -a $path/sorted_uniqTSS.bed -b - -D ref -N | awk '{OFS="\t"}; ($25 > 2000 || $25 < -2000) {print($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12)}' | sort -u | bedtools sort -i - > $path/closest_enhancers.bed
  awk '{OFS="\t"}; {print($1,$2,$3,$4,$5,$6)}' $genesFile | bedtools sort -i - | bedtools closest -a $path/closest_enhancers.bed -b - -D ref -N | awk '{OFS="\t"}; ($19 > 2000 || $19 < -2000) {print($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12)}' | sort -u > $path/closest_genes.bed
  # Removes contig chromosomes.
  awk '(length($1) < 6)' $path/closest_genes.bed > $path/filtered_enhancers/"$fileName"_sorted_uniq.bed 
  wait
  #Remove unnecessary files.
  rm -f $path/*TSS* $path/closest_enhancers.bed $path/closest_genes.bed
done
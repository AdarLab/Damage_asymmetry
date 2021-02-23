#!/bin/bash

# bedtools v2.27.1

#####################################################################################################################
# This script prepares humen genomic elements for plotting RNA poymerase ChIP-seq data and ATAC-seq data across them.
# The first and last exons/introns of every genes are remoev in order to avoid bias.
#####################################################################################################################

dataPath=/RNA_pol_II_data/annotation
chromSizes=hg38.chrom.sizes

function remove_extrem_intervals
{
  elementFile=$1
  element=$2
  fileName=`basename $elementFile .bed`
  awk -F${element} '{OFS="\t"};{print $1,$2}' $elementFile | awk -F'_chr' '{OFS="\t"};{print $1,$2}' | awk '{OFS="\t"}; {print($1,$2,$3,$4,$5,$8)}' > $dataPath/tmp.bed
  # Generate bed file with the first and last intervals
  awk '($4 != lastid)  { if (last) print last; print $0; } { lastid = $4; last = $0  } END {print last }' $dataPath/tmp.bed > $dataPath/first_last.bed
  # Remove from the orignial file the first and last intervals.
  grep -vFf $dataPath/first_last.bed $dataPath/tmp.bed > $dataPath/"$fileName"_no_first_last.bed
  wait
  rm -f  $dataPath/first_last.bed $dataPath/tmp.bed
}

function generate_window_file
{
  elementFile=$1
  length=$2
  left_extension=$3
  right_extestion=$4
  
  
  fileName=`basename $elementFile .bed`  
  bedtools sort -i $elementFile | awk -v var=$length '($3-$2 > var){print;}' | awk '{OFS="\t"} ; {if($6 == "'+'") {print($1,$2,$2,$4,$5,$6)} else {print($1,$3,$3,$4,$5,$6)}}' | bedtools slop -i - -g $chromSizes -l $left_extension -r $right_extestion -s > $dataPath/averageProfilesData/"$fileName"_start_plot.bed

 bedtools sort -i $elementFile | awk -v var=$length '($3-$2 > var){print;}' | awk '{OFS="\t"} ; {if($6 == "'+'") {print($1,$3,$3,$4,$5,$6)} else {print($1,$2,$2,$4,$5,$6)}}' | bedtools slop -i - -g $chromSizes -l $right_extestion -r $left_extension -s > $dataPath/averageProfilesData/"$fileName"_end_plot.bed
}
remove_extrem_intervals $dataPath/chopped_exons.bed _exon_
remove_extrem_intervals $dataPath/chopped_introns.bed _intron_
wait
mkdir -p $dataPath/averageProfilesData

generate_window_file $dataPath/genes.bed 9999 3000 10000
generate_window_file $dataPath/chopped_exons_no_first_last.bed 119 0 120
generate_window_file $dataPath/chopped_introns_no_first_last.bed 1873 0 1874
wait
bedtools sort -i $human_elements_path/hg38_RefSeq_uniq_coding_genes.bed | awk '($3-$2 > 9999){print;}' | awk '{OFS="\t"} ; {if($6 == "'+'") {print($1,$2,$2,$4,$5,$6)} else {print($1,$3,$3,$4,$5,$6)}}' | bedtools slop -i - -g $chromSizes -l 3000 -r 10000 -s > $avgProfilesDir/annotation/hg38_RefSeq_uniq_coding_genes_3Kb_up_plot.bed

bedtools sort -i $human_elements_path/hg38_RefSeq_uniq_coding_genes.bed | awk '($3-$2 > 9999){print;}' | awk '{OFS="\t"} ; {if($6 == "'+'") {print($1,$3,$3,$4,$5,$6)} else {print($1,$2,$2,$4,$5,$6)}}' | bedtools slop -i - -g $chromSizes -l 10000 -r 3000 -s > $avgProfilesDir/annotation/hg38_RefSeq_uniq_coding_genes_3Kb_down_plot.bed

bedtools sort -i $human_elements_path/hg38_RefSeq_uniq_chopped_exons.bed | awk '($3-$2 > 119){print;}' | awk '{OFS="\t"} ; {if($6 == "'+'") {print($1,$2,$2,$4,$5,$6)} else {print($1,$3,$3,$4,$5,$6)}}' | bedtools slop -i - -g $chromSizes -l 0 -r 120 -s > $avgProfilesDir/annotation/hg38_RefSeq_uniq_chopped_exons_start_plot.bed

bedtools sort -i $human_elements_path/hg38_RefSeq_uniq_chopped_introns.bed | awk '($3-$2 > 1873){print;}' | awk '{OFS="\t"} ; {if($6 == "'+'") {print($1,$2,$2,$4,$5,$6)} else {print($1,$3,$3,$4,$5,$6)}}' | bedtools slop -i - -g $chromSizes -l 0 -r 1874 -s > $avgProfilesDir/annotation/hg38_RefSeq_uniq_chopped_introns_start_plot.bed

bedtools sort -i $human_elements_path/hg38_RefSeq_uniq_chopped_exons.bed | awk '($3-$2 > 119){print;}' | awk '{OFS="\t"} ; {if($6 == "'+'") {print($1,$3,$3,$4,$5,$6)} else {print($1,$2,$2,$4,$5,$6)}}' | bedtools slop -i - -g $chromSizes -l 120 -r 0 -s > $avgProfilesDir/annotation/hg38_RefSeq_uniq_chopped_exons_end_plot.bed

bedtools sort -i $human_elements_path/hg38_RefSeq_uniq_chopped_introns.bed | awk '($3-$2 > 1873){print;}' | awk '{OFS="\t"} ; {if($6 == "'+'") {print($1,$3,$3,$4,$5,$6)} else {print($1,$2,$2,$4,$5,$6)}}' | bedtools slop -i - -g $chromSizes -l 1874 -r 0 -s > $avgProfilesDir/annotation/hg38_RefSeq_uniq_chopped_introns_end_plot.bed

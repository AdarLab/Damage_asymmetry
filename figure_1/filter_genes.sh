#!/bin/bash
# bedtools version v2.26.0
#############################################################################################
# This script gets bed files containing genes coordinates defined by RefSeq and filters out
# genes that overlap other genes or at distance of less than 6 Kb upstream of adjacent genes.
#############################################################################################

path=your_path
# Renames the transcripts with unique names.
awk 'BEGIN { OFS = "\t"}; {print;}' $path/hg38_RefSeq_coding_genes_17_12_19.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"_"$1"::"$2":"$3"\t"$5"\t"$6}' > $path/genes_list.bed 
# If it is '+' strand: take the longest one of the same 'TSS'(using the sort comands).
grep '+' $path/genes_list.bed | sort -rn -k5,5 | sort -u -k1,1 -k2,2 > $path/uniqTSS.bed
# '-' strand.
grep '-' $path/genes_list.bed | sort -rn -k5,5 | sort -u -k1,1 -k3,3 >> $path/uniqTSS.bed 
# Sort the bed file.
sortBed -i $path/uniqTSS.bed > $path/sorted_uniqTSS.bed 
# Finds overlaps.
bedtools merge -i $path/sorted_uniqTSS.bed -c 4,4,6 -o collapse,count,distinct > $path/merged_uniqTSS.bed 
# Keep only those that don't overlap.
awk '{if($5 == '1')print($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6)}' $path/merged_uniqTSS.bed > $path/no_overlaps.bed 
#Find the closest genes for the genes that don't overlap:
# -id - ignore downstream closest genes -N - to avoid including the gene itself. 
bedtools closest -a $path/no_overlaps.bed -b $path/sorted_uniqTSS.bed -D ref -id -N > $path/closest_genes.bed
#Keeps only the genes that have no nearby genes 6 kb upstream
awk '{if($13 < '-6000')print($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6)}' $path/closest_genes.bed > $path/TSS_6000.bed
# Removes duplicates - if there are exactly the same genes and random chromosomes.
sort -u -k4,4 $path/TSS_6000.bed | awk '(length($1) < 6)' > $path/hg38_RefSeq_uniq_coding_genes.bed 
#Extracts the list of the specific genes IDs.
awk '{print $4}' $path/hg38_RefSeq_uniq_coding_genes.bed | cut -d '_' -f 1-2 > $path/uniq_coding_genes_list.txt
wait
#Remove unnecessary files.
rm -f $path/*TSS* $path/closest_genes.bed $path/no_overlaps.bed $path/genes_list.bed 

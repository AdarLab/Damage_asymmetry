#!/bin/bash
# bedtools version v2.26.0
###################################################################################################################
# This script gets bed files of genomic elements and generates windows for plotting profiles of CPD/TT frequencies.
###################################################################################################################

chromSizes=hg38.chrom.sizes
dataDir=your_data_directory
human_elements_path=$dataDir/annotation
damage_seq_bedFiles_path=$dataDir/filteredReads/bedFiles

mkdir -p $dataDir/averageProfilesData/annotation
avgProfilesDir=$dataDir/averageProfilesData

# Use bedtools slop for extracting 3 Kb upstream and 10 Kb downstream of the TSS of genes.
bedtools sort -i $human_elements_path/hg38_RefSeq_uniq_coding_genes.bed | awk '($3-$2 > 9999){print;}' | awk '{OFS="\t"} ; {if($6 == "'+'") {print($1,$2,$2,$4,$5,$6)} else {print($1,$3,$3,$4,$5,$6)}}' | bedtools slop -i - -g $chromSizes -l 3000 -r 10000 -s > $avgProfilesDir/annotation/hg38_RefSeq_uniq_coding_genes_3Kb_up_plot.bed # In case the strand is '+', the TSS is the second column of the bed file.

bedtools sort -i $human_elements_path/hg38_RefSeq_uniq_coding_genes.bed | awk '($3-$2 > 9999){print;}' | awk '{OFS="\t"} ; {if($6 == "'+'") {print($1,$3,$3,$4,$5,$6)} else {print($1,$2,$2,$4,$5,$6)}}' | bedtools slop -i - -g $chromSizes -l 10000 -r 3000 -s > $avgProfilesDir/annotation/hg38_RefSeq_uniq_coding_genes_3Kb_down_plot.bed # In case the strand is '-', the TSS is the second column of the bed file.

# Extract 120 bp (median length) downstream the exons' start.
bedtools sort -i $human_elements_path/hg38_RefSeq_uniq_chopped_exons.bed | awk '($3-$2 > 119){print;}' | awk '{OFS="\t"} ; {if($6 == "'+'") {print($1,$2,$2,$4,$5,$6)} else {print($1,$3,$3,$4,$5,$6)}}' | bedtools slop -i - -g $chromSizes -l 0 -r 120 -s > $avgProfilesDir/annotation/hg38_RefSeq_uniq_chopped_exons_start_plot.bed
# Extract 1874 bp (median length) downstream the introns' start.
bedtools sort -i $human_elements_path/hg38_RefSeq_uniq_chopped_introns.bed | awk '($3-$2 > 1873){print;}' | awk '{OFS="\t"} ; {if($6 == "'+'") {print($1,$2,$2,$4,$5,$6)} else {print($1,$3,$3,$4,$5,$6)}}' | bedtools slop -i - -g $chromSizes -l 0 -r 1874 -s > $avgProfilesDir/annotation/hg38_RefSeq_uniq_chopped_introns_start_plot.bed
# Extract 120 bp (median length) upstream the exons' end.
bedtools sort -i $human_elements_path/hg38_RefSeq_uniq_chopped_exons.bed | awk '($3-$2 > 119){print;}' | awk '{OFS="\t"} ; {if($6 == "'+'") {print($1,$3,$3,$4,$5,$6)} else {print($1,$2,$2,$4,$5,$6)}}' | bedtools slop -i - -g $chromSizes -l 120 -r 0 -s > $avgProfilesDir/annotation/hg38_RefSeq_uniq_chopped_exons_end_plot.bed
# Extract 1874 bp (median length) upstream the introns' end.
bedtools sort -i $human_elements_path/hg38_RefSeq_uniq_chopped_introns.bed | awk '($3-$2 > 1873){print;}' | awk '{OFS="\t"} ; {if($6 == "'+'") {print($1,$3,$3,$4,$5,$6)} else {print($1,$2,$2,$4,$5,$6)}}' | bedtools slop -i - -g $chromSizes -l 1874 -r 0 -s > $avgProfilesDir/annotation/hg38_RefSeq_uniq_chopped_introns_end_plot.bed

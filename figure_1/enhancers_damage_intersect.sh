#!/bin/bash
#This script gets damage seq bed files of 4 nucleotides, bed files of enhancers and count the number of appearances of the 4 nucleotides in every element.
#Bedtools version: v2.26.0

dataDir=/vol/sci/bio/data/sheera.adar/elisheva.h/enhancers
elementFile=$dataDir/stranded_enhancers/hg38_enhancers_plot_sorted.bed
damage_seq_bedFiles_path_hg19=/vol/sci/bio/data/sheera.adar/elisheva.h/Damages_Dinucleotides_comparison/filteredReads/bedFiles
damage_seq_bedFiles_path_hg38=/vol/sci/bio/data/sheera.adar/elisheva.h/Damages_Dinucleotides_comparison_hg38_data/filteredReads/bedFiles
damage_seq_naked_DNA_file=/vol/sci/bio/data/sheera.adar/elisheva.h/enhancers/filteredReadsNakedDNA/bedFiles/NCPD_naked_DNA_repA_repB_sorted_4_nt_dipyrimidine.bed
pattern_path_hg19=/vol/sci/bio/data/sheera.adar/elisheva.h/Damages_Dinucleotides_comparison/bioconductorProfiles/patterns/wholeGenome
pattern_path_hg38=/vol/sci/bio/data/sheera.adar/elisheva.h/Damages_Dinucleotides_comparison_hg38_data/patterns/wholeGenome
hg19_genome=/cs/icore/elisheva.h/hg19-ref_genome/hg19.chrom.sizes
hg38_genome=/cs/icore/elisheva.h/hg38-ref_genome/hg38.chrom.sizes


mkdir -p $dataDir/damageIntersect

intersectDir=$dataDir/damageIntersect
declare -a versions=("hg19" "hg38") #An array for damage type

#Scans all the genomic elements, count how many damages are in each element and normalize the reasult to the element length.
for version in "${versions[@]}"
do
  for bedFile in $dataDir/stranded_enhancers/*/*"$version"_enhancers_plot*
  do
    dir=`dirname $bedFile`
    dirName=`basename $dir`
    mkdir -p $intersectDir/$dirName/transcribed/$version
    mkdir -p $intersectDir/$dirName/nonTranscribed/$version 
    
    if [ "$version" == hg19 ]
    then
      sbatch intersect.sh $bedFile $damage_seq_bedFiles_path_hg19/NCPD_con_RepA_RepB_sorted_4nt.bed $intersectDir/$dirName/transcribed/"$version" $intersectDir/$dirName/nonTranscribed/"$version" $hg19_genome
      sbatch intersect.sh $bedFile $damage_seq_bedFiles_path_hg19/NCPD_T0_RepA_RepB_sorted_4nt.bed $intersectDir/$dirName/transcribed/"$version" $intersectDir/$dirName/nonTranscribed/"$version" $hg19_genome
      sbatch intersect.sh $bedFile $pattern_path_hg19/TT.bed $intersectDir/$dirName/transcribed/"$version" $intersectDir/$dirName/nonTranscribed/"$version" $hg19_genome
    fi
    if [ "$version" == hg38 ]
    then
      sbatch intersect.sh $bedFile $damage_seq_bedFiles_path_hg38/NCPD_con_RepA_RepB_sorted_uniq_dipyrimidine.bed $intersectDir/$dirName/transcribed/"$version" $intersectDir/$dirName/nonTranscribed/"$version" $hg38_genome
      sbatch intersect.sh $bedFile $damage_seq_bedFiles_path_hg38/NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed $intersectDir/$dirName/transcribed/"$version" $intersectDir/$dirName/nonTranscribed/"$version" $hg38_genome
      sbatch intersect.sh $bedFile $damage_seq_naked_DNA_file $intersectDir/$dirName/transcribed/"$version" $intersectDir/$dirName/nonTranscribed/"$version" $hg38_genome
      sbatch intersect.sh $bedFile $pattern_path_hg38/TT.bed $intersectDir/$dirName/transcribed/"$version" $intersectDir/$dirName/nonTranscribed/"$version" $hg38_genome
    fi
  done 
done 
  


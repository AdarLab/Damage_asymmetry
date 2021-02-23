#!/bin/bash
# bedtools version v2.26.0
############################################################################################################################
# This script gets Damage-seq bed files of 4 nucleotides, XR-seq files of 4 nucleotides,
# bed files of humen enhancers coordinates and extract the damages the interseceted the given elements for further analyses.
############################################################################################################################


dataDir=your_path
enhancersFile=hg38_enhancers_sorted_uniq.bed
damageFile=NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed
repairPath=/XR_Seq_data/bedFiles
ref_genome_hg38=hg38.fa
tables_script=tables.py

mkdir -p $dataDir/coverageData/damage
mkdir -p $dataDir/coverageData/repair
coverageDir=$dataDir/coverageData
fileName=`basename $enhancersFile _sorted_uniq.bed`

# Extend 500 nt from each side of the eTSSs of the filtered file on the non-transceibed (coding) strand.
awk -v var=500 '{OFS="\t"}; {print($1, $2-var, $2, $4, $5, "'-'")}' $enhancersFile > $coverageDir/"$fileName"_500.bed # Extend 500 bases upstream the TSS and define the strand as minus.
awk -v var=500 '{OFS="\t"}; {print($1, $3, $3+var, $4, $5, "'+'")}' $enhancersFile >> $coverageDir/"$fileName"_500.bed # Extend 500 bases downstream the TSS and define the strand as plus.
wait
bedtools sort -i $coverageDir/"$fileName"_500.bed > $coverageDir/"$fileName"_500_sorted.bed
wait
rm -f $coverageDir/"$fileName"_500.bed # Remove unneccesary files
wait
enhancersStrandedFile=$coverageDir/"$fileName"_500_sorted.bed

################
#--- Damage ---#
################

# Scan all the genomic elements, count how many damages are in each element and normalize the reasult to the element length.
mkdir -p $coverageDir/damage/transcribedStrand
mkdir -p $coverageDir/damage/nonTranscribedStrand 
damageFileName=`basename $damageFile`
echo "$damageFileName"
bedtools coverage -s -a $enhancersStrandedFile -b $damageFile | awk '{print($0=$0"\t"($7/($3-$2)))}' > $coverageDir/damage/nonTranscribedStrand/enhancers_hg38_"$damageFileName" #-s: select only the reads on the same strand of the original element.
bedtools coverage -S -a $enhancersStrandedFile -b $damageFile | awk '{print($0=$0"\t"($7/($3-$2)))}' > $coverageDir/damage/transcribedStrand/enhancers_hg38_"$damageFileName" #-s: select only the reads on the opposite strand.

################
#--- Repair ---#
################
mkdir -p $coverageDir/repair/transcribedStrand
mkdir -p $coverageDir/repair/nonTranscribedStrand
# Scan all the genomic elements, count how many repair reads obtained by XR-seq are in each element and normalize the reasult to the element length.
for repairFile in $repairPath/*.bed
do
  repairFileName=`basename $repairFile`
  bedtools coverage -s -a $enhancersStrandedFile -b $repairFile | awk '{print($0=$0"\t"($7/($3-$2)))}' > $coverageDir/repair/nonTranscribedStrand/enhancers_hg38_"$repairFileName" #-s: select only the reads on the same strand of the original element.
  bedtools coverage -S -a $enhancersStrandedFile -b $repairFile | awk '{print($0=$0"\t"($7/($3-$2)))}' > $coverageDir/repair/transcribedStrand/enhancers_hg38_"$repairFileName" #-s: select only the reads on the opposite strand.
done

mkdir -p $coverageDir/sequence
bedtools getfasta -s -fi $ref_genome_hg38 -bed $enhancersStrandedFile -name -fo $coverageDir/sequence/hg_38_enhancers_seq.fa
wait

# Run the tables script
mkdir -p $coverageDir/tables/counting
for file in $coverageDir/sequence/*.fa
do
  name=`basename $file`
  python $tables_script $coverageDir/tables $file $name
done


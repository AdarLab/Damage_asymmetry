#!/bin/bash
# bedtools version v2.26.0
##################################################################################################
# This script gets Damage-seq bed files of 4 nucleotides and filters only those with dipyrimidins.
##################################################################################################
currentDir=your_path
ref_genome=hg38.fa #hg-38 ref genome.

#Scans all the damage seq given bed files.
mkdir -p $currentDir/reads/fastaFiles
for bedFile in $currentDir/reads/bedFiles/*.bed
do
 name=`basename $bedFile`

# Adjust the fasta file name. 
 faName="${name%.bed*}.fa"    
 bedtools getfasta -s -fi $ref_genome -bed $bedFile -name -fo $currentDir/reads/fastaFiles/$faName #Get the sequence itself.
done
wait

mkdir -p $currentDir/filteredReads/bedFiles
mkdir -p $currentDir/filteredReads/fastaFiles

#Scans all the obtained fasta files.
for faFile in $currentDir/reads/fastaFiles/*.fa
do
  name=`basename $faFile`
  faOutput="${name%.fa*}_dipyrimidine.fa"
  bedOutput="${name%.fa*}_dipyrimidine.bed"
  #Filters the files for getting only those which contain dypirimidnes. 
  cat $faFile  |
  paste - -  |
  awk 'toupper($2) ~ /(TT|TC|CT|CC)/' |
  tr "\t" "\n" |
  tee $currentDir/filteredReads/fastaFiles/$faOutput | 
  awk -F '[>():-]' '/^>/ {printf("%s\t%s\t%s\t%s\t%s\t%s\n",$4,$5,$6,$2,$2,$7);}' | 
  awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "-" }; 1' > $currentDir/filteredReads/bedFiles/$bedOutput 
done

#Scan all the files to get an indication of the number of sequences left after the filter.
for file1 in $currentDir/reads/bedFiles/*.bed
do
  for file2 in $currentDir/filteredReads/bedFiles/*.bed
  do
    name1=`basename $file1`
    name2=`basename $file2`

   if [[ *$name1* == $name2 ]]; then # The original file and the adugsred filtered file.
      raw=$(cat $file1 | wc -l) # The number of sequences in the original files.
      filtered=$(cat $file2 | wc -l) # The number of sequences in the filtered files.
      result=$(awk "BEGIN {printf \"%.2f\",${filtered}/${raw} *100}") # Calculate the percentage of the remained sequences.
      printf '%s%s%s%s\n' $name1 ',' $result '%' >> $currentDir/summary_table.csv # Write the results into a table.
      
    fi
  done
done



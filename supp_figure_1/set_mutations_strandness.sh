#!/bin/bash
####################################################################################
# This script scan the mutations files and merge them to a single stranded bed file.
####################################################################################
mutationPath=mutation_data
C_T_path=$mutationPath/C_T
G_A_path=$mutationPath/G_A

for c_t_file in $C_T_path/*.csv
do
  fileName=`basename $c_t_file .csv`
  cancerName="${fileName/CtoT/}"
  sed '1d' $c_t_file | sed 's/[^[:blank:][:print:]]//g' | tr -d '"' | awk -F, '{OFS="\t"};{print("chr"$2,$3-1,$4,"1","C/T","+")}' > $mutationPath/stranded_tables/"$fileName"_stranded.bed # Define the strand of C->T as '+'.
  g_a_file=$G_A_path/*"$cancerName"*
  sed '1d' $g_a_file | sed 's/[^[:blank:][:print:]]//g' | tr -d '"' | awk -F, '{OFS="\t"};{print("chr"$2,$3-1,$4,"1","G/A","-")}' >> $mutationPath/stranded_tables/"$fileName"_stranded.bed # Define the strand of G->A as '-' (complementary sequence).
done
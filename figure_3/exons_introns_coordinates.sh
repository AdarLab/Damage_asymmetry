#!/bin/bash
# bedtools v2.27.1
################################################################################################################
# This script gets exons and introns coordinates and genartaes couting tables of the filtered exons and introns.
################################################################################################################

path=your_path
ref_genome=hg38.fa
tables_script=tables.py

#This function gets a bed file and returns a FASTA file using bedtools get fasta command.
function get_fasta
{
	bedtools getfasta -s -fi $ref_genome -bed $1 -name -fo $path/sequence/$2
}

# In case of isoforms, keep only the longest vriants.
grep '+' $path/annotation/raw_data/hg38_RefSeq_exons.bed | sort -rn -k5,5 | sort -u -k1,1 -k2,2 > $path/annotation/hg38_RefSeq_uniq_exons.bed
grep '-' $path/annotation/raw_data/hg38_RefSeq_exons.bed | sort -rn -k5,5 | sort -u -k1,1 -k3,3 >> $path/annotation/hg38_RefSeq_uniq_exons.bed

grep '+' $path/annotation/raw_data/hg38_RefSeq_introns.bed | sort -rn -k5,5 | sort -u -k1,1 -k2,2 > $path/annotation/hg38_RefSeq_uniq_introns.bed
grep '-' $path/annotation/raw_data/hg38_RefSeq_introns.bed | sort -rn -k5,5 | sort -u -k1,1 -k3,3 >> $path/annotation/hg38_RefSeq_uniq_introns.bed

wait

#Choppes 100 bases from each side of the introns.
awk '{
	if($3 - $2 > 200)
		print($1"\t"$2 + 100"\t"$3 - 100"\t"$4"\t"$5"\t"$6);
	}' $path/annotation/hg38_RefSeq_uniq_introns.bed > $path/annotation/hg38_RefSeq_uniq_chopped_introns.bed
#Choppes 10 bases from each side of the exons.
awk '{
	if($3 - $2 > 20)
		print($1"\t"$2 + 10"\t"$3 - 10"\t"$4"\t"$5"\t"$6);
	}' $path/annotation/hg38_RefSeq_uniq_exons.bed > $path/annotation/hg38_RefSeq_uniq_chopped_exons.bed
#Delet haploid and randon chromosomes.
for bedFile in $path/annotation/*.bed
do
  gawk -i inplace '(length($1) < 6) {print $0}' $bedFile
done

# Chopped exons
get_fasta $path/annotation/hg38_RefSeq_uniq_chopped_exons.bed "chopped_exons_seq.fa"

# Chopped introns
get_fasta $path/annotation/hg38_RefSeq_uniq_chopped_introns.bed "chopped_introns_seq.fa" 

wait
# Runs the tables script
mkdir -p $path/tables/counting
mkdir -p $path/tables/ratio
for file in $path/sequence/*.fa
do
	name=`basename $file`
	python $tables_script $path/tables $file $name
done



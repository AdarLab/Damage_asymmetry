#!/bin/bash
# bedtools v2.27.1
##############################################################################################
# This script gets expression data by GTEx and define the median TPM value of each transcript.
##############################################################################################
dataDir=your_path
gff_file=$dataDir/gencode.v25.annotation.gtf
exons_intrnos_script=gff3_to_exons_introns.py
refGenome=hg38.fa

mkdir -p $dataDir/introns
mkdir -p $dataDir/exons

# Run the python script for extracting exons and introns from the gtf file.
python $exons_intrnos_script -i $gff_file -o $dataDir/introns/gencode.v25.introns.annotation.gtf $dataDir/exons/gencode.v25.exons.annotation.gtf
wait
function get_fasta
{  
  bedFile=$1
  dirName=`dirname $bedFile`
  fileName=`basename $bedFile`
  faName="${fileName%.bed}.fa"
  mkdir -p $dirName/tables/counting
  mkdir -p $dirName/tables/ratio
  mkdir -p $dirName/sequence 
  
  # Extract the sequence
  bedtools getfasta -fi $refGenome -bed $bedFile -s -name > $dirName/sequence/$faName" 
  wait
  python tables.py $dirName/tables $dirName/sequence/$faName $faName"
  # In case the elements are exons or introns, merge them first to get all the exons/introns of the same gene as one entry.
   merge_elements.sh $dirName $faName
   wait 
   python tables.py $dirName/tables $dirName/sequence/continuous_$faName continuous_$faName"

}

awk '{print $1}' $dataDir/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm_Testis.txt > $dataDir/GTEx_Analysis_Testis_genes_names.tsv # Extract the gene Ids.

# Extract from the gtf file only the protein coding genes.
grep -wFf $dataDir/GTEx_Analysis_Testis_genes_names.tsv $gff_file | awk '($3 == "gene"){print;}' | grep 'protein_coding' | awk '{OFS="\t"} {print ( $1,$4,$5,$10,$16,$7)}' | tr -d '"' | tr -d ';' > $dataDir/gencode.v25.coding_genes.annotation.bed

# In case several genes have the same TSS on the same strand select only the longest gene.
awk '($6 == "'+'")' $dataDir/gencode.v25.coding_genes.annotation.bed | sort -n -k2,2 | sort -u -k1,1 -k2,2 > $dataDir/gencode.v25.coding_genes.annotation_uniq_TSS.bed
awk '($6 == "'-'")' $dataDir/gencode.v25.coding_genes.annotation.bed | sort -r -n -k3,3 | sort -u -k1,1 -k3,3 >> $dataDir/gencode.v25.coding_genes.annotation_uniq_TSS.bed

# Remove overlapping genes using bedtools cluster.
bedtools sort -i $dataDir/gencode.v25.coding_genes.annotation_uniq_TSS.bed | bedtools cluster -i - -d 0 | sort -u -k7,7 | awk '{OFS="\t"}; {print($1,$2,$3,$4,$5,$6)}' > $dataDir/uniq_coding_genes.bed

rm -f $dataDir/gencode.v25.coding_genes.annotation_uniq_TSS.bed

awk '{print $5"\t"$4}' $dataDir/uniq_coding_genes.bed > $dataDir/uniq_coding_genes_symbols.txt

sort -k1,1 $dataDir/uniq_coding_genes_symbols.txt > $dataDir/uniq_coding_genes_symbols_sorted.txt
awk '{OFS=FS="\t"}; {print($2,$1,$3)}' $dataDir/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm_Testis.txt | sort -k1,1 > $dataDir/GTEx_Analysis_gene_median_tpm_Testis_sorted.txt
# Merge the unique Ids and the expression values based on the genes Ids.
join -j 1 -o 1.1 1.2 2.3 $dataDir/uniq_coding_genes_symbols_sorted.txt $dataDir/GTEx_Analysis_gene_median_tpm_Testis_sorted.txt | awk '{OFS="\t"};{print($1,$2,$3)}' > $dataDir/GTEx_Analysis_gene_median_tpm_Testis_joined.txt
# Extract the unique Ids.
awk '{print $2}' $dataDir/GTEx_Analysis_gene_median_tpm_Testis_joined.txt > $dataDir/uniq_coding_genes_Ids.txt
rm -f $dataDir/uniq_coding_genes_symbols.txt $dataDir/GTEx_Analysis_gene_median_tpm_Testis_sorted.txt
# Keeps only the longest transcripst.
grep -wFf $dataDir/uniq_coding_genes_Ids.txt $gff_file | awk '{OFS="\t"};($3 == "transcript") {print ( $1,$4,$5,$10,$12,$7,$5-$4)}' | tr -d '"' | tr -d ';' | sort -r -n -k7,7 | sort -u -k4,4 | awk '{print $5}' > $dataDir/uniq_coding_trnascripts_Ids.txt

wait

# Extract the matching exons and introns.
grep -wFf $dataDir/uniq_coding_trnascripts_Ids.txt $dataDir/introns/gencode.v25.introns.annotation.gtf | awk '{OFS="\t"}; {print($1,$4,$5,$8,$9,$7)}' > $dataDir/introns/uniq_coding_introns.bed
grep -wFf $dataDir/uniq_coding_trnascripts_Ids.txt $dataDir/exons/gencode.v25.exons.annotation.gtf |  awk '{OFS="\t"}; {print($1,$4,$5,$8,$9,$7)}' > $dataDir/exons/uniq_coding_exons.bed

get_fasta $dataDir/uniq_coding_genes.bed
get_fasta $dataDir/introns/uniq_coding_introns.bed
get_fasta $dataDir/exons/uniq_coding_exons.bed



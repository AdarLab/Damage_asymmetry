#!/bin/bash
# bedtools Version: v2.26.0

###################################################################################################################
# This script gets human PGC cells expression data and geneartes counting tables matching to the given transcripts.
###################################################################################################################

dataDir=sperm
gff_file=$dataDir/gencode.v19.annotation.gtf
exons_intrnos_script=gff3_to_exons_introns.py
refGenome=hg_19_ref_genome.fa

function get_fasta
{  
  bedFile=$1
  dirName=`dirname $bedFile`
  fileName=`basename $bedFile`
  faName="${fileName%.bed}.fa"
  mkdir -p $dirName/tables/counting
  mkdir -p $dirName/sequence 
  
  # Extract the sequence
  bedtools getfasta -fi $refGenome -bed $bedFile -s -name > $dirName/sequence/$faName" 
  wait
  tables.py $dirName/tables $dirName/sequence/continuous_$faName continuous_$faName
}

mkdir -p $dataDir/introns
mkdir -p $dataDir/exons

# Run the python script for extracting exons and introns from the gtf file.
python $exons_intrnos_script -i $gff_file -o $dataDir/introns/gencode.v19.introns.annotation.gtf $dataDir/exons/gencode.v19.exons.annotation.gtf
wait
awk '{print $1}' $dataDir/PGC_bulk_avg_FPKM.tsv > $dataDir/PGC_bulk_genes_names.tsv # Extract the gene Ids.

# Extract from the gtf file only the protein coding genes.
grep -wFf $dataDir/PGC_bulk_genes_names.tsv $dataDir/gencode.v19.annotation.gtf | awk '($3 == "gene"){print;}' | grep 'protein_coding' | awk '{OFS="\t"} {print ( $1,$4,$5,$10,$18,$7)}' | tr -d '"' | tr -d ';' > $dataDir/gencode.v19.coding_genes.annotation.bed

# In case several genes have the same TSS on the same strand select only the longest gene.
awk '($6 == "'+'")' $dataDir/gencode.v19.coding_genes.annotation.bed | sort -n -k2,2 | sort -u -k1,1 -k2,2 > $dataDir/gencode.v19.coding_genes.annotation_uniq_TSS.bed
awk '($6 == "'-'")' $dataDir/gencode.v19.coding_genes.annotation.bed | sort -r -n -k3,3 | sort -u -k1,1 -k3,3 >> $dataDir/gencode.v19.coding_genes.annotation_uniq_TSS.bed

# Remove overlapping genes using bedtools cluster.
bedtools sort -i $dataDir/gencode.v19.coding_genes.annotation_uniq_TSS.bed | bedtools cluster -i - -d 0 | sort -u -k7,7 | awk '{OFS="\t"}; {print($1,$2,$3,$4,$5,$6)}' > $dataDir/uniq_coding_genes.bed

rm -f $dataDir/gencode.v19.coding_genes.annotation_uniq_TSS.bed

awk '{print $5"\t"$4}' $dataDir/uniq_coding_genes.bed > $dataDir/uniq_coding_genes_symbols.txt
sort -k1,1 $dataDir/uniq_coding_genes_symbols.txt > $dataDir/uniq_coding_genes_symbols_sorted.txt
sort -k1,1 $dataDir/PGC_bulk_avg_FPKM.tsv > $dataDir/PGC_bulk_avg_FPKM_sorted.tsv
join -j 1 -o 1.1 1.2 2.2 $dataDir/uniq_coding_genes_symbols_sorted.txt $dataDir/PGC_bulk_avg_FPKM_sorted.tsv | awk '{OFS="\t"};{print($1,$2,$3)}' > $dataDir/PGC_bulk_avg_FPKM_joined.tsv
awk '{print $2}' $dataDir/PGC_bulk_avg_FPKM_joined.tsv > $dataDir/uniq_coding_genes_Ids.txt
rm -f $dataDir/uniq_coding_genes_symbols.txt $dataDir/PGC_bulk_avg_FPKM_sorted.tsv $dataDir/PGC_bulk_avg_FPKM_sorted.tsv
# Keeps only the longest transcripst.
grep -wFf $dataDir/uniq_coding_genes_Ids.txt $dataDir/gencode.v19.annotation.gtf | awk '{OFS="\t"};($3 == "transcript") {print ( $1,$4,$5,$10,$12,$7, $5-$4)}' | tr -d '"' | tr -d ';' | sort -r -n -k7,7 | sort -u -k4,4 | awk '{print $5}' > $dataDir/uniq_coding_trnascripts_Ids.txt

wait

# Extract the matching exons and introns.
grep -wFf $dataDir/uniq_coding_trnascripts_Ids.txt $dataDir/introns/gencode.v19.introns.annotation.gtf | awk '{OFS="\t"}; {print($1,$4,$5,$8,$9,$7)}' > $dataDir/introns/uniq_coding_introns.bed
grep -wFf $dataDir/uniq_coding_trnascripts_Ids.txt $dataDir/exons/gencode.v19.exons.annotation.gtf |  awk '{OFS="\t"}; {print($1,$4,$5,$8,$9,$7)}' > $dataDir/exons/uniq_coding_exons.bed
wait

get_fasta $dataDir/introns/uniq_coding_introns.bed
get_fasta $dataDir/exons/uniq_coding_exons.bed

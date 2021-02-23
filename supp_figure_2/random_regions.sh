#!/bin/bash 
# bedtools Version: v2.26.0

#############################################################################################################
# This script gets human genomic elements bed files and generate bed files containing randm coordinates
# at the same size of the orginal elements but at different genomic poistions using bedtools shuffle command.
#############################################################################################################

dataDir=your_path
beds_path=$dataDir/annotation #The bed files path
chromSizes=hg38.chrom.sizes
refGenome=hg38.fa
rawGenesFile=hg38_RefSeq_genes_17_12_19.bed
awk '{OFS=FS="\t"}; {print($1,$2,$3,$4,$5,$6)}' $rawGenesFile > $dataDir/random_regions/transcribed_regions.bed
genesFile=$dataDir/random_regions/transcribed_regions.bed
NsFile=$dataDir/random_regions/unknown_sequences.bed
cat $genesFile $NsFile > $dataDir/random_regions/excluded_regions.bed # Define the genes and unknown sequences as excluded regions.
excludedRegionsFile=$dataDir/random_regions/excluded_regions.bed

mkdir -p $dataDir/random_regions/annotation
mkdir -p $dataDir/random_regions/sequence
mkdir -p $dataDir/random_regions/tables/counting
# Extract random regions that are not transcribed at the same number and length of the given elements. 
shuffile_job_Id=""
for file in $beds_path/*.bed
do
	raw_name=`basename $file`
	name=${raw_name::-4}
  echo "$name"
bedtools shuffle -i $file -g $chromSizes -excl $excludedRegionsFile -maxTries 10000000  > $dataDir/random_regions/annotation/"${name}_rand.bed"
done
wait
shuffile_job_Id="${shuffile_job_Id::-1}"
get_fasta_job_Id=""
# Get the sequence itself.
for file in $dataDir/random_regions/annotation/*.bed
do
	raw_name=`basename $file`
	name=${raw_name::-4}
	bedtools getfasta -s -fi $refGenome -bed $file -name -fo $dataDir/random_regions/sequence/"${name}.fa"
done
wait
get_fasta_job_Id="${get_fasta_job_Id::-1}"
# Run the python scripts
for file in $dataDir/random_regions/sequence/*.fa
do
	name=`basename $file`
  python tables.py $dataDir/random_regions/tables $file $name
done
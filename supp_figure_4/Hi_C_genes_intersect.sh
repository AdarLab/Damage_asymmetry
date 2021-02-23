#!/bin/bash
# Bedtools version: v2.26.0

######################################################################################
# This script gets A/B compartment score of Hi-C data of H9 ESC,
# human genes coordinates and find which genes are relies in the A/B genomic windows.
######################################################################################

data_dir=Hi_C_data
intersected_dir=$data_dir/intersected_data
A_B_c_file=$data_dir/H9_ESC/H9_ESC_Hi_C_A_B_compartment.bed
genes_file=hg38_RefSeq_uniq_coding_genes.bed

bedtools intersect -a $A_B_c_file -b $genes_file -wa -wb > $intersected_dir/H9_ESC_Hi_C_A_B_compartment_overlapped_genes.bed # -wa & -wb: write the entries in both genomic windows and genes for each overlap.
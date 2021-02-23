#!/usr/bin/Rscript 
# R version 3.5.2
library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(ggpointdensity)
library(tidyverse)

#############################################################################################################################################
# This script gets tables of A/B compartments scores based on Hi-C data across 50 kb genomic windows,
# Classified human elements based on the A/B score and generates box plots representing the asymmetry score distribution in each compartment.
#############################################################################################################################################

intersectedfilePath <- "H9_ESC_Hi_C_A_B_compartment_overlapped_genes.bed"
exonsPath <- "continuous_chopped_exons_seq_counting.csv"
intronsPath <- "continuous_chopped_introns_seq_counting.csv"
intersectedFile <- read.table(intersectedfilePath, header = F, sep = "\t")
exonsCountingTable <- read.csv(exonsPath, header = T, sep = ",")
intronsCountingTable <- read.csv(intronsPath, header = T, sep = ",")
figurePath <- "figures"

adjust_p_value <- function(data)
{
  "This function gets data of 2 populations and compare them by using Wilcoxon test according to the givven alternative hypothesis"
  customSymnum <- list(cutpoints = c(0, c(0.0001, 0.001, 0.01, 0.05)/(nrow(data)/2), 1.01), symbols = c("*", "***", "**", "*", "n.s.")) # Define the p values symbols based on Bonferroni correction.
  
  A_compartment <- subset(data, Compartment == "A")$avg_asymmetry_score # Define the group of the transcribed strand.
  B_compartment <- subset(data, Compartment == "B")$avg_asymmetry_score # "" non transcribed strand.
  wilcoxon_stats <- wilcox.test(A_compartment, B_compartment, paired = F, alternative = "less")
  p_value <- wilcoxon_stats$p.value
  # Scan the customSymnum list in order to adjust the matching symbol to the calculated p value.
  for(i in 1:(length(customSymnum$cutpoints)-1)){
    if(customSymnum$cutpoints[i] <= p_value && p_value < customSymnum$cutpoints[i+1]){
      symbol <- customSymnum$symbols[i]
    }
  }
  return(symbol)
}

asymmetry_score <- function(counting_table) {
  asymmetryScore <- (counting_table[, "AA"] - counting_table[, "TT"]) / (counting_table[, "AA"] + counting_table[, "TT"])
  asymmetryTable <- cbind(counting_table[, "Id"], as.data.frame(asymmetryScore))
  names(asymmetryTable) <- c("Id", "asymmetryScore")


  return(asymmetryTable)
}

mergeTables <- function(hi_c_table, asymmetry_table) {
  hi_c_table[, "V10"] <- gsub("_chr.*", "", hi_c_table[, "V10"])
  hi_c_table <- hi_c_table[, c(1, 2, 3, 4, 5, 10)]
  names(hi_c_table) <- c("chr", "start", "end", "DpnII", "HindIII", "Id")
  mergedTable <- merge(hi_c_table, asymmetry_table, by = "Id")
  return(mergedTable)
}

make_box_plot <- function(asymmetry_compartment_table, enzyme, symbols, genomic_element)
{
  "This function gets df for plotting and dinucleotide name for the label
  and creates a boxplot graph using ggplot package"

    p <- ggplot(asymmetry_compartment_table, aes(x = Compartment, y = avg_asymmetry_score, fill = Compartment)) + #Define the elements for plotting - group by "strandness".
    geom_boxplot(outlier.shape = NA, colour = "black") +
    scale_fill_manual(values = c("purple3","palevioletred1")) + theme_classic() +
    coord_cartesian(ylim = c(-0.42, 0.43)) +
    stat_summary(fun=mean, colour="black", geom ="point", shape=18, size=4 ,show.legend = FALSE, position = position_dodge(0.75)) +
    geom_hline(yintercept = 0) +
    theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title.x = element_text(color="black", size=24), axis.title.y = element_text(color="black", size=28, margin = margin(0,20,0,0)),
          axis.text.x = element_text(color="black", size=24), axis.text.y = element_text(color="black", size=24), legend.position="none") +
    labs(x = "Compartment", y = "Average asymmetry score") +
    ggtitle(enzyme) +
    theme(plot.title = element_text(size = 24, hjust = 0.5)) +
      annotate("text", label = symbols, x = 1.5, y = Inf, vjust = 2, size = 6)
  
  pdf(paste(figurePath, paste(genomic_element, enzyme, "compartments_asymmetry_box_plot.pdf", sep = "_"), sep="/"), height = 8, width = 6)  
  print(p)
  dev.off()
}


arrangeTablesUniqGenes <- function(mergedTable, genomic_element, enzyme)
{
  
  mergedTable <- mergedTable[, c("Id", "DpnII", "HindIII", "asymmetryScore")]
  asymmetry_compartment_table <- mergedTable %>%  mutate(Compartment = if_else(get(enzyme) > 0, "A", "B")) # Define the compartment based on A/B compartment score.  
  
  # Scan the genes, in case one gene is related to both, compatment A and B, remove it from the table.
  for(i in 1:nrow(asymmetry_compartment_table)-1)
  {
    if(identical(asymmetry_compartment_table[i,"Id"], asymmetry_compartment_table[i+1,"Id"]) && !(identical(asymmetry_compartment_table[i,"Compartment"], asymmetry_compartment_table[i+1,"Compartment"]))){
    asymmetry_compartment_table <- asymmetry_compartment_table[-c(i,i+1),] # Remove the gene from the table.
    }
  }
  
  asymmetry_compartment_table <- asymmetry_compartment_table[complete.cases(asymmetry_compartment_table),] # Remove rows that contain NA.
  agg_asymmetry_compartment_table <- aggregate(asymmetry_compartment_table, by = list(asymmetry_compartment_table$Id), FUN = mean) # Calculate the average asymmetry score for each gene.
  agg_asymmetry_compartment_table <- agg_asymmetry_compartment_table[,c(1,3,4,5)] # Extract the non NAs columns.
  names(agg_asymmetry_compartment_table) <- c("Id", "DpnII", "HindIII", "avg_asymmetry_score")
  agg_asymmetry_compartment_table <- agg_asymmetry_compartment_table %>%  mutate(Compartment = if_else(get(enzyme) > 0, "A", "B")) # Define the compartment based on A/B compartment score.

  symbols <- adjust_p_value(agg_asymmetry_compartment_table)
  make_box_plot(agg_asymmetry_compartment_table, enzyme, symbols, genomic_element)
 
  
}

asymmetryTableExons <- asymmetry_score(exonsCountingTable)
mergedTableExons <- mergeTables(intersectedFile, asymmetryTableExons)

asymmetryTableIntrons <- asymmetry_score(intronsCountingTable)
mergedTableIntrons <- mergeTables(intersectedFile, asymmetryTableIntrons)

arrangeTablesUniqGenes(mergedTableExons, "uniq_exons", "DpnII")
arrangeTablesUniqGenes(mergedTableExons, "uniq_exons", "HindIII")

arrangeTablesUniqGenes(mergedTableIntrons, "uniq_introns", "DpnII")
arrangeTablesUniqGenes(mergedTableIntrons, "uniq_introns", "HindIII")
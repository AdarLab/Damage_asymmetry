#!/usr/bin/Rscript 
# R version 3.5.2
library(ggplot2)
library(reshape2)
library(ggpubr)
library(plyr)
###########################################################################
# This script gets C->T mutation coverge data of Melanoma across enhancers
# and generate bar plot representing the mutation counts for each strand.
###########################################################################

path <- "mutationsCoverage"
figure_path <- "barPlots"

adjust_p_value <- function(data, alternate)
{
  "This function gets data of 2 populations and compare them by using Wilcoxon test according to the givven alternative hypothesis"
  customSymnum <- list(cutpoints = c(0, c(0.0001, 0.001, 0.01, 0.05)/(nrow(data)/2), 1.01), symbols = c("**", "*****", "**", "*", "n.s.")) # Define the p values symbols based on Bonferroni correction.
  
  transcribed_data <- subset(data, group == "TS")$value # Define the group of the transcribed strand.
  non_transcribed_data <- subset(data, group == "NTS")$value # "" non transcribed strand.
  wilcoxon_stats <- wilcox.test(transcribed_data, non_transcribed_data, paired = T, alternative = alternate)
  p_value <- wilcoxon_stats$p.value
  print(wilcoxon_stats)
  # Scan the customSymnum list in order to adjust the matching symbol to the calculated p value.
  for(i in 1:(length(customSymnum$cutpoints)-1)){
    if(customSymnum$cutpoints[i] <= p_value && p_value < customSymnum$cutpoints[i+1]){
      symbol <- customSymnum$symbols[i]
    }
  }
  return(symbol)
}

generatePlotData <- function(transcribedTable, nonTranscribedTable, elementName)
{
  "This function gets the element counting table,
  the element name and the specific dinucleotide for plotting.
  and returns data frame of all the relevent data for plotting"
  
  a <- data.frame(group = "TS", value = transcribedTable[,11], element = rep(elementName, length(transcribedTable[,8]))) 
  b <- data.frame(group = "NTS", value = nonTranscribedTable[,11] ,element = rep(elementName, length(nonTranscribedTable[,8])))
 
  Data <- rbind(a, b) # Concatnate the transcribed and non-transcribed data frames into one new df.

  return(Data)
}

makeBarPlot <- function(plot.data, symbols, element, cancerType)
{
  "This function gets df for plotting and dinucleotide name for the label
  and creates a boxplot graph using ggplot package"

    p <- ggplot(plot.data, aes(x = group, y = value, fill = group)) + # Define the elements for plotting - group by "strandness".
    geom_bar(stat="identity") +
    scale_fill_manual(values = c("goldenrod","darkgreen")) + 
    scale_color_manual(values = c("goldenrod","darkgreen")) +
    theme_classic() +
    labs(x = "", y = "") +
    theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title.x = element_text(color="black", size=24), axis.title.y = element_text(color="black", size=28, margin = margin(0,20,0,0)),
          axis.text.x = element_text(color="black", size=24), axis.text.y = element_text(color="black", size=24), legend.position="none") +
    labs(y = "C>T mutations counts (per Kb)") +
      annotate("text", label = symbols, x = 1.5, y = Inf, vjust = 2, size = 6) +
      ggtitle(cancerType)
  
  pdf(paste(figure_path, paste(element, cancerType, "C_T_mutations_bar_plot.pdf", sep = "_"), sep="/"), height = 8, width = 6)  
  print(p)
  dev.off()
}


for (file in Sys.glob(paste0(path, "/transcribedStrand/*hg19*"))) {
  file_name <- basename(file)
  cancerType <- gsub("hg19_enhancers_CtoT_", "", file_name)
  cancerType <- gsub("_stranded.bed","",cancerType)
  
  TableEnhancersTranscribed <- read.table(file = paste(path, "transcribedStrand", file_name ,sep = '/'), sep = '\t',header=F)  
  TableEnhancersNonTranscribed <- read.table(file = paste(path, "nonTranscribedStrand", file_name, sep = '/'), sep = '\t',header=F)
  element.data <- generatePlotData(TableEnhancersTranscribed, TableEnhancersNonTranscribed, "Enhancers")
  symbols <- c(symbols, adjust_p_value(element.data, "less"))
  makeBarPlot(element.data, symbols, elements[i], cancerType)  
  
}
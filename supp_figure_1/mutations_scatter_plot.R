#!/usr/bin/Rscript 
# R version 3.5.2
library(ggplot2)
library(reshape2)
library(ggpubr)
library(plyr)
################################################################################################################################
# This script gets C->T mutation coverge data of different cancer types across enhancers
# and compare the strand asymmetry score vs the siginficance of the comparison of the mutation frequencies between the strands.
################################################################################################################################

path <- "mutationsCoverage"
figure_path <- "figures"

generatePlotData <- function(transcribedTable, nonTranscribedTable, elementName)
{
  "This function gets the element counting table,
  the element name and the specific dinucleotide for plotting.
  and returns data frame of all the relevent data for plotting"
  
  a <- data.frame(group = "TS", value = transcribedTable[,11], element = rep(elementName, length(transcribedTable[,8]))) 
  b <- data.frame(group = "NTS", value = nonTranscribedTable[,11] ,element = rep(elementName, length(nonTranscribedTable[,8])))
  Data <- rbind(a, b) # Concatenate the transcribed and non-transcribed data frames into one new df.
  return(Data)
}

adjust_p_value <- function(data, alternate)
{
  "This function gets data of 2 populations and compare them by using Wilcoxon test according to the givven alternative hypothesis"
  customSymnum <- list(cutpoints = c(0, c(0.0001, 0.001, 0.01, 0.05)/(nrow(data)/2), 1.01), symbols = c("**", "***", "**", "*", "n.s.")) # Define the p values symbols based on Bonferroni correction.
  
  transcribed_data <- subset(data, group == "TS")$value # Define the group of the transcribed strand.
  non_transcribed_data <- subset(data, group == "NTS")$value # "" non transcribed strand.
  wilcoxon_stats <- wilcox.test(transcribed_data, non_transcribed_data, paired = T, alternative = alternate)
  p_value <- wilcoxon_stats$p.value
  #print(wilcoxon_stats)
  # Scan the customSymnum list in order to adjust the matching symbol to the calculated p value.
  for(i in 1:(length(customSymnum$cutpoints)-1)){
    if(customSymnum$cutpoints[i] <= p_value && p_value < customSymnum$cutpoints[i+1]){
      symbol <- customSymnum$symbols[i]
    }
  }
  p_adj <- p_value*(nrow(data)/2) # Adjust the p-value according to Bonferroni correction.

  return(-log(p_adj)) # log10
}


asymmetry_score <- function(data)
{
  transcribed_data <- subset(data, group == "TS")$value # Define the group of the transcribed strand.
  non_transcribed_data <- subset(data, group == "NTS")$value # "" non transcribed strand.
  asymmetry_score <- (transcribed_data - non_transcribed_data)/(transcribed_data + non_transcribed_data)
  avg_asymmetry_score <- mean(na.omit(asymmetry_score)) # Extract the average asymmetry score.
  return(avg_asymmetry_score)
}

makeScatterPlot <- function(plot.data, figurePath)
{
  "This function plot the average asymetry score across enhancers vs the
   adjusted p-value for comparison of C->T mutation counts between the strands."
  p <- ggplot(plot.data, aes(x = asymmetry_score, y = log_p_adj , colour = log_p_adj > 1.3, label = cancer)) + geom_point(size = 2.5, alpha = 0.8) + #xlim(0,axis.lim) + ylim(0,axis.lim) +
    geom_hline(yintercept=1.3, linetype="dashed") + scale_colour_manual(name = 'log_p_adj > 1.3', values = setNames(c('#9F1B32','gray70'),c(T, F))) +
    theme(legend.title=element_blank(),
          legend.text = element_text(size=16),
          axis.text.x = element_text(color = "black", size = 36),
          axis.text.y = element_text(color = "black", size = 36),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          legend.key = element_rect(colour = "transparent", fill = "white")) +
    theme(plot.title = element_text(lineheight=.8, hjust = 0.5, size = 20),
          axis.title.y = element_text(size = 28, angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.x = element_text(size = 28, angle = 0, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
    labs(y = "q value(-log10)" , x = "Average asymmetry score") +
    theme(legend.position = "none") +
    geom_text(aes(label = cancer),hjust=0.5, vjust=-0.5, size = 14)
  pdf(paste(figurePath, "mutations_scatter_plot.pdf", sep = "/"), height = 18, width = 20)
  print(p)
  dev.off()

}

plot_table <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(plot_table) <- c("cancer", "asymmetry_score", "log_p_adj")

for (file in Sys.glob(paste0(path, "/transcribedStrand/*hg19*"))) {
  file_name <- basename(file)
  cancerType <- gsub("hg19_enhancers_CtoT_", "", file_name)
  cancerType <- gsub("_stranded.bed","",cancerType)
  print(cancerType)
  
  TableEnhancersTranscribed <- read.table(file = paste(path, "transcribedStrand", file_name ,sep = '/'), sep = '\t',header=F)  
  TableEnhancersNonTranscribed <- read.table(file = paste(path, "nonTranscribedStrand", file_name, sep = '/'), sep = '\t',header=F)

  element.data <- generatePlotData(TableEnhancersTranscribed, TableEnhancersNonTranscribed, "Enhancers")
  log_p_adj <- adjust_p_value(element.data, "less")
  asymmetryScore <- asymmetry_score(element.data)
  plot_table <- rbind(plot_table, data.frame(cancer = cancerType, asymmetry_score = asymmetryScore, log_p_adj = log_p_adj))
}

makeScatterPlot(plot_table, figure_path)

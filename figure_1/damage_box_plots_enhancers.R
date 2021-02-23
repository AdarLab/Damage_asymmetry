#!/usr/bin/Rscript 
# R version 3.5.2
library(ggplot2)
library(reshape2)
library(ggpubr)

########################################################################################
# This script gets CPD counts across human enhancers and generate boxplots representing
# the distribution of CPDs over the transcribd and non-transcribed strands.
########################################################################################

path <- "your_path"
figure_path <- "figures"

TableEnhancersTranscribed <- read.table(file = paste(path, "transcribedStrand", "enhancers_hg38_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed" ,sep = '/'), sep = '\t',header=F)
TableEnhancersNonTranscribed <- read.table(file = paste(path, "nonTranscribedStrand", "enhancers_hg38_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed" ,sep = '/'), sep = '\t',header=F)

adjust_p_value <- function(data, alternate)
{
  "This function gets data of 2 populations and compare them by using Wilcoxon test according to the givven alternative hypothesis"
  customSymnum <- list(cutpoints = c(0, c(0.0001, 0.001, 0.01, 0.05)/(nrow(data)/2), 1.01), symbols = c("**", "***", "*", "****", "ns")) # Define the p values symbols based on Bonferroni correction.
  
  transcribed_data <- subset(data, group == "TS")$value # Define the group of the transcribed strand.
  non_transcribed_data <- subset(data, group == "NTS")$value # "" non transcribed strand.
  wilcoxon_stats <- wilcox.test(transcribed_data, non_transcribed_data, paired = T, alternative = alternate)
  p_value <- wilcoxon_stats$p.value
  # Scan the customSymnum list in order to adjust the matching symbol to the calculated p value.
  for(i in 1:(length(customSymnum$cutpoints)-1)){
    if(customSymnum$cutpoints[i] <= p_value && p_value < customSymnum$cutpoints[i+1]){
      symbol <- customSymnum$symbols[i]
    }
  }
  print(wilcoxon_stats)
  return(symbol)
}

generatePlotData <- function(transcribedTable, nonTranscribedTable, elementName)
{
  "This function gets the element counting table,
  the element name and the specific dinucleotide for plotting.
  and returns data frame of all the relevent data for plotting"
  a <- data.frame(group = "TS", value = transcribedTable[,7]/transcribedTable[,9]*10^3, element = rep(elementName, length(transcribedTable[,10]))) #Set new data frame with given columns names.
  b <- data.frame(group = "NTS", value = nonTranscribedTable[,7]/nonTranscribedTable[,9]*10^3, element = rep(elementName, length(nonTranscribedTable[,10])))
  Data <- rbind(a, b) # Concatenate the transcribed and non-transcribed data frames into one new df.
  return(Data)
}

makeBoxPlot <- function(plot.data, symbols)
{
  "This function gets df for plotting and dinucleotide name for the label
  and creates a boxplot graph using ggplot package"

    p <- ggplot(plot.data, aes(x = group, y = value, fill = group)) + #Define the elements for plotting - group by "strandness".
    geom_boxplot(outlier.shape = NA, colour = "black") +
    scale_fill_manual(values = c("goldenrod","darkgreen")) + theme_classic() +
    coord_cartesian(ylim = c(0, yLim)) +
    stat_summary(fun.y = mean, colour="black", geom ="point", shape=18, size=4 ,show.legend = FALSE, position = position_dodge(0.75)) +
    labs(x = "", y = "") +
    theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title.x = element_text(color="black", size=24), axis.title.y = element_text(color="black", size=28, margin = margin(0,20,0,0)),
          axis.text.x = element_text(color="black", size=24), axis.text.y = element_text(color="black", size=24), legend.position="none") +
    labs(y = paste("CPD Damage-seq read counts\n(per Kb)")) +
    annotate("text", label = symbols, x = 1.5, y = Inf, vjust = 2, size = 6)
  
  pdf(paste(figure_path, paste0("enhancers_eTSS_500.pdf"), sep="/"), height = 8, width = 6)  
  print(p)
  dev.off()
}


element.data <- generatePlotData(TableEnhancersTranscribed, TableEnhancersNonTranscribed, "Enhancers")
symbols <- adjust_p_value(element.data, "less")
makeBoxPlot(element.data, symbols)
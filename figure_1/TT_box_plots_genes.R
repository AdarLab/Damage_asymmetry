#!/usr/bin/Rscript 
# R version 3.5.2
######################################################################################################
# This script gets human genes counting tables of dinucleodites frequency and generate boxplots
# representing the distribution of UV terget sequences on the transcribd and non-transcribed strands.
######################################################################################################
library(ggplot2)
library(reshape2)
library(ggpubr)

# Elements tables
path <- "your_path"
transDinuc <- "TT"
nonTransDinuc <- "AA"
figure_path <- "your_figures_path"


TableGenes <- read.csv(paste(path, "genes_seq_counting.csv", sep = "/"), row.names = 1, header = T)

adjust_p_value <- function(data, alternate)
{
  "This function gets data of 2 populations and compare them by using Wilcoxon test according to the givven alternative hypothesis"
  customSymnum <- list(cutpoints = c(0, c(0.0001, 0.001, 0.01, 0.05)/(nrow(data)/2), 1.01), symbols = c("***", "**", "*", "?", "n.s.")) # Define the p values symbols based on Bonferroni correction.
  
  transcribed_data <- subset(data, group == "TS")$value # Define the group of the transcribed strand.
  non_transcribed_data <- subset(data, group == "NTS")$value # "" non transcribed strand.
  stats <- wilcox.test(transcribed_data, non_transcribed_data, paired = T, alternative = alternate)
  p_value <- stats$p.value
  # Scan the customSymnum list in order to adjust the matching symbol to the calculated p value.
  for(i in 1:(length(customSymnum$cutpoints)-1)){
    if(customSymnum$cutpoints[i] <= p_value && p_value < customSymnum$cutpoints[i+1]){
      symbol <- customSymnum$symbols[i]
    }
  }
  return(symbol)
}

generatePlotData <- function(elementTable, elementName, transDinuc, nonTransDinuc) {
  "This function gets the element counting table,
  the element name and the specific dinucleotide for plotting.
  and returns data frame of all the relevent data for plotting"
  a <- data.frame(group = "TS", value = elementTable[[transDinuc]], element = rep(elementName, length(elementTable[transDinuc]))) 
  b <- data.frame(group = "NTS", value = elementTable[[nonTransDinuc]], element = rep(elementName, length(elementTable[transDinuc])))
  Data <- rbind(a, b) # Concatenate the transcribed and non-transcribed data frames into one new df.
  return(Data)
}
makeBoxPlot <- function(plot.data, dinuc, symbols, element) {
  "This function gets df for plotting and dinucleotide name for the label
  and creates a boxplot graph using ggplot package"
  p <- ggplot(plot.data, aes(x = group, y = value, fill = group)) + # Define the elements for plotting - group by "strandness".
    geom_boxplot(outlier.shape = NA, colour = "black") +
    scale_fill_manual(values = c("goldenrod", "darkgreen")) + theme_classic() +
    stat_summary(fun = mean, colour = "black", geom = "point", shape = 18, size = 4, show.legend = FALSE, position = position_dodge(0.75)) +
    theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title.x = element_text(color = "black", size = 28), axis.title.y = element_text(color = "black", size = 28, margin = margin(0, 20, 0, 0)),
      axis.text.x = element_text(color = "black", size = 24), axis.text.y = element_text(color = "black", size = 24), legend.position = "none") +
    labs(y = paste(dinuc, "frequency (per Kb)", sep = " "), x = "") +
    annotate("text", label = symbols, x = 1.5, y = Inf, vjust = 2, size = 6)
   pdf(paste(figure_path,(paste(element, dinuc, "frequency.pdf", sep = "_")), sep = "/"), height = 8, width = 6)
  print(p)
   dev.off()
}


TableGenes[,-c(1,2)] <- (TableGenes[,-c(1,2)]/TableGenes$length)*10^3 # Normalize the sequence counts to 1 Kb.
element.data <- generatePlotData(TableGenes, "Genes", nonTransDinuc, transDinuc)
symbols <- c(symbols, adjust_p_value(element.data, "less"))
makeBoxPlot(element.data, transDinuc, symbols, "Genes")
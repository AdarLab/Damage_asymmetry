#!/usr/bin/Rscript 
library(ggplot2)
library(reshape2)
library(ggpubr)

path <- "ATAC_seq_data/coverageData"

figure_path <- "figures"

tableExons_NO_UV <- read.table(file = paste(path, "exons", "ATAC_seq_NO_UV_chopped_exons_no_first_last.bed", sep = "/"), sep = '\t',header = F)
tableIntrons_NO_UV <- read.table(file = paste(path, "introns", "ATAC_seq_NO_UV_chopped_introns_no_first_last.bed", sep = "/"), sep = '\t',header = F)

tableExons_2_h <- read.table(file = paste(path, "exons", "ATAC_seq_2h_15J_chopped_exons_no_first_last.bed", sep = "/"), sep = '\t',header = F)
tableIntrons_2_h <- read.table(file = paste(path, "introns", "ATAC_seq_2h_15J_chopped_introns_no_first_last.bed", sep = "/"), sep = '\t',header = F)


adjust_p_value <- function(data, alternate)
{
  "This function gets data of 2 populations and compare them by using Wilcoxon test according to the givven alternative hypothesis"
  customSymnum <- list(cutpoints = c(0, c(0.0001, 0.001, 0.01, 0.05)/(nrow(data)/2), 1.01), symbols = c("**", "***", "****", "*", "n.s.")) # Define the p values symbols based on Bonferroni correction.
  
  introns <- subset(data, group == "Introns")$value # Define the group of the transcribed strand.
  exons <- subset(data, group == "Exons")$value # "" non transcribed strand.
  wilcoxon_stats <- wilcox.test(introns, exons, paired = F, alternative = alternate)
  p_value <- wilcoxon_stats$p.value
  # Scan the customSymnum list in order to adjust the matching symbol to the calculated p value.
  for(i in 1:(length(customSymnum$cutpoints)-1)){
    if(customSymnum$cutpoints[i] <= p_value && p_value < customSymnum$cutpoints[i+1]){
      symbol <- customSymnum$symbols[i]
    }
  }
  return(symbol)
}

generatePlotData <- function(tables_introns, tables_exons, condition)
{
  "This function gets the element counting table,
  the element name and the specific dinucleotide for plotting.
  and returns data frame of all the relevent data for plotting"
  
  a <- data.frame(group = "Introns", value = tables_introns[,7]/(tables_introns[,3] - tables_introns[,2])*1000) # Normalize the counts to 1Kb.
  b <- data.frame(group = "Exons", value = tables_exons[,7]/(tables_exons[,3] - tables_exons[,2])*1000)
  Data <- rbind(a, b) # Concatnate the exons and introns data frames into one new df.
  return(Data)
}

makeBoxPlot <- function(plot.data, symbols, condition)
{
  "This function gets df for plotting and dinucleotide name for the label
  and creates a boxplot graph using ggplot package"

    p <- ggplot(plot.data, aes(x = group, y = value, fill = group)) + # Define the elements for plotting - group by element type.
    geom_boxplot(outlier.shape = NA, colour = "black") +
    scale_fill_manual(values = c("gray48","mediumpurple3")) + theme_classic() +
    coord_cartesian(ylim = c(0, 175)) +
    stat_summary(fun=mean, colour="black", geom ="point", shape=18, size=4 ,show.legend = FALSE, position = position_dodge(0.75)) +
    labs(x = "", y = "") +
    theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title.x = element_text(color="black", size=24), axis.title.y = element_text(color="black", size=28, margin = margin(0,20,0,0)),
          axis.text.x = element_text(color="black", size=24), axis.text.y = element_text(color="black", size=24), legend.position="none") +
    labs(y = "ATAC-seq read counts (per Kb)") +
      annotate("text", label = symbols, x = 1.5, y = Inf, vjust = 2, size = 6)
  
  pdf(paste(figure_path, paste("Introns_exons", condition, "ATAC_seq.pdf", sep = "_"), sep="/"), height = 8, width = 6.5)  
  print(p)
  dev.off()
}

tables_exons <- c("tableExons_NO_UV", "tableExons_2_h")
tables_introns <- c("tableIntrons_NO_UV", "tableIntrons_2_h")
conditions <- c("no_UV", "15J_2h")
plot.data <- data.frame() #Intialize data frame.
symbols <- c() # A vector for the p-values symbols
#Scans all the tables, update them and generate df for the boxplot function.
for(i in 1:length(tables_exons))
{
  element.data <- generatePlotData(get(tables_introns[i]), get(tables_exons[i]), conditions[i])
  symbols <- c(symbols, adjust_p_value(element.data, "less"))
  makeBoxPlot(element.data, symbols, conditions[i])
  #plot.data <- rbind(plot.data, element.data ) #Concate the current element's data. 
}
#makeBoxPlot(plot.data, symbols)

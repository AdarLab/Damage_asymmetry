#!/usr/bin/Rscript 
# R version 3.5.2
library(ggplot2)
library(reshape2)
library(ggpubr)

###########################################################################################################################################
# This script gets XR-Seq data across enhancers, normalized it to TT counts and check the differences in repair levels between the strands.
###########################################################################################################################################

figure_path <- "figures"
covergePath <- "coverageData"
repair_path <- paste(covergePath, "repair", sep = "/")
TT_path <- paste(covergePath, "tables/counting", sep = "/")

TT_file_name <- "hg_38_enhancers_seq_counting.csv"

########################
#-----Num of reads-----#
########################

numOfReadsRepairNHF1 <- 11100235 
numOfReadsRepairXPC <- 14856828 
numOfReadsRepairCSB <- 19152949 
numOfReadsDamage <- 68855889
numOfReadsTT <- 576715193 


generatePlotDataRatioTT <- function(enhancersTableTT, enhancersTableRepairTranscribed, enhancersTableRepairNonTranscribed, numOfReadsRepair) {
  "This function gets the element counting table,
  the element name and the specific dinucleotide for plotting.
  and returns data frame of all the relevent data for plotting"
  tableNames <- c("Strand", "Ratio")
  
  # Extract the TT frequency.
  table_TT_transcribed <- enhancersTableTT[,"AA"] / enhancersTableTT[,"length"]
  table_TT_non_transcribed <- enhancersTableTT[,"TT"] / enhancersTableTT[,"length"]
  
  # Normlaize the repair counts to reads per 10^6
  table_repair_transcribed <- (enhancersTableRepairTranscribed[, 7] / enhancersTableRepairTranscribed[, 9])/numOfReadsRepair*10^6
  table_repair_non_transcribed <- (enhancersTableRepairNonTranscribed[, 7] / enhancersTableRepairNonTranscribed[, 9])/numOfReadsRepair*10^6
  
  # Normalize the repair levels to TT frequency.
  table_transcribed_ratio <- as.data.frame(table_repair_transcribed / table_TT_transcribed)
  table_non_transcribed_ratio <- as.data.frame(table_repair_non_transcribed / table_TT_non_transcribed)
  
  table_transcribed_plot <- cbind(as.data.frame(rep("TS", nrow(table_transcribed_ratio))), table_transcribed_ratio) # Define column of transcribed labels
  table_non_transcribed_plot <- cbind(as.data.frame(rep("NTS", nrow(table_non_transcribed_ratio))), table_non_transcribed_ratio) # Define column of transcribed labels
  
  # Set the names of the 2 dataframes to an identical name for enabling rbind function.
  names(table_transcribed_plot) <- tableNames
  names(table_non_transcribed_plot) <- tableNames
  tableForPlot <- rbind(table_transcribed_plot, table_non_transcribed_plot)
  names(tableForPlot) <- tableNames
  
  return(tableForPlot)
}

makeBoxPlot <- function(plot.data, title, figureName, yLab, yLim, stars_pos) {

  "This function gets df for plotting and dinucleotide name for the label
  and creates a boxplot graph using ggplot package"

  # yLim = boxplot.stats(plot.data$value)$stats[c(1, 5)]
  customSymnum <- list(cutpoints = c(0, c(0.0001, 0.001, 0.01, 0.05) / (nrow(plot.data) / 2), 1.01), symbols = c("***", "****", "**", "*?", "n.s.")) # Define the p values symbols based on Bonferroni correction.
  p <- ggplot(plot.data, aes(x = Strand, y = value, fill = Strand)) + # Define the elements for plotting - group by "strandness".
    stat_compare_means(aes(group = Strand), label = "p.signif", method.args = list(paired = T, alternative = "less"), label.x = 1.5, label.y = stars_pos, symnum.args = customSymnum, size = 10) +
    geom_boxplot(outlier.shape = NA, colour = "black") +
    scale_fill_manual(values = c("red", "blue")) + theme_classic() +
    coord_cartesian(ylim = yLim) +
    stat_summary(fun = mean, colour = "black", geom = "point", shape = 18, size = 4, show.legend = FALSE, position = position_dodge(0.75)) +
    theme(
      plot.title = element_text(size = 22, hjust = 0.5), axis.title.x = element_text(color = "black", size = 28), axis.title.y = element_text(color = "black", size = 28, margin = margin(0, 20, 0, 0)),
      axis.text.x = element_text(color = "black", size = 24), axis.text.y = element_text(color = "black", size = 24), legend.position = "none"
    ) +
    labs(y = yLab, x = "") +
    ggtitle(title)
  pdf(paste(figure_path, figureName, sep = "/"), height = 8, width = 6)
  print(p)
  dev.off()
}


make_box_plot_repair_TT <- function(repair_path, numOfReadsRepair, cellType)
{
  repairFileName <- paste("enhancers_hg38", cellType, "CPD_4nt.bed", sep = "_")
  enhancers_repair_transcribed <- read.table(paste(repair_path, "transcribedStrand", repairFileName, sep = "/"), sep = "\t", header = F)
  enhancers_repair_non_transcribed <- read.table(paste(repair_path, "nonTranscribedStrand", repairFileName, sep = "/"), sep = "\t", header = F)
  
  
  damage_repair_both_strands_ratio_data <- generatePlotDataRatioTT(enhancers_TT_table, enhancers_repair_transcribed, enhancers_repair_non_transcribed, numOfReadsRepair)
  damage_repair_both_strands_ratio_data_plot <- melt(damage_repair_both_strands_ratio_data, id.var = "Strand")
  makeBoxPlot(damage_repair_both_strands_ratio_data_plot, cellType, paste0(cellType, "_repair_TT.pdf"),  "XR-seq/TT counts", c(0,0.017), 0.016)
  
}

##################
#---- TT ----#
##################

enhancers_TT_table <- read.csv(paste(TT_path, TT_file_name, sep = "/"), header = T)

make_box_plot_repair_TT(repair_path, numOfReadsRepairNHF1, "NHF1")
make_box_plot_repair_TT(repair_path, numOfReadsRepairXPC, "XPC")
make_box_plot_repair_TT(repair_path, numOfReadsRepairCSB, "CSB")
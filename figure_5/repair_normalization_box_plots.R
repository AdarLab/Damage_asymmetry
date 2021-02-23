#!/usr/bin/Rscript 
# R version 3.5.2

library(ggplot2)
library(reshape2)
library(ggpubr)

########################################################################################################################################
# This script gets XR-Seq and Damage-Seq data and check the differences in repair levels normalized to TT frequency in exons and intons.
########################################################################################################################################

figure_path <- "figures"
repair_path <- "XR_Seq_data"
TT_path <- "tables/counting"

exons_repair_file_name <- "CPD_hg38_RefSeq_uniq_chopped_exons_cov.bed"
introns_repair_file_name <- "CPD_hg38_RefSeq_uniq_chopped_introns_cov.bed"

exons_TT_file_name <- "chopped_exons_seq_counting.csv"
introns_TT_file_name <- "chopped_introns_seq_counting.csv"

numOfReadsRepairNHF1 <- 11100235 
numOfReadsRepairXPC <- 14856828 
numOfReadsRepairCSB <- 19152949 

generatePlotDataRatioTT <- function(exonsTableTT, exonsTableRepair, intronTableTT, intronTableRepair, numOfReadsRepair, dinuc) {
  "This function gets the element counting table,
  the element name and the specific dinucleotide for plotting.
  and returns data frame of all the relevent data for plotting"
  tableNames <- c("Type", "Ratio")
  if(dinuc == "Both strands")
  {
    exons_table_TT <- ((exonsTableTT[, "AA"] +  exonsTableTT[, "TT"]) / exonsTableTT[, "length"]) # Define column of transcribed labels
    introns_table_TT <- ((intronTableTT[, "AA"] + intronTableTT[, "TT"]) / intronTableTT[, "length"]) # Define column of transcribed labels
  } else
  {  
    exons_table_TT <- ((exonsTableTT[, dinuc]) / exonsTableTT[, "length"]) # Define column of transcribed labels
    introns_table_TT <- ((intronTableTT[, dinuc]) / intronTableTT[, "length"]) # Define column of transcribed labels
  }
  
  exons_table_repair <- (exonsTableRepair[, 7] / (exonsTableRepair[, 3] - exonsTableRepair[, 2])) /numOfReadsRepair * 10^6# Define column of transcribed labels
  introns_table_repair <- (intronTableRepair[, 7] / (intronTableRepair[, 3] - intronTableRepair[, 2])) /numOfReadsRepair * 10^6  # Define column of transcribed labels

  exons_ratio <- as.data.frame(exons_table_repair / exons_table_TT)
  introns_ratio <- as.data.frame(introns_table_repair / introns_table_TT)

  exons_table_plot <- cbind(as.data.frame(rep("Exons", nrow(exonsTableTT))), exons_ratio) # Define column of transcribed labels
  introns_table_plot <- cbind(as.data.frame(rep("Introns", nrow(intronTableTT))), introns_ratio) # Define column of transcribed labels

  # Set the names of the 2 dataframes to an identical name for enabling rbind function.
  names(exons_table_plot) <- tableNames
  names(introns_table_plot) <- tableNames
  tableForPlot <- rbind(introns_table_plot, exons_table_plot)
  names(tableForPlot) <- tableNames

  return(tableForPlot)
}

makeBoxPlot <- function(plot.data, title, figureName, yLab, yLim, stars_pos) {

  "This function gets df for plotting and dinucleotide name for the label
  and creates a boxplot graph using ggplot package"

  # yLim = boxplot.stats(plot.data$value)$stats[c(1, 5)]
  customSymnum <- list(cutpoints = c(0, c(0.0001, 0.001, 0.01, 0.05) / (nrow(plot.data) / 2), 1.01), symbols = c("**", "****", "***", "*", "n.s.")) # Define the p values symbols based on Bonferroni correction.
  p <- ggplot(plot.data, aes(x = Type, y = value, fill = Type)) + # Define the elements for plotting - group by "strandness".
    stat_compare_means(aes(group = Type), label = "p.signif", method.args = list(paired = FALSE, alternative = "less"), label.x = 1.5, label.y = stars_pos, symnum.args = customSymnum, size = 10) +
    geom_boxplot(outlier.shape = NA, colour = "black") +
    scale_fill_manual(values = c("red", "blue")) + theme_classic() +
    coord_cartesian(ylim = c(0,yLim)) +
    stat_summary(fun = mean, colour = "black", geom = "point", shape = 18, size = 4, show.legend = FALSE, position = position_dodge(0.75)) +
    theme(
      plot.title = element_text(size = 22, hjust = 0.5), axis.title.x = element_text(color = "black", size = 28), axis.title.y = element_text(color = "black", size = 28, margin = margin(0, 20, 0, 0)),
      axis.text.x = element_text(color = "black", size = 24), axis.text.y = element_text(color = "black", size = 24), legend.position = "none"
    ) +
    labs(y = yLab, x = "") +
    ggtitle(title)
  pdf(paste(figure_path, figureName, sep = "/"), height = 8, width = 6.5)
  print(p)
  dev.off()
}


make_box_plot_repair_TT <- function(repair_path, numOfReadsDamage, numOfReadsRepair, cellType, yLim)
{
  exons_repair_transcribed <- read.table(paste(repair_path, cellType, "CPD", "transcribedStrand", paste(cellType, exons_repair_file_name, sep = "_"), sep = "/"), sep = "\t", header = F)
  exons_repair_non_transcribed <- read.table(paste(repair_path,  cellType, "CPD", "nonTranscribedStrand",  paste(cellType, exons_repair_file_name, sep = "_"), sep = "/"), sep = "\t", header = F)
  
  introns_repair_transcribed <- read.table(paste(repair_path,  cellType, "CPD", "transcribedStrand",  paste(cellType, introns_repair_file_name, sep = "_"), sep = "/"), sep = "\t", header = F)
  introns_repair_non_transcribed <- read.table(paste(repair_path,  cellType, "CPD", "nonTranscribedStrand",  paste(cellType, introns_repair_file_name, sep = "_"), sep = "/"), sep = "\t", header = F)
  
  # Both strands
  exons_repair <- exons_repair_transcribed
  exons_repair$V7 <- exons_repair$V7 + exons_repair_non_transcribed$V7
  
  introns_repair <- introns_repair_transcribed
  introns_repair$V7 <- introns_repair$V7 + introns_repair_non_transcribed$V7

  
  damage_repair_both_strands_ratio_data <- generatePlotDataRatioTT(exons_TT, exons_repair, introns_TT, introns_repair, numOfReadsRepair, "Both strands")
  damage_repair_both_strands_ratio_data_plot <- melt(damage_repair_both_strands_ratio_data, id.var = "Type")
  
  damage_repair_transcribed_strand_ratio_data <- generatePlotDataRatioTT(exons_TT, exons_repair_transcribed, introns_TT, introns_repair_transcribed, numOfReadsRepair, "AA")
  
  damage_repair_transcribed_strand_ratio_data_plot <- melt(damage_repair_transcribed_strand_ratio_data, id.var = "Type")
  print(head(damage_repair_transcribed_strand_ratio_data_plot))
  makeBoxPlot(damage_repair_transcribed_strand_ratio_data_plot, cellType, paste0(cellType, "_repair_TT_transcribed_strand.pdf"), "XR-seq/TT counts", yLim, yLim-0.002)
  
  damage_repair_non_transcribed_strand_ratio_data <- generatePlotDataRatioTT(exons_TT, exons_repair_non_transcribed, introns_TT, introns_repair_non_transcribed, numOfReadsRepair, "TT")
  damage_repair_non_transcribed_strand_ratio_data_plot <- melt(damage_repair_non_transcribed_strand_ratio_data, id.var = "Type")
  makeBoxPlot(damage_repair_non_transcribed_strand_ratio_data_plot, cellType, paste0(cellType, "_repair_TT_non_transcribed_strand.pdf"), "XR-seq/TT counts", yLim, yLim-0.002)
  
}

exons_TT <- read.table(paste(TT_path, exons_TT_file_name, sep = "/"), sep = ",", header = T)
introns_TT <- read.table(paste(TT_path, introns_TT_file_name, sep = "/"), sep = ",", header = T)

make_box_plot_repair_TT(repair_path, numOfReadsTT, numOfReadsRepairNHF1, "NHF1", 0.035)
make_box_plot_repair_TT(repair_path, numOfReadsTT, numOfReadsRepairXPC, "XPC", 0.08)
make_box_plot_repair_TT(repair_path, numOfReadsTT, numOfReadsRepairCSB, "CSB", 0.018)

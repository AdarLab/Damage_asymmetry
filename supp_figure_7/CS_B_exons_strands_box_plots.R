#############################################################################################################
# This script gets XR-Seq data and TT counts table, normalized the repair to TT frequency
# and generates box plot of normalized repair in CS-B cells (GGR only) on the transcribed and non-transcribed strands of exons.
###############################################################################################################################
library(ggplot2)
library(reshape2)
library(ggpubr)

figure_path <- "figures"
repair_path <- "XR_Seq_data"
TT_path <- "tables/counting"

exons_repair_file_name <- "CPD_hg38_RefSeq_uniq_chopped_exons_cov.bed"
exons_TT_file_name <- "chopped_exons_seq_counting.csv"

numOfReadsRepairCSB <- 19152949 


generatePlotDataRatioTT <- function(tableTT, tableRepairTranscribed, tableRepairNonTranscribed, numOfReadsRepair) {
  "This function gets the element counting table,
  the element name and the specific dinucleotide for plotting.
  and returns data frame of all the relevent data for plotting"
  tableNames <- c("Type", "Ratio")
  
  # Normalize the TT counts to the length.
  table_TT_trans <- as.data.frame(tableTT[, "AA"] / tableTT[, "length"]) 
  table_TT_non_trans <- as.data.frame(tableTT[, "TT"] / tableTT[, "length"])
  
  # Normalize the reads to Mb.
  table_repair_trans <- (tableRepairTranscribed[, 7] / (tableRepairTranscribed[, 3] - tableRepairTranscribed[, 2])) /numOfReadsRepair * 10^6
  table_repair_non <- (tableRepairNonTranscribed[, 7] / (tableRepairNonTranscribed[, 3] - tableRepairNonTranscribed[, 2])) /numOfReadsRepair * 10^6  # Define column of transcribed labels

  transcribed_ratio <- as.data.frame(table_repair_trans / table_TT_trans)
  non_trans_ratio <- as.data.frame(table_repair_non / table_TT_non_trans)
  
  transcribed_table_plot <- cbind(as.data.frame(rep("Transcribed", nrow(table_TT_trans))), transcribed_ratio) # Define column of transcribed labels
  non_trans_table_plot <- cbind(as.data.frame(rep("Non-transcribed", nrow(table_TT_non_trans))), non_trans_ratio) # Define column of transcribed labels

  # Set the names of the 2 dataframes to an identical name for enabling rbind function.
  names(transcribed_table_plot) <- tableNames
  names(non_trans_table_plot) <- tableNames
  tableForPlot <- rbind(transcribed_table_plot, non_trans_table_plot)
  names(tableForPlot) <- tableNames

  return(tableForPlot)
}

makeBoxPlot <- function(plot.data, title, figureName, yLab, yLim, stars_pos) {

  "This function gets df for plotting and dinucleotide name for the label
  and creates a boxplot graph using ggplot package"
  print(head(plot.data))
  # yLim = boxplot.stats(plot.data$value)$stats[c(1, 5)]
  customSymnum <- list(cutpoints = c(0, c(0.0001, 0.001, 0.01, 0.05) / (nrow(plot.data) / 2), 1.01), symbols = c("**", "****", "***", "*?", "ns")) # Define the p values symbols based on Bonferroni correction.
  p <- ggplot(plot.data, aes(x = Type, y = value, fill = Type)) + # Define the elements for plotting - group by "strandness".
    stat_compare_means(comparisons = list(c("Transcribed", "Non-transcribed")), label = "p.signif", method.args = list(paired = TRUE, alternative = "greater"), label.x = 1.5, label.y = stars_pos, symnum.args = customSymnum, size = 10, tip.length = 0) +
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
  
  repair_exons_ratio_data <- generatePlotDataRatioTT(exons_TT, exons_repair_transcribed, exons_repair_non_transcribed, numOfReadsRepair)
  repair_exons_ratio_data_plot <- melt(repair_exons_ratio_data, id.var = "Type")
  

  makeBoxPlot(repair_exons_ratio_data_plot, paste(cellType, "Exons"), paste0(cellType, "_repair_TT_exons.pdf"), "XR-seq/TT counts", yLim, yLim-0.002)
}

exons_TT <- read.table(paste(TT_path, exons_TT_file_name, sep = "/"), sep = ",", header = T)
make_box_plot_repair_TT(repair_path, numOfReadsTT, numOfReadsRepairCSB, "CSB", 0.018)
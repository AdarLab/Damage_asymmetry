#!/usr/bin/Rscript 
# R version 3.5.2
library(dplyr)
library(ggplot2)
library(stringr)
library(ggsignif)
library(ggpubr)
warnings()
####################################################################################
# This script gets expression data and counting tables of exons and introns
# and generates box plot representing the expression values vs the asymmetry score.
####################################################################################
path <- "your_path"
figurePath <- "figures"

adjust_p_value <- function(data, alternate)
{
  "This function gets data of 2 populations and compare them by using Wilcoxon test according to the givven alternative hypothesis"
  customSymnum <- list(cutpoints = c(0, c(0.0001, 0.001, 0.01, 0.05)/(nrow(data)/2), 1.01), symbols = c("*", "***", "**", "****", "ns")) # Define the p values symbols based on Bonferroni correction.
  first_quartile_data <- subset(data, quartile_FPKM == "1")$AA.TT # Define the group of the transcribed strand.
  fourth_quartile_data <- subset(data, quartile_FPKM == "4")$AA.TT # "" non transcribed strand.
  t_stats <- t.test(first_quartile_data, fourth_quartile_data, paired = T, alternative = alternate, exact = F)
  print(t_stats)
  p_value <- t_stats$p.value
  # Scan the customSymnum list in order to adjust the matching symbol to the calculated p value.
  for(i in 1:(length(customSymnum$cutpoints)-1)){
    if(customSymnum$cutpoints[i] <= p_value && p_value < customSymnum$cutpoints[i+1]){
      symbol <- customSymnum$symbols[i]
    }
  }
  return(wilcoxon_stats)
}

makeBoxPlots <- function(ntile_table, yLim, yLabels, p_value, title, figureName, p_anove_pose)
{
  xLabs <- c("Lowest", "Low", "High", "Highest")
  my_comparisons <- list( c("1", "2"), c("2", "3"), c("3", "4"))
  bonferroni_symnum.args <- list(cutpoints = c(0, c(0.0001, 0.001, 0.01, 0.05)/(nrow(ntile_table)/4), 1), symbols = c("****", "***", "**", "*", "ns"))
    boxPlot <- ggplot(ntile_table, aes(x=factor(quartile_FPKM), y = AA.TT, fill = factor(quartile_FPKM))) + geom_boxplot(outlier.shape = NA, colour="black") + 
    stat_compare_means(method = "anova", method.args = list(alternative = "greater"), label.y = p_anove_pose) +
    stat_compare_means(method = "t.test", comparisons = my_comparisons, symnum.args = bonferroni_symnum.args, label = "p.adj", method.args = list(alternative = "greater"), label.y = yLabels, tip.length = 0.01) +
    theme_classic() + 
    theme(plot.title = element_text(size = 28, hjust = 0.5), axis.title.x = element_text(color="black", size=28, vjust = -1), axis.title.y = element_text(color="black", size=28, margin = margin(0,20,0,0)),
          axis.text.x = element_text(color="black", size=20), axis.text.y = element_text(color="black", size=20), legend.position="none") +
    scale_x_discrete(labels = xLabs) +
    labs( x = paste("FPKM quartiles"), y = "TT asymmetry score") +
    stat_summary(fun=mean, geom="point", shape=18, size=4, colour = "black") +
    ggtitle(title) +
    scale_fill_manual(values = c("#E87280", "#2D8A26", "#4693C2", "#6763A6")) +
    coord_cartesian(ylim = yLim)
  pdf(paste(figurePath, paste0(figureName, "_expression_vs_asymmetry.pdf"), sep = "/"), height = 8, width = 10)
  print(boxPlot)
  dev.off()
}



creatDataTestites <- function(ratioTable, expressionTable, yLim, yLabels, title)
{
  ratioTable[["Id"]] <- gsub("\\(.*","",ratioTable[["Id"]])
  ratioTableTT <- ratioTable[,c("Id", "length")]
  
  transcribed_tmp <- as.data.frame((ratioTable[,"AA"]/ratioTable$length)*1000) # Etract the expected dinu frequencies on the transcribed strand.
  non_transcribed_tmp <- as.data.frame((ratioTable[,"TT"]/ratioTable$length)*1000) # "" non transcribed.
  # Divide the frequencies in the sum of the dinucleotide on both strands.
  asymmetry <- (transcribed_tmp - non_transcribed_tmp)/(transcribed_tmp + non_transcribed_tmp)
  ratioTableTT <- cbind(ratioTableTT, asymmetry)
  names(ratioTableTT) <- c("Id", "length", "AA.TT")
  
  # ratioTableTT$AA.TT <- log2(ratioTableTT$AA.TT)
  expressionTable <- expressionTable[,-1] # Remove the first columns - Genes symblos
  names(expressionTable) <- c("Id", "FPKM")
  mergedTable <- merge(ratioTableTT, expressionTable, by = "Id" )
  
  mergedTable$quartile_FPKM <- ntile(mergedTable["FPKM"], 4)

  figureName <- paste("testites", title, sep = "_")
  makeBoxPlots(mergedTable, yLim, yLabels, "", title, figureName)
  
  
  return(mergedTable)
}


##########################
#-------- Testites---------#
##########################

testitesRatioTableExons <- read.csv(paste(path, "germCells/sperm/GTEx/exons/tables/counting/continuous_uniq_coding_exons_counting.csv", sep = "/"))
testitesRatioTableIntrons <- read.csv(paste(path, "germCells/sperm/GTEx/introns/tables/counting/continuous_uniq_coding_introns_counting.csv", sep = "/"))
testitesExpressionTable <- read.table(paste(path, "germCells/sperm/GTEx/GTEx_Analysis_gene_median_tpm_Testis_joined.txt", sep = "/"), sep = "\t", header = F)

testites_stats <- creatDataTestites(testitesRatioTableExons, testitesExpressionTable, c(-0.37,0.58), c(0.48, 0.52, 0.56), "Exons")
testites_stats <- creatDataTestites(testitesRatioTableIntrons, testitesExpressionTable, c(-0.37,0.33), c(0.23, 0.27, 0.31), "Introns")


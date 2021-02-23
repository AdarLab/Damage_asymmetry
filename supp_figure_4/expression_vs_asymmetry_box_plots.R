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
path <- "RNA_seq"
figurePath <- "figures"

makeBoxPlots <- function(ntile_table, yLim, yLabels, title, figureName, p_anova_pose)
{
  xLabs <- c("Lowest", "Low", "High", "Highest")
  my_comparisons <- list( c("1", "2"), c("2", "3"), c("3", "4"))
  bonferroni_symnum.args <- list(cutpoints = c(0, c(0.0001, 0.001, 0.01, 0.05)/(nrow(ntile_table)/4), 1), symbols = c("****", "***", "**", "*", "ns"))
    boxPlot <- ggplot(ntile_table, aes(x=factor(quartile_TPM), y = AA.TT, fill = factor(quartile_TPM))) + geom_boxplot(outlier.shape = NA, colour="black") + 
    stat_compare_means(method = "anova", method.args = list(alternative = "greater"), label.y = p_anova_pose) +
    stat_compare_means(method = "t.test", comparisons = my_comparisons, symnum.args = bonferroni_symnum.args, label = "p.adj", method.args = list(alternative = "greater"), label.y = yLabels, tip.length = 0.01) +
    theme_classic() + 
    theme(plot.title = element_text(size = 28, hjust = 0.5), axis.title.x = element_text(color="black", size=28, vjust = -1), axis.title.y = element_text(color="black", size=28, margin = margin(0,20,0,0)),
          axis.text.x = element_text(color="black", size=20), axis.text.y = element_text(color="black", size=20), legend.position="none") +
    scale_x_discrete(labels = xLabs) +
    labs(x = "TPM quartiles", y = "TT asymmetry score") +
    stat_summary(fun=mean, geom="point", shape=18, size=4, colour = "black") +
    ggtitle(title) +
    scale_fill_manual(values = c("#E87280", "#2D8A26", "#4693C2", "#6763A6")) +
    coord_cartesian(ylim = yLim)
  pdf(paste(figurePath, paste0(figureName, "_expression_vs_asymmetry.pdf"), sep = "/"), height = 8, width = 10)
  print(boxPlot)
  dev.off()
}

creatDataSperm <- function(ratioTable, expressionTable, yLim, yLabels,title)
{
  ratioTable[["Id"]] <- gsub("\\(.*","",ratioTable[["Id"]])
  ratioTableTT <- ratioTable[,c("Id", "length")]
  
  transcribed_tmp <- as.data.frame((ratioTable[,"AA"]/ratioTable$length)*1000) # Etract the expected dinu frequencies on the transcribed strand.
  non_transcribed_tmp <- as.data.frame((ratioTable[,"TT"]/ratioTable$length)*1000) # "" non transcribed.
  # Divide the frequencies in the sum of the dinucleotide on both strands.
  asymmetry <- (transcribed_tmp - non_transcribed_tmp)/(transcribed_tmp + non_transcribed_tmp)
  ratioTableTT <- cbind(ratioTableTT, asymmetry)
  names(ratioTableTT) <- c("Id", "length", "AA.TT")
  
  expressionTable <- expressionTable[,-1] # Remove the first columns - Genes symblos
  names(expressionTable) <- c("Id", "FPKM")

  mergedTable <- merge(ratioTableTT, expressionTable, by = "Id" )
  
  mergedTable$quartile_FPKM <- ntile(mergedTable["FPKM"], 4)
  
  figureName <- paste("sperm", title, sep = "_")
  makeBoxPlots(mergedTable, yLim, yLabels, "", title, figureName)
  
  
  return(mergedTable)
}

creatDataESC <- function(countingTable, expressionTable, yLim, yLabels, title, figureName, p_anova_pose)
{
  countingTableTT <- countingTable[,c("Id", "length")]
  countingTableTT$Id <- sub("\\..*", "", countingTableTT$Id)
  
  transcribed_tmp <- as.data.frame((countingTable[,"AA"]/countingTable$length)*1000) # Extract the dinucleotide frequencies on the transcribed strand.
  non_transcribed_tmp <- as.data.frame((countingTable[,"TT"]/countingTable$length)*1000) # "" non transcribed.
  # Divide the frequencies in the sum of the dinucleotide on both strands.
  asymmetry <- (transcribed_tmp - non_transcribed_tmp)/(transcribed_tmp + non_transcribed_tmp)
  countingTableTT <- cbind(countingTableTT, asymmetry)
  names(countingTableTT) <- c("Id", "length", "AA.TT")
  
  # ratioTableTT$AA.TT <- log2(ratioTableTT$AA.TT)
  expressionTableAvg <- cbind(as.data.frame(expressionTable[,"V1"]), (expressionTable[,"V8"]+expressionTable[,"V9"])/2) # Remove the first columns - Genes symblos
  names(expressionTableAvg) <- c("Id", "TPM")
  mergedTable <- merge(countingTableTT, expressionTableAvg, by = "Id" )
  mergedTable <- unique(mergedTable)  
  mergedTable$quartile_TPM <- ntile(mergedTable["TPM"], 4) # Define quartiles based on the TPM values.

  figureName <- paste("H9ESC", title, sep = "_")
  makeBoxPlots(mergedTable, yLim, yLabels, title, figureName, p_anova_pose)
  
  return(mergedTable)
}

##########################
#-------- Sperm ---------#
##########################

spermTableExons <- read.csv(paste(path, "continuous_uniq_coding_exons_counting.csv", sep = "/"))
spermTableIntrons <- read.csv(paste(path, "continuous_uniq_coding_introns_counting.csv", sep = "/"))
spermExpressionTable <- read.table(paste(path, "PGC_bulk_avg_FPKM_joined.tsv", sep = "/"), sep = "\t", header = F)

sperm_stats <- creatDataSperm(spermTableExons, spermExpressionTable, c(-0.45,0.6), c(0.58, 0.54, 0.5), "Exons")
sperm_stats <- creatDataSperm(spermTableIntrons, spermExpressionTable, c(-0.35,0.28), c(0.26, 0.22, 0.18), "Introns")

##########################
#-------- H9 ESCs -------#
##########################
tableExons <- read.csv(paste(path, "continuous_chopped_exons_seq_counting.csv", sep = "/"))
tableIntrons <- read.csv(paste(path, "continuous_chopped_introns_seq_counting.csv", sep = "/"))
expressionTable <- read.table(paste(path, "H9_ESC_RNA_seq.txt", sep = "/"), sep = "\t", header = F)

exon_stats <- creatDataESC(tableExons, expressionTable, c(-0.5,0.7), c(0.57, 0.61, 0.65), "Exons", "H9_ESC_expression_vs_asymmetry_exons.pdf", 0.68)
intron_stats <- creatDataESC(tableIntrons, expressionTable, c(-0.35,0.4), c(0.27, 0.31, 0.35), "Introns", "H9_ESC_expression_vs_asymmetry_introns.pdf",0.38)


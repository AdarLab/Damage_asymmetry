#!/usr/bin/Rscript 
###############################################################################################
# This script gets human genomic elements counting tables of dinucleodites frequency
# and generate boxplots of:
# 1) The observed and expected frequencies of all the dinucleotides on both strands of the DNA.
# 2) The observed and expected asymmetry scores.
###############################################################################################
library(ggplot2)
library(reshape2)
library(ggpubr)

path <- "your_path"
figure_path <- "figures"

TableGenes <- read.csv(paste(path, "genes_seq_counting.csv", sep = "/"), row.names = 1, header = T)

###---Random---###
pathRandom <- "random_path"
TableGenesRandom <- read.csv(paste(pathRandom, "hg38_RefSeq_uniq_coding_genes_rand_counting.csv", sep = "/"), row.names = 1, header = T)

generatePlotDataBothStrands <- function(elementTable, filtered_dinucs) {
  "This function gets the element counting table, the element name and the specific dinucleotide for plotting.
  and returns data frame of all the observed and expected frequencies of all the dinucleotide on both strands of the DNA"
  
  transcribed_dinuc <- filtered_dinucs[[1]] # Extract the vector of the dinucleotides on the transcribed strand.
  non_transcribed_dinuc <- filtered_dinucs[[2]] # "" non-transcribed strand.
  
  obs_label <- as.data.frame(rep("Observed",nrow(elementTable))) # Define column of observed labels.
  exp_label <- as.data.frame(rep("Expected",nrow(elementTable))) # "" expected labels.
  
  # Set the names of the 2 dataframes to an identical name for enabling rbind function.
  names(obs_label) <- "Type"
  names(exp_label) <- "Type"
  tableForPlot <- rbind(obs_label, exp_label)
  tableNames <- "Type"
  names(tableForPlot) <- tableNames
  
  # Scan the dinucleotides vectors and assign the frequency values of thetranscribed and non transcribed strands based on the locations of the matching labels in the "tableForPlot" table.
  for(i in 1:length(transcribed_dinuc))
  {
    firstDinucTrans <- strsplit(transcribed_dinuc[i], "")[[1]][1] # The first nucleotide of the dinucleotide on the transcribed strand.
    secondDinucTrans <- strsplit(transcribed_dinuc[i], "")[[1]][2] # The second nucleotide.
    
    firstDinucNonTrans <- strsplit(non_transcribed_dinuc[i], "")[[1]][1] # "" the non-transcribed strand.
    secondDinucNonTrans <- strsplit(non_transcribed_dinuc[i], "")[[1]][2]
    
    obs_tmp <- as.data.frame(elementTable[,transcribed_dinuc[i]]/elementTable$length*1000 + elementTable[,non_transcribed_dinuc[i]]/elementTable$length*1000)  # The current dinuleotide frequencies (counts normalized to 1Kb) on both strands.
    exp_tmp <- as.data.frame((elementTable[,firstDinucTrans]/elementTable$length)*(elementTable[,secondDinucTrans]/elementTable$length)*1000 + (elementTable[,firstDinucNonTrans]/elementTable$length)*(elementTable[,secondDinucNonTrans]/elementTable$length)*1000) # The current dinuleotide expected frequencies (counts normalized to 1Kb) based on the product of the nucleotides that make it up.
    
    # Set the names of the 2 dataframes to an identical name for enabling rbind function.    
    names(obs_tmp) <- non_transcribed_dinuc[i]
    names(exp_tmp) <- non_transcribed_dinuc[i]
    # Combine the frequency values to one column.
    dinuc_tmp <- rbind(obs_tmp, exp_tmp)
    
    # Add the new column to the table.
    tableForPlot <- cbind(tableForPlot, dinuc_tmp)
    # Update the table's names.
    tableNames <- c(tableNames, non_transcribed_dinuc[i])
    names(tableForPlot) <- tableNames
    print(head(tableForPlot))
  }
  
  return(tableForPlot)
}

generatePlotAsymmetryScore <- function(elementTable, filtered_dinucs)
{
  "This function gets the element counting table, the element name and the specific dinucleotide for plotting.
  and returns data frame of all the observed and expected strand asymmetry scores (transcribed vs non-transcribed) of all the dinucleotide"
  
  transcribed_dinuc <- filtered_dinucs[[1]]
  non_transcribed_dinuc <- filtered_dinucs[[2]]
  obs_label <- as.data.frame(rep("Observed",nrow(elementTable))) # Define column of observed labels
  exp_label <- as.data.frame(rep("Expected",nrow(elementTable))) # "" expected.
  
  # Set the names of the 2 dataframes to an identical name for enabling rbind function.
  names(obs_label) <- "Type"
  names(exp_label) <- "Type"
  tableForPlot <- rbind(obs_label, exp_label)
  tableNames <- "Type"
  names(tableForPlot) <- tableNames
  
  # Scan the dinucleotides vectors and assign the frequency values of the
  # transcribed and non transcribed strands based on the locations of the matching labels in the "tableForPlot" table.
  for(i in 1:length(transcribed_dinuc))
  {
    firstDinucTrans <- strsplit(transcribed_dinuc[i], "")[[1]][1]
    secondDinucTrans <- strsplit(transcribed_dinuc[i], "")[[1]][2]
    
    firstDinucNonTrnas <- strsplit(non_transcribed_dinuc[i], "")[[1]][1]
    secondDinucNonTrnas <- strsplit(non_transcribed_dinuc[i], "")[[1]][2]
    
    
    # Observed
    transcribed_tmp_obs <- as.data.frame((elementTable[,transcribed_dinuc[i]]/elementTable$length)*1000) # Etract the expected dinu frequencies on the transcribed strand.
    non_transcribed_tmp_obs <- as.data.frame((elementTable[,non_transcribed_dinuc[i]]/elementTable$length)*1000) # "" non transcribed.
    
    # Divide the frequencies in the sum of the dinucleotide on both strands - the asymmetry score.
    asymmetry_obs <- (transcribed_tmp_obs - non_transcribed_tmp_obs)/(transcribed_tmp_obs + non_transcribed_tmp_obs)
    
    # Expected
    transcribed_tmp_exp <- as.data.frame((elementTable[,firstDinucTrans]/elementTable$length)*(elementTable[,secondDinucTrans]/elementTable$length)*1000) # The expected dinu frequencies on the transcribed strand.
    non_transcribed_tmp_exp <- as.data.frame((elementTable[,firstDinucNonTrnas]/elementTable$length)*(elementTable[,secondDinucNonTrnas]/elementTable$length)*1000) # "" non transcribed.
    # Divide the frequencies in the sum of the dinucleotide on both strands.
    asymmetry_exp <- (transcribed_tmp_exp - non_transcribed_tmp_exp)/(transcribed_tmp_exp + non_transcribed_tmp_exp)
    
    
    # Set the names of the dataframe to the name of the currrent dinucleotide.
    names(asymmetry_obs) <- non_transcribed_dinuc[i]
    names(asymmetry_exp) <- non_transcribed_dinuc[i]
    
    # In case there are no occurrences of the given di/nucleotide in the sequence, the asymmetry score will be NaN (since the denominator is equal to 0). 
    nan_obs <- which(is.nan(asymmetry_obs[,non_transcribed_dinuc[i]]))
    nan_exp <- which(is.nan(asymmetry_exp[,non_transcribed_dinuc[i]]))

    nan_rows <- c(nan_obs,nan_exp)
    nan_rows <- unique(nan_rows) # Extract the indices of the NaN rows.
    # If there are values of NaN in the data frame set them to 0 - no occurrences.
    if(length(nan_rows) > 0)
    {
      asymmetry_obs[nan_rows, non_transcribed_dinuc[i]] <- 0
      asymmetry_exp [nan_rows, non_transcribed_dinuc[i]] <- 0
    }

    # Combine the frequency values to one column.
    dinuc_tmp <- rbind(asymmetry_obs, asymmetry_exp)

    # Add the new column to the table.
    tableForPlot <- cbind(tableForPlot, dinuc_tmp)
    # Update the table's names.
    tableNames <- c(tableNames, non_transcribed_dinuc[i])
    names(tableForPlot) <- tableNames
  }
  
  return(tableForPlot)
}

makeBoxPlotBothStrands <- function(plot.data, element, yLim, yLab, figureName) {
  "This function gets df for plotting and dinucleotide name for the label
  and creates a boxplot graph using ggplot package"
  
  customSymnum <- list(cutpoints = c(0, c(0.0001, 0.001, 0.01, 0.05)/(nrow(plot.data)/2), 1.01), symbols = c("**", "*", "***", "?", "n.s.")) # Define the p values symbols based on Bonferroni correction.
  p <- ggplot(plot.data, aes(x = variable, y = value, fill = Type)) + # Define the elements for plotting - group by "strandness".
    stat_compare_means(aes(group = Type), label = "p.signif", paired = TRUE, label.y = rep(280, 10), symnum.args = customSymnum, size = 10) +
    geom_boxplot(outlier.shape = NA, colour = "black") +
    scale_fill_manual(values = c("#A21829","#3989BD")) + theme_classic() +
    coord_cartesian(ylim = yLim)) +
    stat_summary(fun = mean, colour = "black", geom = "point", shape = 18, size = 4, show.legend = FALSE, position = position_dodge(0.75)) +
    theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title.x = element_text(color = "black", size = 28), axis.title.y = element_text(color = "black", size = 28, margin = margin(0, 20, 0, 0)),
          axis.text.x = element_text(color = "black", size = 24), axis.text.y = element_text(color = "black", size = 24), legend.position = "none") +
     labs(y = yLab, x = "")
  pdf(paste(figure_path,(paste(element, figureName, sep = "_")), sep = "/"), height = 8, width = 10)
  print(p)
  dev.off()
}

# Define the nucleotides on the transcribed and the non-transcribed strands.
filtered_dinucs_both_strands <- list(c("AA", "GA", "AG", "CA", "AC", "CC"), c("TT", "TC", "CT", "TG", "GT", "GG"))
filtered_dinucs <- list(c("AA", "GA", "AG", "CA", "AC", "CC", "AT", "TA", "CG", "GC"), c("TT", "TC", "CT", "TG", "GT", "GG", "AT", "TA", "CG", "GC"))
Tables <- list(TableGenes, TableGenesRandom) # The tables' names.
elements <- c("Genes", "Randome_genes") # The element names.

# Scan all the tables, update them and generate df for the boxplot function.
for (i in 1:length(Tables))
{
   bothStrandsTable <- generatePlotDataBothStrands(Tables[[i]], filtered_dinucs)
   both.strands.data <- melt(bothStrandsTable, id.var = "Type")
   makeBoxPlotBothStrands(both.strands.data, elements[i], c(0,300), "Dinucleotide frequencies (per Kb)", "all_dinucs_obs_vs_exp_both_strands.pdf")
  
   asymmetyTable <- generatePlotAsymmetryScore(Tables[[i]], filtered_dinucs_both_strands)
   asymmetry.data <- melt(asymmetyTable, id.var = "Type")
   makeBoxPlotBothStrands(asymmetry.data, elements[i], c(-0.4,0.2), "Asymmetry score", "all_dinucs_obs_vs_exp_asymmetry_score.pdf")
}

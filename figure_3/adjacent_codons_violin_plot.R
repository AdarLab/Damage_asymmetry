#!/usr/bin/Rscript 
# R version 3.5.2
library(ggplot2)
theme_set(theme_classic())
########################################################################################################
# This script gets table of adjacent codon counts and generate violin plots representing the 
# observed/expected frequency of all codons that start after a given nucleotide with a given nucleotide.
########################################################################################################
NUCLEOTIDES <- c("A", "T", "C", "G")
stopCodons <- c("TAA", "TAG", "TGA")

path <- "codon_usage_data"
figurePath <- paste(path, "figures/violinPlots", sep = "/")

path <- "Arabidopsis"
figurePath <- paste("codon_usage_data", "figures/violinPlots", sep = "/")

observedCodonsTable <- read.csv(paste(path, "adjacentCodonsTable.csv", sep = "/"), row.names = 1)
expectedCodonsTable <- read.csv(paste(path, "codonsTable.csv", sep = "/"), row.names = 1)

# Initialize a matrix for the Observed/Expected analysis.
normalizedCodonsTable <- as.data.frame(matrix(0L, nrow = nrow(observedCodonsTable), ncol = ncol(observedCodonsTable)))
colnames(normalizedCodonsTable) <- colnames(observedCodonsTable)
rownames(normalizedCodonsTable) <- row.names(observedCodonsTable)

# Extract the relative frequency of each codn: The number of its occurenced divided by the total number of occurences of all codons.
normalizedExpectedCodonsTable <- expectedCodonsTable / colSums(expectedCodonsTable)

# Scan the observed table and divide each value in the expected value: the relative frequency of the second codon.
for (codon_1 in row.names(observedCodonsTable)) {
  for (codon_2 in colnames(observedCodonsTable)) {
    normalizedCodonsTable[codon_1, codon_2] <- (as.numeric(observedCodonsTable[codon_1, codon_2])/sum(as.numeric(observedCodonsTable[codon_1,]))) / as.numeric(normalizedExpectedCodonsTable[codon_2, "frequency"])
    }
}

# Remove stop codons from the table rows, since they can't be first codons at all.
normalizedCodonsTable <- normalizedCodonsTable[-which(names(normalizedCodonsTable) %in% stopCodons), ]

specificNucCodonsTable <- function(normalizedCodonsTable, nuc) {
  " This function extract all the occurances of adjacent codons where the last codon of the first codon is the given nuc "
  Table_plot <- as.data.frame(matrix(0L, nrow = 0, ncol = 3)) # Initalize a table for the specific nucleotide comparison.
  # Scan the table in order to extract only pairs which the last codon of their first codon.  
  for (codon_1 in row.names(normalizedCodonsTable)) {
    for (codon_2 in colnames(normalizedCodonsTable)) {
      # Check if the last codon of the first nucleotide is indeed the given nucleotide.
      if (unlist(strsplit(codon_1, ""))[3] == nuc) {
        Group <- paste0(unlist(strsplit(codon_2, ""))[1], "XX") # Extract the first nucleotide of the second codon, and define it as the category of the codons.
        Pair <- paste(codon_1, codon_2, sep = "_") # Define the codons pair. 
        Value <- as.numeric(normalizedCodonsTable[codon_1, codon_2]) # Extract the relative frequency of the pair of codons.
        Table_plot <- rbind(Table_plot, cbind(Group, Pair, Value)) # concatenate the values to the data frame.
      }
    }
  }
  names(Table_plot) <- c("Group", "Pair", "Value")
  Table_plot$Value <- as.numeric(as.character(Table_plot$Value))
  return(Table_plot)
}

makeViolinPlot <- function(nuc) {
  " This function generate a violin plot of representing the relative frequencies of adjancent codons
    The last nucleotide of the first codon is the given nuc and each violin representing the distribution
    of the pairs where the first nucleotide of the second codon is differ "
  tablePlotNuc <- specificNucCodonsTable(normalizedCodonsTable, nuc)
  figureName <- paste0("Arabidopsis_XX", nuc, ".pdf")
  g <- ggplot(tablePlotNuc, aes(as.factor(Group), Value)) +
    geom_violin(aes(fill = Group), trim = FALSE) +
    geom_boxplot(width = 0.2, colour = "black") +
    scale_fill_manual(values = c("#a21a29", "steelblue3", "darkseagreen3", "#D76327")) +
    theme(legend.position = "none") +
    ggtitle(paste(paste0("XX", nuc), "Adjacent Codons Distribution")) +
    theme(plot.title = element_text(size = 28, hjust = 0.5), axis.text.x = element_text(size = 24, colour = "black"), axis.text.y = element_text(size = 24, colour = "black"), axis.title.x = element_text(size = 28, margin = margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y = element_text(size = 28, margin = margin(t = 0, r = 20, b = 0, l = 0)), legend.position = "none") +
    labs(x = "Second codon", y = "Relative frequencies (Obsereved/Expected)") +
    geom_hline(yintercept=1) +
    ylim(0,3.5)
  pdf(paste(figurePath, figureName, sep = "/"), height = 8, width = 10)
  print(g)
  dev.off()

}

for(nuc in NUCLEOTIDES){
  makeViolinPlot(nuc)
}

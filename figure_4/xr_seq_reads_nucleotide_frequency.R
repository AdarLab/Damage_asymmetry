#!/usr/bin/Rscript 
# R version 3.5.2
library(Biostrings)
library(ggplot2)
library(reshape2)
################################################################################
# This script gets XR-seq read sequences and generate bar plots representing
# the nucleotide frequencies along the reads in order to locate the dimer site. 
###############################################################################

globalPath <- "fastaFiles"
figuresPath <- "figures"

makeBarPlot <- function(faFilePath, figurePath, figureName) {
  
  faFile <- readDNAStringSet(faFilePath, "fasta")
  freqMat <- consensusMatrix(faFile, baseOnly=T,as.prob = T)
  tFreqMat <- t(freqMat)
  rFreqMat <- melt(tFreqMat)
  colnames(rFreqMat) <- c("Position","Base","Percentage")
  rFreqMat["Percentage"] <- rFreqMat["Percentage"]*100 # Extract the frequency of the nucleotide in percentes.
  rFreqMat <- subset(rFreqMat, Base!="other")
  freqBar <- ggplot(rFreqMat, aes(x = Position, y = Percentage ,fill=Base)) +
    geom_bar(stat = "identity", colour = "black") +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"), panel.border = element_rect(colour = "black", fill=NA, size=0.5)
          ,axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
          axis.title.x = element_text(size = 20, angle = 0, margin = margin(t = 20, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size = 20, angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
          legend.key = element_rect(colour = "transparent", fill = "white")) +
    scale_x_continuous(breaks = 1:ncol(freqMat), labels = as.character(1:ncol(freqMat))) +
    labs(y = "nucleotide frequency (%)", x = "Position") 
  pdf(paste(figurePath, figureName, sep = "/"), height = 8, width = 12)
  print(freqBar)
  dev.off()
}


for(file in Sys.glob(paste0(globalPath, "/*/*.fa")))
{
  fileName <- basename(file)
  fileName <- gsub(".fa", "", fileName)
  dirName <- basename(dirname(file))
  dir.create(paste(figuresPath, dirName, sep = "/"))
  figurePath <- paste(figuresPath, dirName, sep = "/")
  figureName <- paste0(gsub("length", "", paste0(dirName, fileName)), "_nt.pdf")
  makeBarPlot(file, figurePath, figureName)
}

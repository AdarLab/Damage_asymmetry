#!/usr/bin/Rscript 
# R version 3.5.2
library(coRdon)
library(seqinr)
library(Biostrings)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
#################################################################################
# This script gets table of codon counts and generate bar plot representing the 
# observed/expected frequency of all codons which contain given dinucleotide.
###############################################################################

nucleotides <- c("A","T","C","G") 
dinucleotides <- c(paste0("A", nucleotides),paste0("T", nucleotides), paste0("C", nucleotides), paste0("G", nucleotides)) # Define the nucleotides on the transcribed strand
getPalette = colorRampPalette(brewer.pal(4, "Paired")) # Color pallet

path <- "your_path"
figurePath <- paste(path, "figures", sep = "/")

make_bar_plot <- function(data)
{
  p <- ggplot(data, aes(x=dinucleotide, y=as.numeric(as.character(obs_vs_exp)), fill=dinucleotide)) + geom_bar(stat="identity") +
    geom_hline(yintercept=1) + ggtitle("") + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 28),
          axis.title = element_text(size = 28, colour = "black"),
          axis.text = element_text(size = 24, colour = "black"),
          axis.title.x = element_text(vjust = -1),
          axis.title.y = element_text(margin = margin(0,20,0,0)),
          legend.position="none") +
    #ylim(c(0,3)) +
    xlab("Dinucleotide contained in codons") +
    ylab("Frequency (Observed/Expected)") +
    scale_fill_manual(values = getPalette(16))
  # pdf(paste(figurePath, "codons_freq_obs_vs_exp_Arabidopsis.pdf", sep = "/"), height = 8, width = 12)
  print(p)
  # dev.off()
}

cdsFileHuman <- paste(path, "sequence", "RefSeq_hg38_CDS_13_01_20.fa", sep = "/")
cdsFileArabidopsis <- paste(path, "Arabidopsis", "sequence", "Arabidopsis_thaliana_cds_ensembl_release_47_01_06_2020.fa", sep = "/")

cdsSequence <- readSet(file = cdsFileArabidopsis) # Read the fasta file of the cds.
codons <- codonTable(cdsSequence) # 
codonsCounts <- codonCounts(codons)
codonsCounts <- as.data.frame(codonsCounts) # Convert the codonsCounts matrix into data frame.

illegal_indices <- as.numeric(row.names(codonsCounts[(codonsCounts[,"TAG"] > 1 | codonsCounts[,"TGA"] > 1 | codonsCounts[,"TAA"] > 1 | codonsCounts[,"ATG"] == 0), ]))
codonsCounts <- codonsCounts[-illegal_indices,]

codonsTable <- as.data.frame(colSums(codonsCounts)) # Calculate the total occurences of each codon.
totalCodons <- colSums(codonsTable)  # Calculate the total occurences of all codons.

observedExpectedTable <- data.frame()

for(i in 1:length(dinucleotides)){
  codonsdinuc <- subset(codonsTable, grepl(dinucleotides[i], rownames(codonsTable)) == T) # Extract only codons that contain the current dinuc.
  codonsObservedDinuc<- colSums(codonsdinuc)/totalCodons # Calculate the relative frequency of all codons that contain the current dinuc compared to all codons.
  codonsExpectedDinuc <- nrow(codonsdinuc)/nrow(codonsTable) # Calculate the number of all codons that coontain the current dinuc compared to all codons (64).
  observedVsExpectedDinuc<- codonsObservedDinuc/codonsExpectedDinuc
  observedExpectedTable <- rbind(observedExpectedTable, cbind(dinucleotides[i], observedVsExpectedDinuc))
}

names(observedExpectedTable) <- c("dinucleotide", "obs_vs_exp")

make_bar_plot(observedExpectedTable)

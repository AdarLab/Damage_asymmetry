#!/usr/bin/Rscript 
# R version 3.5.2
library(ggplot2)

##########################################################################################################################################
# This script gets multiple organisms' exons and introns counting tables of dinucleodites frequency and generate  
# scatter plots representing the avergae frequencyof the dinucleotides on the transcribed vs non-transcribed strands of exons and introns.
##########################################################################################################################################

exons <- Sys.glob("exons/tables/counting/*.csv")
introns <- Sys.glob("introns/tables/counting/*.csv")

figurePath <- "figures"


ExtractMeans <- function(exons, introns, exon_trans_vs_non, intron_trans_vs_non, dinucTrans, dinucNonTrans)
{
  for (exon in exons)
  {
    organTableExons <- read.table(exon, header = T, sep = ',')
    tmp_exons <- data.frame(basename(exon), mean((organTableExons[,dinucTrans]/organTableExons$length)*1000, na.rm = T), mean((organTableExons[,dinucNonTrans]/organTableExons$length)*1000, na.rm = T), "Exons")
    colnames(tmp_exons) <- data.frame.names
    exon_trans_vs_non <- rbind(exon_trans_vs_non,tmp_exons)
    
  }
  
  for(intron in introns)
  {
    organTableIntrons <- read.table(intron, header = T, sep = ',')
    tmp_introns <- data.frame(basename(intron), mean((organTableIntrons[,dinucTrans]/organTableIntrons$length)*1000, na.rm = T), mean((organTableIntrons[,dinucNonTrans]/organTableIntrons$length)*1000, na.rm = T), "Introns")
    colnames(tmp_introns) <- data.frame.names
    intron_trans_vs_non <- rbind(intron_trans_vs_non,tmp_introns)
  }
  
  organisms.data <- rbind(exon_trans_vs_non, intron_trans_vs_non)
  return(organisms.data)
  write.csv(organisms.data, "organisms_table.csv")
}

GenerateGraph <- function(plot.data, dinuc, axis.lim, figurePath)
{
  p <- ggplot(plot.data, aes(x = NonTranscribed, y = Transcribed , colour = Type)) + geom_point(size = 2.5, alpha = 0.8) + xlim(0,axis.lim) + ylim(0,axis.lim) +
    geom_abline(slope = 1, intercept = 0) + scale_color_manual(values = c("#D76327","steelblue3")) +
    theme(legend.title=element_blank(),
          legend.text = element_text(size=16),
          axis.text.x = element_text(color = "black", size = 14),
          axis.text.y = element_text(color = "black", size = 14),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          legend.key = element_rect(colour = "transparent", fill = "white")) +
    theme(plot.title = element_text(lineheight=.8, hjust = 0.5, size = 20),
          axis.title.y = element_text(size = 28, angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.x = element_text(size = 28, angle = 0, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
    labs(y = "Transcribed strand" , x = "Non transcribed strand") +
    ggtitle(paste("Average", dinuc, "frequency (per Kb)"))  +
    theme(legend.position = "none")
  return(p) 

  
}


data.frame.names <<- c("Organism", "Transcribed", "NonTranscribed", "Type")
exon_trans_vs_non <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(exon_trans_vs_non) <- data.frame.names
intron_trans_vs_non <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(intron_trans_vs_non) <- data.frame.names


nucleotides_trans <- c("A","T","C","G") 
nucleotides_non_trans <- c("T","A","G","C")
transcribed_dinuc <- c(paste0("A", nucleotides_trans),paste0("T", nucleotides_trans), paste0("C", nucleotides_trans), paste0("G", nucleotides_trans)) # Define the nucleotides on the transcribed strand
non_transcribed_dinuc <- c(paste0(nucleotides_non_trans,"T"),paste0(nucleotides_non_trans,"A"), paste0(nucleotides_non_trans,"G"), paste0(nucleotides_non_trans, "C")) # "" non-transcribed strand (complemantry to those on the transcribed strand respectively)

transcribed_dinuc <- c("GA", "CA", "CC", "AC", "AG")
non_transcribed_dinuc <- c("TC", "TG", "GG", "GT", "CT")

for(i in 1:length(transcribed_dinuc)){
  
  plot.data <- ExtractMeans(exons, introns, exon_trans_vs_non, intron_trans_vs_non, transcribed_dinuc[i], non_transcribed_dinuc[i])
  figureName <- paste0(non_transcribed_dinuc[i], "_exons_introns_scatter_plot.pdf")
  pdf(paste(figurePath, figureName, sep = "/"))
  print(GenerateGraph(plot.data, non_transcribed_dinuc[i], limitis))
  dev.off()
}

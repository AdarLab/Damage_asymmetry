#!/usr/bin/Rscript 
# R version 3.5.2
library(ggplot2)

##############################################################################################################################
# This script gets multiple organisms' exons and introns counting tables of di/nucleodite frequencies 
# and generate scatter plots of averge TT strand asymmetry score vs average GC content of each organism for exons and introns.
##############################################################################################################################

exons <- Sys.glob("exons/tables/counting/*.csv")
introns <- Sys.glob("introns/tables/counting/*.csv")

figurePath <- "figures"

lm_eqn <- function(df){
    m <- lm(GC_content ~ Asymmetry_score, df);    
    eq <- substitute(~~italic(r)^2~"="~r2, 
         list(r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}               

ExtractMeans <- function(exons, introns, exon_trans_vs_non, intron_trans_vs_non, dinucTrans, dinucNonTrans)
{
  for (exon in exons)
  {
    organTableExons <- read.table(exon, header = T, sep = ',')
    # Concatnate the data frame the average asymmetry score and the asymmetry GC content(%).
    tmp_exons <- data.frame(basename(exon), mean((organTableExons[,dinucTrans] - organTableExons[,dinucNonTrans])/(organTableExons[,dinucTrans] + organTableExons[,dinucNonTrans]), na.rm = T), mean(((organTableExons[,"G"] + organTableExons[,"C"])/(organTableExons[,"A"] + organTableExons[,"T"] + organTableExons[,"G"] + organTableExons[,"C"]))*100, na.rm = T), "Exons")
    colnames(tmp_exons) <- data.frame.names
    exon_trans_vs_non <- rbind(exon_trans_vs_non,tmp_exons)
    
  }
  
  for(intron in introns)
  {
    organTableIntrons <- read.table(intron, header = T, sep = ',')
    tmp_introns <- data.frame(basename(intron), mean((organTableIntrons[,dinucTrans] - organTableIntrons[,dinucNonTrans])/(organTableIntrons[,dinucTrans] + organTableIntrons[,dinucNonTrans]), na.rm = T), mean(((organTableIntrons[,"G"] +organTableIntrons[,"C"])/(organTableIntrons[,"A"] +organTableIntrons[,"T"]+organTableIntrons[,"G"] +organTableIntrons[,"C"]))*100, na.rm = T), "Introns")
    colnames(tmp_introns) <- data.frame.names
    intron_trans_vs_non <- rbind(intron_trans_vs_non,tmp_introns)
  }
  
  organisms.data <- rbind(exon_trans_vs_non, intron_trans_vs_non) 
  
  write.csv("table_exons.csv")
  write.csv("table_introns.csv")
  write.csv(organisms.data, "organisms_table.csv")
  
  return(list(exon_trans_vs_non, intron_trans_vs_non))
}

GenerateGraph <- function(plot.data, dinuc, color, element)
{
  p <- ggplot(plot.data, aes(x = Asymmetry_score, y = GC_content)) + geom_point(size = 2.5, alpha = 0.8, color = color) + 
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
    labs(x = "Average TT asymmetry score" , y = "Average GC content(%)") +
    ggtitle(paste("Plants'", element))  +
    theme(legend.position = "none") +
    annotate("text", -Inf+1, Inf-5, label = lm_eqn(plot.data), hjust = 0, vjust = 1, parse = TRUE, size = 6) 

  return(p)     
}


data.frame.names <<- c("Organism", "Asymmetry_score", "GC_content", "Type")
exon_trans_vs_non <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(exon_trans_vs_non) <- data.frame.names
intron_trans_vs_non <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(intron_trans_vs_non) <- data.frame.names


transcribed_dinuc <-"AA"
non_transcribed_dinuc <- "TT"

plot.data <- ExtractMeans(exons, introns, exon_trans_vs_non, intron_trans_vs_non, transcribed_dinuc, non_transcribed_dinuc)
exonsTable <- plot.data[[1]]
intronsTable <- plot.data[[2]]
figureNameExons <- paste0(non_transcribed_dinuc[i], "_exons_asymmetry_GC_content.pdf")
figureNameIntrons <- paste0(non_transcribed_dinuc[i], "_introns_asymmetry_GC_content.pdf")
# Exons
pdf(paste(figurePath, figureNameExons, sep = "/"))
print(GenerateGraph(exonsTable, non_transcribed_dinuc[i], "#1b7837", "exons"))
dev.off()
# Introns
pdf(paste(figurePath, figureNameIntrons, sep = "/"))
print(GenerateGraph(intronsTable, non_transcribed_dinuc[i], "#762a83", "introns"))
dev.off()


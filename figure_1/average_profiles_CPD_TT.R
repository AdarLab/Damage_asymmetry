#!/usr/bin/Rscript 
library(genomation)
library(BSgenome)
# R version 3.5.2

#######################################################################################################################################
# This script gets genomic coordinates and plot meta-gene profiles of CPDs and TT across them using genomation package of Bioconductor.
#######################################################################################################################################

########################
#--------Pathes--------#
########################

dataPah <- "your_data_path"
dinucPath <- paste(dataPah, "patterns", sep = "/")
damagePath <-  paste(dataPah, "damageIntersect", sep = "/")
elementPath <- paste(dataPah, "annotation", sep = "/")
figurePath <- paste(dataPah, "figures", sep = "/")

########################
#-----File names-------#
########################

# Element files
genes_3Kb_up <- "hg38_RefSeq_uniq_coding_genes_3Kb_up_plot.bed"
genes_3Kb_down <- "hg38_RefSeq_uniq_coding_genes_3Kb_down_plot.bed"

exons_start <- "hg38_RefSeq_uniq_chopped_exons_start_plot.bed"
exons_end <- "hg38_RefSeq_uniq_chopped_exons_end_plot.bed"

introns_start <- "hg38_RefSeq_uniq_chopped_introns_start_plot.bed"
introns_end <- "hg38_RefSeq_uniq_chopped_introns_end_plot.bed"

# Dinucs files
TT_genes_up <- "hg38_RefSeq_uniq_coding_genes_3Kb_up_plot_TT.bed"
TT_genes_down <- "hg38_RefSeq_uniq_coding_genes_3Kb_down_plot_TT.bed"
TT_exons_start <- "hg38_RefSeq_uniq_chopped_exons_start_plot_TT.bed"
TT_exons_end <- "hg38_RefSeq_uniq_chopped_exons_end_plot_TT.bed"
TT_introns_start <- "hg38_RefSeq_uniq_chopped_introns_start_plot_TT.bed"
TT_introns_end <- "hg38_RefSeq_uniq_chopped_introns_end_plot_TT.bed"

# Damage files
NCPD_genes_up <- "hg38_RefSeq_uniq_coding_genes_3Kb_up_plot_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
NCPD_genes_down <- "hg38_RefSeq_uniq_coding_genes_3Kb_down_plot_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
NCPD_exons_start <- "hg38_RefSeq_uniq_chopped_exons_start_plot_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
NCPD_exons_end <- "hg38_RefSeq_uniq_chopped_exons_end_plot_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
NCPD_introns_start <- "hg38_RefSeq_uniq_chopped_introns_start_plot_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
NCPD_introns_end <- "hg38_RefSeq_uniq_chopped_introns_end_plot_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"



make_avg_profile <- function(elementFile, transcribedFile, nonTranscribedFile, title, coordinates, labels, xLab, yLab, yLim, figureName) {

  "This function gets element coordinates
  sm1 = ScoreMatrixBin(target = transcribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 325) # Define the Damage/TT as the target to plot over the given element.
  sm2 = ScoreMatrixBin(target = nonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 325)
  
  
  plot.data = new("ScoreMatrixList",list(sm1,sm2))
 
  y_coordinates <- seq(yLim[1], yLim[2], length.out=4)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.

  pdf(paste(figurePath, figureName, sep = "/"),height=14,width=20)
  op <- par(mar = c(10,15,10,15) + 0.1)
  plotMeta(mat = plot.data, overlay = TRUE, line.col = c("goldenrod","darkgreen"), lwd = 8, xaxt = 'n', yaxt = "n",
           main = title, cex.main = 4, font.main = 1, ylim = yLim,
           ylab = "", xlab = "", cex.lab = 4, cex.axis = 2.5)
  axis(1, at = coordinates, labels = labels, cex.axis = 4))
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(ylab =  yLab, cex.lab=4, line = 8, xlab = xLab, main = "", cex = 8)
  dev.off()
}

make_avg_exons <- function(elementFile, transcribedFile, nonTranscribedFile, title, labels, xLab, yLab, yLim, figureName) {

  sm1 = ScoreMatrixBin(target = transcribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 120)
  sm2 = ScoreMatrixBin(target = nonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 120)

  plot.data = new("ScoreMatrixList",list(sm1,sm2))
  
  y_coordinates <- seq(yLim[1], yLim[2], length.out=4)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.
  
  pdf(paste(figurePath, figureName, sep = "/"),height=14,width=20)
  op <- par(mar = c(10,15,10,15) + 0.1)
  plotMeta(mat = plot.data, overlay = TRUE, line.col = c("goldenrod","darkgreen"), lwd = 8, xaxt = 'n', yaxt = "n",
           main = title, cex.main = 4, font.main = 1, ylim = yLim,
           ylab = "", xlab = "", cex.lab = 4, cex.axis = 2.5)
  axis(1, at = c(0,120), labels = labels, cex.axis = 4))
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(ylab =  yLab, cex.lab=4, line = 8, xlab = xLab, main = "", cex = 8)
  dev.off()
}

make_avg_introns <- function(elementFile, transcribedFile, nonTranscribedFile, title, labels, xLab, yLab, yLim, figureName) {

  sm1 = ScoreMatrixBin(target = transcribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 180)
  sm2 = ScoreMatrixBin(target = nonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 180)

  plot.data = new("ScoreMatrixList",list(sm1,sm2))
  
  y_coordinates <- seq(yLim[1], yLim[2], length.out=4)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.

  pdf(paste(figurePath, figureName, sep = "/"),height=14,width=20)
  op <- par(mar = c(10,15,10,15) + 0.1)
  plotMeta(mat = plot.data, overlay = TRUE, line.col = c("goldenrod","darkgreen"), lwd = 8, xaxt = 'n', yaxt = "n",
           main = title, cex.main = 4, font.main = 1, ylim = yLim,
           ylab = "", xlab = "", cex.lab = 4, cex.axis = 2.5)
  axis(1, at = c(0,180), labels = labels, cex.axis = 4))
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(ylab =  yLab, cex.lab=4, line = 8, xlab = xLab, main = "", cex = 8)
  dev.off()
}

##############
# ---Genes---#
##############

 elementFileStart <- readGeneric(paste(elementPath, genes_3Kb_up, sep = "/"), strand = 6) # Use readGeneric function to read the bed file containing the coordinates of the 3Kb to 10 Kb downstream of the TSS.
 elementFileEnd <- readGeneric(paste(elementPath, genes_3Kb_down, sep = "/"), strand = 6) # 10 Kb upstream and 3 Kb downstream of the TES.
 
 #  CPD
 transcribedFileStart <- readGeneric(paste0(damagePath, "/transcribed/NCPD/", NCPD_genes_up), strand = 6) # CPD reads overlapping the 3 Kb upstream to 10 Kb downstream of the TSS on the transcribed strand.
 nonTranscribedFileStart <- readGeneric(paste0(damagePath, "/nonTranscribed/NCPD/", NCPD_genes_up), strand = 6) # "" the non-transcribed strand.
 
 transcribedFileEnd <- readGeneric(paste0(damagePath, "/transcribed/NCPD/", NCPD_genes_down), strand = 6)
 nonTranscribedFileEnd <- readGeneric(paste0(damagePath, "/nonTranscribed/NCPD/", NCPD_genes_down), strand = 6)
 
 
 make_avg_profile(elementFileStart, transcribedFileStart, nonTranscribedFileStart, "CPD frequency", c(0,75,325), c("-3Kb", "TSS", "10Kb"), "Relative distance from TSS", "Average read count", c(0.025, 0.08), "genes_up_CPD.pdf")
 make_avg_profile(elementFileEnd, transcribedFileEnd, nonTranscribedFileEnd, "CPD frequency", c(0,250,325), c("-10Kb", "TES", "3Kb"), "Relative distance from TES", "Average read count", c(0.025, 0.08), "genes_down_CPD.pdf")
 
# TT
 
 elementFileStart <- readGeneric(paste(elementPath, genes_3Kb_up, sep = "/"), strand = 6)
 transcribedFileStart <- readGeneric(paste0(dinucPath, "/transcribedStrand/", TT_genes_up), strand = 6)
 nonTranscribedFileStart <- readGeneric(paste0(dinucPath, "/nonTranscribedStrand/", TT_genes_up), strand = 6)
 
 elementFileEnd <- readGeneric(paste(elementPath, genes_3Kb_down, sep = "/"), strand = 6)
 transcribedFileEnd <- readGeneric(paste0(dinucPath, "/transcribedStrand/", TT_genes_down), strand = 6)
 nonTranscribedFileEnd <- readGeneric(paste0(dinucPath, "/nonTranscribedStrand/", TT_genes_down), strand = 6)
 
 
 make_avg_profile(elementFileStart, transcribedFileStart, nonTranscribedFileStart, "TT frequency", c(0,75,325), c("-3Kb", "TSS", "10Kb"), "Relative distance from TSS", "Average frequency", c(0.1,0.55), "genes_up_TT.pdf")
 make_avg_profile(elementFileEnd, transcribedFileEnd, nonTranscribedFileEnd, "TT frequency", c(0,250,325), c("-10Kb", "TES", "3Kb"), "Relative distance from TES", "Average frequency", c(0.1,0.55), "genes_down_TT.pdf")

#############
#---Exons---#
#############
elementFileStart <- readGeneric(paste(elementPath, exons_start, sep = "/"), strand = 6)
elementFileEnd <- readGeneric(paste(elementPath, exons_end, sep = "/"), strand = 6)

# CPD
transcribedFileStart <- readGeneric(paste0(damagePath, "/transcribed/NCPD/", NCPD_exons_start), strand = 6)
nonTranscribedFileStart <- readGeneric(paste0(damagePath, "/nonTranscribed/NCPD/", NCPD_exons_start), strand = 6)

transcribedFileEnd <- readGeneric(paste0(damagePath, "/transcribed/NCPD/", NCPD_exons_end), strand = 6)
nonTranscribedFileEnd <- readGeneric(paste0(damagePath, "/nonTranscribed/NCPD/", NCPD_exons_end), strand = 6)


make_avg_exons(elementFileStart, transcribedFileStart, nonTranscribedFileStart, "CPD frequency", c("Exon start", "120b"), "Relative distance from exon start", "Average read count", c(0.03, 0.065), "exons_start_CPD.pdf")
make_avg_exons(elementFileEnd, transcribedFileEnd, nonTranscribedFileEnd, "CPD frequency", c("-120b", "Exon end"), "Relative distance from exon end", "Average read count", c(0.03, 0.065), "exons_end_CPD.pdf")

# TT
transcribedFileStart <- readGeneric(paste0(dinucPath, "/transcribedStrand/", TT_exons_start), strand = 6)
nonTranscribedFileStart <- readGeneric(paste0(dinucPath, "/nonTranscribedStrand/", TT_exons_start), strand = 6)

transcribedFileEnd <- readGeneric(paste0(dinucPath, "/transcribedStrand/", TT_exons_end), strand = 6)
nonTranscribedFileEnd <- readGeneric(paste0(dinucPath, "/nonTranscribedStrand/", TT_exons_end), strand = 6)


make_avg_exons(elementFileStart, transcribedFileStart, nonTranscribedFileStart, "TT frequency", c("Exon start", "120b"), "Relative distance from exon start", "Average frequency", c(0.1, 0.36), "exons_start_TT.pdf")
make_avg_exons(elementFileEnd, transcribedFileEnd, nonTranscribedFileEnd, "TT frequency", c("-120b", "Exon end"), "Relative distance from exon end", "Average frequency", c(0.1, 0.36), "exons_end_TT.pdf")

#############
#--Introns--#
#############

elementFileStart <- readGeneric(paste(elementPath, introns_start, sep = "/"), strand = 6)
elementFileEnd <- readGeneric(paste(elementPath, introns_end, sep = "/"), strand = 6)

# CPD
transcribedFileStart <- readGeneric(paste0(damagePath, "/transcribed/NCPD/", NCPD_introns_start), strand = 6)
nonTranscribedFileStart <- readGeneric(paste0(damagePath, "/nonTranscribed/NCPD/", NCPD_introns_start), strand = 6)
transcribedFileEnd <- readGeneric(paste0(damagePath, "/transcribed/NCPD/", NCPD_introns_end), strand = 6)
nonTranscribedFileEnd <- readGeneric(paste0(damagePath, "/nonTranscribed/NCPD/", NCPD_introns_end), strand = 6)


make_avg_introns(elementFileStart, transcribedFileStart, nonTranscribedFileStart, "CPD frequency", c("Intron start", "1.874Kb"), "Relative distance from intron start", "Average read count", c(0.03, 0.065), "introns_start_CPD.pdf")
make_avg_introns(elementFileEnd, transcribedFileEnd, nonTranscribedFileEnd, "CPD frequency", c("-1.874Kb", "Intron end"), "Relative distance from intron end", "Average read count",c(0.03, 0.065), "introns_end_CPD.pdf")

# TT
transcribedFileStart <- readGeneric(paste0(dinucPath, "/transcribedStrand/", TT_introns_start), strand = 6)
nonTranscribedFileStart <- readGeneric(paste0(dinucPath, "/nonTranscribedStrand/", TT_introns_start), strand = 6)

transcribedFileEnd <- readGeneric(paste0(dinucPath, "/transcribedStrand/", TT_introns_end), strand = 6)
nonTranscribedFileEnd <- readGeneric(paste0(dinucPath, "/nonTranscribedStrand/", TT_introns_end), strand = 6)


make_avg_introns(elementFileStart, transcribedFileStart, nonTranscribedFileStart, "TT frequency", c("Intron start", "1.874Kb"), "Relative distance from intron start", "Average frequency", c(0.1, 0.36), "introns_start_TT.pdf")
make_avg_introns(elementFileEnd, transcribedFileEnd, nonTranscribedFileEnd, "TT frequency", c("-1.874Kb", "Intron end"), "Relative distance from intron end", "Average frequency", c(0.1, 0.36), "introns_end_TT.pdf")

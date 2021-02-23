#!/usr/bin/Rscript 
library(genomation)
library(BSgenome)
# R version 3.5.2

##########################################################################################################################################
# This script gets enhancer's coordinates and plot meta-gene profiles of CPDs and TT across them using genomation package of Bioconductor.
##########################################################################################################################################

path <- "your_path"
transcribedPath <- paste(path, "damageIntersect/transcribed", sep = "/")
nonTranscribedPath <- paste(path, "damageIntersect/nonTranscribed", sep = "/")
figurePath <- paste0(path, "/figures")

make_avg_profile <- function(elementFile, transcribedFile, nonTranscribedFile, yLim, title, figureName, colors, yLab, coordinates, xLabels) {
  
  elementFile <- readGeneric(elementFile, strand = 6) #
  transcribedFile <- readGeneric(transcribedFile, strand = 6)
  nonTranscribedFile <- readGeneric(nonTranscribedFile, strand = 6)
  
  sm1 = ScoreMatrixBin(target = transcribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 75)
  sm2 = ScoreMatrixBin(target = nonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 75)

  plot.data = new("ScoreMatrixList",list(sm1,sm2))
  
  y_coordinates <- seq(yLim[1], yLim[2], length.out=4)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.
  
  pdf(paste(figurePath, figureName, sep = "/"),height=14,width=20)
  op <- par(mar = c(10,15,10,15) + 0.1)
  plotMeta(mat = plot.data, overlay = TRUE, line.col = colors, lwd = 10, xaxt = 'n', yaxt = 'n', ylim = yLim,
           main = paste(title, "frequency"), cex.main = 4, font.main = 1, ylab = "", xlab = "", cex.lab = 4, cex.axis = 2.5))
  axis(1, at = coordinates, labels = xLabels, cex.axis = 4))
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(ylab =  yLab, cex.lab=4, line = 8, xlab = "Relative distance from eTSS", main = "", cex = 8)
  dev.off()
}

# Enhancer's coordinates.
enhancers_right <- "hg38_enhancers_750_right_sorted.bed" # 750bp upstream the eTSS.
enhancers_left <- "hg38_enhancers_750_left_sorted.bed" # 750bp downstream the eTSS.

# CPD coordinates obtained by Damage-seq across enhancers.
CPD_plus_strand_right <- "hg38_enhancers_750_right_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
CPD_minus_strand_right <- "hg38_enhancers_750_right_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"

CPD_plus_strand_left <- "hg38_enhancers_750_left_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
CPD_minus_strand_left <- "hg38_enhancers_750_left_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"

# TT coordinates across enhancers.
TT_plus_strand_right <- "hg38_enhancers_750_right_TT.bed"
TT_minus_strand_right <- "hg38_enhancers_750_right_TT.bed"

TT_plus_strand_left <- "hg38_enhancers_750_left_TT.bed"
TT_minus_strand_left <- "hg38_enhancers_750_left_TT.bed"

# Avergae profile CPD.
make_avg_profile(enhancers_right, CPD_plus_strand_right, CPD_minus_strand_right, c(0.046, 0.072), "CPD", "hg38_enhancers_plot_CPD_strands_750_right.pdf", c("goldenrod", "darkgreen"), "Average read count", c(0,50,75), c("eTSS", "500b", "750b"))
make_avg_profile(enhancers_left, CPD_plus_strand_left, CPD_minus_strand_left, c(0.046, 0.072), "CPD", "hg38_enhancers_plot_CPD_strands_750_left.pdf", c("goldenrod","darkgreen"), "Average read count", c(0,25,75), c("-750b", "-500b", "eTSS"))

# Avergae profile TT.
make_avg_profile(enhancers_right, TT_plus_strand_right, TT_minus_strand_right, c(0.18, 0.31), "TT", "hg38_enhancers_plot_TT_strands_750_right.pdf", c("goldenrod","darkgreen"), "TT frequency", c(0,50,75), c("eTSS", "500b", "750b"))
make_avg_profile(enhancers_left, TT_plus_strand_left, TT_minus_strand_left, c(0.18, 0.31), "TT", "hg38_enhancers_plot_TT_strands_750_left.pdf", c("goldenrod","darkgreen"), "TT frequency", c(0,25,75), c("-750b", "-500b", "eTSS"))




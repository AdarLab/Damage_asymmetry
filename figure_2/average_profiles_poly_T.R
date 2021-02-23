#!/usr/bin/Rscript 
library(genomation)
library(BSgenome)
library(Biostrings)
library(formattable)
# R version 3.5.2

############################################################################################################################
# This script gets genes coordinates and plot meta-gene profiles of TT across them using genomation package of Bioconductor.
############################################################################################################################
path <- "patterns"
figurePath <<- "figures"
Ts <- c("T", "TT", "TTT", "TTTT", "TTTTT", "TTTTTT", "TTTTTTT", "TTTTTTTT", "TTTTTTTTT", "TTTTTTTTTT")

elementFilesPath <- "annotation"

make_profiles <- function(path, element, file, elementFile, Ts, bins, coordinates, labels, x_label) {
  fileName <- basename(file)
  y_label <- paste("T", nchar(Ts), sep = "_")

  transcribedFile <- readGeneric(paste0(path, "/transcribedStrand/", fileName), strand = 6)
  nonTranscribedFile <- readGeneric(paste0(path, "/nonTranscribedStrand/", fileName), strand = 6)
  elementFile <- readGeneric(elementFile, strand = 6)

  sm1 <- ScoreMatrixBin(target = transcribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  sm2 <- ScoreMatrixBin(target = nonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)

  plot.data <- new("ScoreMatrixList", list(sm1, sm2))
  
  sm1_min <- min(colMeans(sm1))
  sm2_min <- min(colMeans(sm2))
  
  y_min <- min(sm1_min, sm2_min) # Lower Y limit
  
  sm1_max <- max(colMeans(sm1))
  sm2_max <- max(colMeans(sm2))

  y_max <- max(sm1_max, sm2_max) # Upper Y limit.
  
  y_coordinates <- seq(y_min, y_max, length.out=3)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.
  

  pdf(paste0(figurePath, "/", element, "_", y_label, ".pdf"),height=14,width=20)
  op <- par(mar = c(7, 10, 7, 1) + 0.1)
  plotMeta(
    mat = plot.data, overlay = TRUE, line.col = c("goldenrod", "darkgreen"), lwd = 8, xaxt = "n", yaxt = 'n', ylim = c(y_min, y_max),
    ylab = "", xlab = "", cex.lab = 4, cex.axis = 2.5
  )
  axis(1, at = coordinates, labels = labels, cex.axis = 4)
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(ylab = y_label, cex.lab = 4, line = 6, xlab = x_label)
  dev.off()
}

# Scan all the poly T sequences and send them to the average profile function.
for (i in 1:length(Ts))
{
   # Genes up
   genes_file <- paste(path, paste0("hg38_RefSeq_uniq_coding_genes_3Kb_up_plot_", Ts[i], ".bed"), sep = "/")
   make_profiles(path, "genes_up", genes_file, paste(elementFilesPath, "hg38_RefSeq_uniq_coding_genes_3Kb_up_plot.bed", sep = "/"), Ts[i], 13000/(40*(1+0.3*i)), c(0, 3000/(40*(1+0.3*i)), 13000/(40*(1+0.3*i))), c("-3Kb", "TSS", "10Kb"), "Relative distance from TSS")
}


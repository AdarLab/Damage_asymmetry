#!/usr/bin/Rscript 
library(genomation)
library(BSgenome)
library(Biostrings)
###################################################################################################
# This script generates average profiles of ATAC-Seq data across human genomic elements.
###################################################################################################

########################
#-----Num of reads-----#
########################

numOfreads_NO_UV <- 81763487 
numOfreads_2h <- 70850480 
########################
#--------Pathes--------#
########################

ATAC_seq_path <- "ATAC_seq_data"
elementPath <- "annotation"
figurePath <- "figures"

########################
#-----File names-------#
########################

# Element files
genes_3Kb_up <- "genes_start_plot.bed"
genes_3Kb_down <- "genes_end_plot.bed"

exons_start <- "chopped_exons_no_first_last_start_plot.bed"
exons_end <- "chopped_exons_no_first_last_end_plot.bed"

introns_start <- "chopped_introns_no_first_last_start_plot.bed"
introns_end <- "chopped_introns_no_first_last_end_plot.bed"


# Repair file Names have the same name as the element files.

make_profile <- function(element, exp, fileName, elementFile, bins, xLab, coordinates, labels, counts, title, figureName, color, yLim)
{
  # This function generates average profiles of damage/repair on the transcribed and non-transcribed strand 3 Kb upstream
  # and 10 Kb downstream of the TSS of human genes.
  
  ATAC_seq_file <- readGeneric(paste(ATAC_seq_path, exp, fileName, sep = "/"))
  elementFile <- readGeneric(paste(elementPath, elementFile, sep = "/"), strand = 6)
  
  sm = ScoreMatrixBin(target = ATAC_seq_file, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  y_min <- min(colMeans(sm))
  y_max <- max(colMeans(sm))

  y_coordinates <- seq(y_min, y_max, length.out=4)
  
  y_coordinates <- seq(yLim[1], yLim[2], length.out=4)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.
  print(y_coordinates)
  
  y_coordinates <- seq(yLim[1], yLim[2], length.out=4)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.
  
  pdf(paste0(figurePath, "/", element, "/", figureName, ".pdf"),height=14,width=20)
  op <- par(mar = c(10,15,10,15) + 0.1)
  plotMeta(mat = sm, line.col = color, lwd = 8, xaxt = 'n', yaxt = 'n', ylim = yLim,
           ylab = "", xlab = "", cex.lab = 4, cex.axis = 2.5, font.main = 1)
  axis(1,at = coordinates, labels = labels, cex.axis = 4)
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(main = title, cex.main=4, font.main = 1, ylab = "Average ATAC-seq read count\n", cex.lab=4, line = 6, xlab = paste("Relative distance from", xLab))
  dev.off()  
}

make_profile_exons_introns <- function(exonsFile, intronsFile, exonaATAC, intronsATAC, bins, xLab, coordinates, labels, counts, title, figureName, yLim)
{
  # This function generates average profiles of ATAC-seq data across exons and introns.
  exons_ATAC_seq <- readGeneric(exonsATAC)
  introns_ATAC_seq <- readGeneric(intronsATAC)
  exonsFile <- readGeneric(paste(elementPath, exonsFile, sep = "/"), strand = 6)
  intronsFile <- readGeneric(paste(elementPath, intronsFile, sep = "/"), strand = 6)

  sm_exons = ScoreMatrixBin(target = exons_ATAC_seq, windows = exonsFile, strand.aware = TRUE, bin.num = bins)
  sm_introns = ScoreMatrixBin(target = introns_ATAC_seq, windows = intronsFile, strand.aware = TRUE, bin.num = bins)


  y_coordinates <- seq(yLim[1], yLim[2], length.out=4)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.
  print(y_coordinates)
  
  plot.data <- new("ScoreMatrixList", list(sm_exons, sm_introns))

  pdf(paste0(figurePath, "/", element, "/", figureName, ".pdf"),height=14,width=20)
  op <- par(mar = c(10,15,10,15) + 0.1)
  plotMeta(mat = plot.data, line.col = c(pink, gray), lwd = 10, xaxt = 'n', yaxt = 'n', ylim = yLim,
           ylab = "", xlab = "", cex.lab = 4, cex.axis = 2.5, font.main = 1)
  axis(1,at = coordinates, labels = labels, cex.axis = 4)
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(main = title, cex.main=4, font.main = 1, ylab = "Average ATAC-seq read count\n", cex.lab=4, line = 6, xlab = "Relative distance from element start")
  dev.off()
}

exps <- c("NO_UV", "2h")
elementFiles <- c(genes_3Kb_up, genes_3Kb_down, exons_start, exons_end, introns_start, introns_end)
elements <- c("genes", "genes", "exons", "exons", "introns", "introns")
bins <- c(325, 325, 120, 120, 180, 180)
coordinates <- list(c(0,75,325),c(0,250,325), c(0,120), c(0,120), c(0,180), c(0,180))
labels <- c("TSS", "TES", "Exon start", "Exon end", "Intron start", "Intron end")
x_labels <- list(c("-3Kb", "TSS", "10Kb"), c("-10Kb", "TES", "3Kb"), c("Exon start", "120b"), c("-120b", "Exon End"),  c("Intron start", "1.874Kb"), c("-1.874Kb", "Intron End"))
numOfreads <- c(numOfreads_NO_UV, numOfreads_2h)
y_limitis <- list(c(5,15),c(5,15), c(0.51,0.7),  c(6,9),  c(0.51,0.7),  c(6,9))
colors <- c("dodgerblue4", "red2")

for(i in 1:length(exps)){
  for(j in 1:length(elementFiles))
  {
    readFile <- Sys.glob(paste(ATAC_seq_path, exps[i], paste0("*", elementFiles[j]), sep = "/"))
    fileName <- basename(readFile)
    print(elementFiles[j])
    print(fileName)
    make_profile(elements[j], exps[i], fileName, elementFiles[j], bins[j], labels[j], coordinates[[j]], x_labels[[j]], numOfreads[i], "ATAC-seq reads", paste(exps[i], gsub(".bed", "", elementFiles[j]), "ATAC_seq_profile", sep = "_"), colors[i], y_limitis[[j]])
    
    } 
}  

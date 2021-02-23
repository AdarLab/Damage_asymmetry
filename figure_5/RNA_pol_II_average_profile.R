#!/usr/bin/Rscript 
# R version 3.5.2
library(genomation)
library(BSgenome)
library(Biostrings)

###################################################################################################
# This script generates average profiles of RNA pol II ChIP-Seq data across human genomic elements.
###################################################################################################

########################
#-----Num of reads-----#
########################
# repair
numOfreads_NHF1_CPD <- 11100235 # NHF1 CPD.
numOfreads_CSB_CPD <- 19152949 # CSB CPD
numOfreads_XPC_CPD <- 14856828 # XPC CPD

# pol II
numOfreads_NO_UV <- 29442869 
numOfreads_1h <- 14186642 
numOfreads_2h <- 14354300 
########################
#--------Pathes--------#
########################

repairPath <- "xr_seq_data"
pol_II_path <- "ChIP_seq_data"
elementPath <- "/annotation/averageProfilesData"
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

make_profile <- function(element, cell, type, fileName, elementFile, bins, xLab, coordinates, labels, counts, title, figureName, colors, yLim)
{
  # This function generates average profiles of damage/repair on the transcribed and non-transcribed strand 3 Kb upstream
   
  transcribedFile <- readGeneric(paste(repairPath, "transcribedStrand", cell, fileName, sep = "/"), strand = 6)
  nonTranscribedFile <- readGeneric(paste(repairPath,"nonTranscribedStrand", cell, fileName, sep = "/"), strand = 6)
  elementFile <- readGeneric(paste(elementPath, elementFile, sep = "/"), strand = 6)

  sm1 = ScoreMatrixBin(target = transcribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  sm2 = ScoreMatrixBin(target = nonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  
  #Normalize the data
  sm1 <- (sm1/counts)*10^9
  sm2 <- (sm2/counts)*10^9
  
  plot.data = new("ScoreMatrixList",list(sm1,sm2))
  
  y_coordinates <- seq(yLim[1], yLim[2], length.out=4)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.
  
  pdf(paste0(figurePath, "/repair/", element, "/", figureName, ".pdf"),height=12,width=20)
  op <- par(mar = c(10,15,10,15) + 0.1)
  plotMeta(mat = plot.data, overlay = TRUE, line.col = colors, lwd = 8, xaxt = 'n', yaxt = 'n',ylim = yLim,
           ylab = "", xlab = "", cex.lab = 4, cex.axis = 2.5, font.main = 1)
  axis(1,at = coordinates, labels = labels, cex.axis = 4, lwd = 0, line = 1.5)
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(main = title, cex.main=4, font.main = 1, ylab = "Average counts per bil reads\n", cex.lab=4, line = 6, xlab = paste("Relative distance from", xLab))
  dev.off()
}

make_profile_both_strands <- function(element, cell, type, fileName, elementFile, bins, xLab, coordinates, labels, counts, title, figureName, color, yLim)
{
  # This function generates average profiles of RNA pol II data across exons and introns.
  
  transcribedFile <- readGeneric(paste(repairPath, "transcribedStrand", cell, fileName, sep = "/"), strand = 6)
  nonTranscribedFile <- readGeneric(paste(repairPath,"nonTranscribedStrand", cell, fileName, sep = "/"), strand = 6)
  elementFile <- readGeneric(paste(elementPath, elementFile, sep = "/"), strand = 6)
  
  sm1 = ScoreMatrixBin(target = transcribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  sm2 = ScoreMatrixBin(target = nonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)  
  sm_both_strands <- (colMeans(sm1)+ colMeans(sm2))/2
  #Normalize the data
  sm_both_strands <- (sm_both_strands/counts)*10^9
    
 
  y_coordinates <- seq(yLim[1], yLim[2], length.out=4)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.
  
  pdf(paste0(figurePath, "/repair_both_strands/", element, "/", figureName, ".pdf"),height=12,width=20)
  op <- par(mar = c(10,15,10,15) + 0.1)
  
    matplot(sm_both_strands,type="l",col = color,lty = c(1,1), lwd = 8,
          xaxt = 'n', ylab = "", ylim = yLim,
          cex.lab = 4, cex.axis = 2.5, main = "", font.main = 1, cex.main = 5)
  axis(1,at = coordinates, labels = labels, cex.axis = 2.5)
  title(main = title, cex.main=4, font.main = 1, ylab = "Average counts per bil reads\n", cex.lab=4, line = 6, xlab = paste("Relative distance from", xLab))

  dev.off()
}

make_profile_pol <- function(element, exp, fileName, elementFile, bins, xLab, coordinates, labels, counts, title, figureName, color, yLim)
{
  # This function generates average profiles of damage/repair on the transcribed and non-transcribed strand 3 Kb upstream
  # and 10 Kb downstream of the TSS of human genes.
  
  chipFile <- readGeneric(paste(pol_II_path, exp, fileName, sep = "/"))
  elementFile <- readGeneric(paste(elementPath, elementFile, sep = "/"), strand = 6)
  
  sm = ScoreMatrixBin(target = chipFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  
  #Normalize the data
  #sm <- (sm/counts)
  y_min <- min(colMeans(sm))
  y_max <- max(colMeans(sm))

  y_coordinates <- seq(y_min, y_max, length.out=4)
  
  y_coordinates <- seq(yLim[1], yLim[2], length.out=4)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.
  #print(y_coordinates)
  
  pdf(paste0(figurePath, "/pol_II/", element, "/", figureName, ".pdf"),height=14,width=20)
  op <- par(mar = c(10,15,10,15) + 0.1)
  plotMeta(mat = sm, line.col = color, lwd = 8, xaxt = 'n', yaxt = 'n', ylim = yLim,
           ylab = "", xlab = "", cex.lab = 4, cex.axis = 2.5, font.main = 1)
  axis(1,at = coordinates, labels = labels, cex.axis = 4)
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(main = title, cex.main=4, font.main = 1, ylab = "Average ChIP-seq read count\n", cex.lab=4, line = 6, xlab = paste("Relative distance from", xLab))
  dev.off()  
}

cellTypes <- c("NHF1", "CSB", "XPC")
numOfreads <- c(numOfreads_NHF1_CPD, numOfreads_CSB_CPD, numOfreads_XPC_CPD)
elementFiles <- c(genes_3Kb_up, genes_3Kb_down, exons_start, exons_end, introns_start, introns_end)
elements <- c("genes", "genes", "exons", "exons", "introns", "introns")
bins <- c(325, 325, 120, 120, 180, 180)
coordinates <- list(c(0,75,325),c(0,250,325), c(0,120), c(0,120), c(0,180), c(0,180))
labels <- c("TSS", "TES", "Exon start", "Exon end", "Intron start", "Intron end")
x_labels <- list(c("-3Kb", "TSS", "10Kb"), c("-10Kb", "TES", "3Kb"), c("Exon start", "120b"), c("-120b", "Exon End"),  c("Intron start", "1.874Kb"), c("-1.874Kb", "Intron End"))
y_limitis <- list(c(0,50), c(0,50), c(1.05e-8,1.25e-8), c(0,36), c(1.05e-8,1.25e-8), c(0,36))


exps <- c("NO_UV", "1h", "2h")
numOfreadsChIP <- c(numOfreads_NO_UV, numOfreads_1h, numOfreads_2h)
y_limitis_pol_II <- list(c(7,34), c(7,34), c(0.31,0.37),  c(9,15),  c(0.31,0.37),  c(9,15))
colors <- c("dodgerblue", "forestgreen", "darkorchid")

for(i in 1:length(exps)){
  for(j in 1:length(elementFiles))
  {
    chipFile <- Sys.glob(paste(pol_II_path, exps[i], paste0("*", elementFiles[j]), sep = "/"))
    fileName <- basename(chipFile)
    print(elementFiles[j])
    print(fileName)
    make_profile_pol(elements[j], exps[i], fileName, elementFiles[j], bins[j], labels[j], coordinates[[j]], x_labels[[j]], numOfreadsChIP[i], "RNA pol II ChIP-seq reads", paste(exps[i], gsub(".bed", "", elementFiles[j]), "RNA_pol_II_profile", sep = "_"), colors[i], y_limitis_pol_II[[j]])
    
    } 
} 

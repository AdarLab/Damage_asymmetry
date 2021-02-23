#!/usr/bin/Rscript 
# R version 3.5.2
library(genomation)

##############################################################################
# This script gets genomic coordinates of enhancers and plot:
# 1) Average profiles of XR-seq reads.
# 2) Average profiles of repair normalized to damage (CPD Damage-seq reads).
# 3) Average profiles of repair normalized to TT counts.
##############################################################################

########################
#-----Num of reads-----#
########################

# Dinucs
numOfreadsTT <- 576715193 # TT in the genome

# Damage
numOfreadsNCPD <- 68855889 # CPD dypirimidines

# Repair
numOfreads_NHF1_CPD <- 11100235 # NHF1 CPD.
numOfreads_CSB_CPD <- 19152949 # CSB CPD
numOfreads_XPC_CPD <- 14856828 # XPC CPD

########################
#--------Pathes--------#
########################

path <- "/vol/sci/bio/data/sheera.adar/elisheva.h/enhancers/eTSS_analysis_data"
elementFile <- paste(path, "stranded_enhancers/given_midpoint/hg38_enhancers_plot_plus_1500_sorted.bed", sep = "/")
enhancers_right <- paste(path, "stranded_enhancers/hg38_enhancers_750_right_sorted.bed", sep = "/")
enhancers_left <- paste(path, "stranded_enhancers/hg38_enhancers_750_left_sorted.bed", sep = "/")

# transcribedPath <- paste(path, "damageIntersect/given_midpoint/transcribed/hg38", sep = "/")
# nonTranscribedPath <- paste(path, "damageIntersect/given_midpoint/nonTranscribed/hg38", sep = "/")
repairPlusPath <- paste(path, "intersectedData/repair/plus/hg38", sep = "/")
repairMinusPath <- paste(path, "intersectedData/repair/minus/hg38", sep = "/")
figurePath <- "/vol/sci/bio/data/sheera.adar/elisheva.h/enhancers/figure0s/averageProfiles/eTSS"


# Repair file Names have the same name as the element files.

make_avg_profile <- function(elementFile, transcribedFile, nonTranscribedFile, numOfreads, title, figureName, yLim, x_coordinates, x_labels) {
  
  elementFile <- readGeneric(elementFile, strand = 6)
  transcribedFile <- readGeneric(transcribedFile, strand = 6)
  nonTranscribedFile <- readGeneric(nonTranscribedFile, strand = 6)
  
  sm1 = ScoreMatrixBin(target = transcribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 75)
  sm2 = ScoreMatrixBin(target = nonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 75)
  
  sm1 <- sm1/numOfreads
  sm2 <- sm2/numOfreads
  
  plot.data = new("ScoreMatrixList",list(sm1,sm2))
  
  y_coordinates <- seq(yLim[1], yLim[2], length.out=4)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.
  
  pdf(paste(figurePath, figureName, sep = "/"),height=14,width=20)
  op <- par(mar = c(7,10,7,1) + 0.1)
  plotMeta(mat = plot.data, overlay = TRUE, line.col = c("palevioletred","seagreen"), lwd = 8, xaxt = 'n',
           main = title, cex.main = 4, font.main = 1, yaxt = 'n', ylim = yLim,
           ylab = "", xlab = "", cex.lab = 4, cex.axis = 2.5)
  axis(1, at = x_coordinates, labels = x_labels, cex.axis = 4, lwd = 0, line = 2 )
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(ylab =  "Average normalized XR-seq read counts", cex.lab=4, line = 6, xlab = "Relative distance from eTSS", main = "", cex = 6)
  dev.off()
}


make_profile_ratio <- function(elementFile, damageTranscribedFile, damageNonTranscribedFile, repairTranscribedFile, repairNonTranscribedFile, repairNumOfreads, title, figureName, yLim, x_coordinates, x_labels) {
  # This function generates average ratio profiles of damage & repair on the transcribed and non-transcribed strands.

  elementFile <- readGeneric(elementFile, strand = 6)
  damageTranscribedFile <- readGeneric(damageTranscribedFile, strand = 6)
  damageNonTranscribedFile <- readGeneric(damageNonTranscribedFile, strand = 6)
  repairTranscribedFile <- readGeneric(repairTranscribedFile, strand = 6)
  repairNonTranscribedFile <- readGeneric(repairNonTranscribedFile, strand = 6)

  sm1 <- ScoreMatrixBin(target = damageTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 75)
  sm2 <- ScoreMatrixBin(target = damageNonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 75)
  sm3 <- ScoreMatrixBin(target = repairTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 75)
  sm4 <- ScoreMatrixBin(target = repairNonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 75)

  # Normalize the data
  sm3 <- sm3/repairNumOfreads
  sm4 <- sm4/repairNumOfreads

  # Divide the mean value of damage in the mean value of repair.
  ratio_transcribed <- colMeans(sm3) / colMeans(sm1)
  ratio_non_transcribed <- colMeans(sm4) / colMeans(sm2)
  
  y_coordinates <- seq(yLim[1], yLim[2], length.out=4)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.

  pdf(paste(figurePath, paste0("damage_repair_", figureName), sep = "/"), height = 14, width = 20)
  op <- par(mar = c(10, 15, 10, 10) + 0.1)
  matplot(cbind(ratio_transcribed, ratio_non_transcribed),
    type = "l", col = c("lightskyblue", "gray0"), lty = c(1, 1), lwd = 8,xaxt = "n", ylab = "", yaxt = "n", ylim = yLim,
    cex.lab = 4, cex.axis = 2.5, main = "", font.main = 1, cex.main = 5
  )
  axis(1, at = x_coordinates,labels = x_labels, cex.axis = 4, lwd = 0, line = 2)
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(ylab =  "XR-seq/Damage-seq reads", cex.lab=4, line = 6, xlab = "Relative distance from eTSS", main = "", cex = 6)
  dev.off()
}

make_profile_ratio_dinuc <- function(elementFile, dinucTranscribedFile, dinucNonTranscribedFile, repairTranscribedFile, repairNonTranscribedFile, repairNumOfreads, title, figureName, yLim, x_coordinates, x_labels) {
  # This function generates average ratio profiles of damage & repair on the transcribed and non-transcribed strands.
  elementFile <- readGeneric(elementFile, strand = 6)
  dinucTranscribedFile <- readGeneric(dinucTranscribedFile, strand = 6)
  dinucNonTranscribedFile <- readGeneric(dinucNonTranscribedFile, strand = 6)
  repairTranscribedFile <- readGeneric(repairTranscribedFile, strand = 6)
  repairNonTranscribedFile <- readGeneric(repairNonTranscribedFile, strand = 6)


  sm1 <- ScoreMatrixBin(target = dinucTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 75)
  sm2 <- ScoreMatrixBin(target = dinucNonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 75)
  sm3 <- ScoreMatrixBin(target = repairTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 75)
  sm4 <- ScoreMatrixBin(target = repairNonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = 75)

  # Normalize the data
  sm3 <- sm3/repairNumOfreads
  sm4 <- sm4/repairNumOfreads

  # Divide the mean value of damage in the mean value of repair.
  ratio_transcribed <- colMeans(sm3) / colMeans(sm1)
  ratio_non_transcribed <- colMeans(sm4) / colMeans(sm2)

  sm1_min <- min(ratio_transcribed)
  sm2_min <- min(ratio_non_transcribed)
  
  y_min <- min(sm1_min, sm2_min)
  
  sm1_max <- max(ratio_transcribed)
  sm2_max <- max(ratio_non_transcribed)
  
  y_max <- max(sm1_max, sm2_max)

  y_min <- yLim[1]
  y_max <- yLim[2]
  
  y_coordinates <- seq(y_min, y_max, length.out=3)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point. 

  pdf(paste(figurePath, paste0("dinuc_repair_",figureName), sep = "/"), height = 14, width = 20)
  op <- par(mar = c(10, 15, 10, 10) + 0.1)

  matplot(cbind(ratio_transcribed, ratio_non_transcribed),
    type = "l", col = c("orangered", "gray47"), lty = c(1, 1), lwd = 8, xaxt = "n", ylab = "", ylim = c(y_min, y_max), yaxt = "n",
    cex.lab = 4, cex.axis = 2.5, main = "", font.main = 1, cex.main = 5
  )
  axis(1, at = x_coordinates, labels = x_labels, cex.axis = 4)
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(ylab =  "XR-seq/TT counts", cex.lab=4, line = 6, xlab = "Relative distance from center", main = "", cex = 6)
  dev.off()
}
###################################
#----Average repair profiles------#
###################################
# NHF1
NHF1_plus_strand_right <- paste0(repairPlusPath, "/hg38_enhancers_750_right_NHF1_CPD_4nt.bed")
NHF1_minus_strand_right <- paste0(repairMinusPath, "/hg38_enhancers_750_right_NHF1_CPD_4nt.bed")
make_avg_profile(enhancers_right, NHF1_plus_strand_right, NHF1_minus_strand_right, numOfreads_NHF1_CPD, "NHF1", "hg38_enhancers_plot_NHF1_strands_eTSS_right.pdf", c(8e-10, 3.4e-9), c(0,75), c("eTSS", "750b"))

NHF1_plus_strand_left <- paste0(repairPlusPath, "/hg38_enhancers_750_left_NHF1_CPD_4nt.bed")
NHF1_minus_strand_left <- paste0(repairMinusPath, "/hg38_enhancers_750_left_NHF1_CPD_4nt.bed")
make_avg_profile(enhancers_left, NHF1_plus_strand_left, NHF1_minus_strand_left, numOfreads_NHF1_CPD, "NHF1", "hg38_enhancers_plot_NHF1_strands_eTSS_left.pdf", c(8e-10, 3.4e-9), c(0,75), c("750b", "eTSS"))

# CSB
CSB_plus_strand_right <- paste0(repairPlusPath, "/hg38_enhancers_750_right_CSB_CPD_4nt.bed")
CSB_minus_strand_right <- paste0(repairMinusPath, "/hg38_enhancers_750_right_CSB_CPD_4nt.bed")
make_avg_profile(enhancers_right, CSB_plus_strand_right, CSB_minus_strand_right, numOfreads_CSB_CPD, "CS-B", "hg38_enhancers_plot_CSB_strands_eTSS_right.pdf", c(1e-9, 3.2e-9), c(0,75), c("eTSS", "750b"))

CSB_plus_strand_left <- paste0(repairPlusPath, "/hg38_enhancers_750_left_CSB_CPD_4nt.bed")
CSB_minus_strand_left <- paste0(repairMinusPath, "/hg38_enhancers_750_left_CSB_CPD_4nt.bed")
make_avg_profile(enhancers_left, CSB_plus_strand_left, CSB_minus_strand_left, numOfreads_CSB_CPD, "CS-B", "hg38_enhancers_plot_CSB_strands_eTSS_left.pdf", c(1e-9, 3.2e-9), c(0,75), c("750b", "eTSS"))

# XPC
XPC_plus_strand_right <- paste0(repairPlusPath, "/hg38_enhancers_750_right_XPC_CPD_4nt.bed")
XPC_minus_strand_right <- paste0(repairMinusPath, "/hg38_enhancers_750_right_XPC_CPD_4nt.bed")
make_avg_profile(enhancers_right, XPC_plus_strand_right, XPC_minus_strand_right, numOfreads_XPC_CPD, "XP-C", "hg38_enhancers_plot_XPC_strands_eTSS_right.pdf", c(3.5e-10, 1.8e-9), c(0,75), c("eTSS", "750b"))

XPC_plus_strand_left <- paste0(repairPlusPath, "/hg38_enhancers_750_left_XPC_CPD_4nt.bed")
XPC_minus_strand_left <- paste0(repairMinusPath, "/hg38_enhancers_750_left_XPC_CPD_4nt.bed")
make_avg_profile(enhancers_left, XPC_plus_strand_left, XPC_minus_strand_left, numOfreads_XPC_CPD, "XP-C", "hg38_enhancers_plot_XPC_strands_eTSS_left.pdf", c(3.5e-10,1.8e-9), c(0,75), c("750b", "eTSS"))

 
#######################################################
#----Average repair normalized to damage profiles-----#
#######################################################
 
# Damage-seq files.
CPD_plus_strand_right <- "/vol/sci/bio/data/sheera.adar/elisheva.h/enhancers/eTSS_analysis_data/intersectedData/damage/plus/hg38/hg38_enhancers_750_right_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
CPD_minus_strand_right <- "/vol/sci/bio/data/sheera.adar/elisheva.h/enhancers/eTSS_analysis_data/intersectedData/damage/minus/hg38/hg38_enhancers_750_right_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"

CPD_plus_strand_left <- "/vol/sci/bio/data/sheera.adar/elisheva.h/enhancers/eTSS_analysis_data/intersectedData/damage/plus/hg38/hg38_enhancers_750_left_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
CPD_minus_strand_left <- "/vol/sci/bio/data/sheera.adar/elisheva.h/enhancers/eTSS_analysis_data/intersectedData/damage/minus/hg38/hg38_enhancers_750_left_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"


# NHF1
make_profile_ratio(enhancers_right, CPD_plus_strand_right, CPD_minus_strand_right, NHF1_plus_strand_right, NHF1_minus_strand_right, numOfreads_NHF1_CPD, "NHF1", "hg38_enhancers_plot_NHF1_eTSS_right_750.pdf", c(1.5e-8,5.5e-8), c(0,75), c("eTSS","750b"))
make_profile_ratio(enhancers_left, CPD_plus_strand_left, CPD_minus_strand_left, NHF1_plus_strand_left, NHF1_minus_strand_left, numOfreads_NHF1_CPD, "NHF1", "hg38_enhancers_plot_NHF1_eTSS_left_750.pdf", c(1.5e-8,5.5e-8), c(0,75), c("-750b","eTSS"))

# CSB
make_profile_ratio(enhancers_right, CPD_plus_strand_right, CPD_minus_strand_right, CSB_plus_strand_right, CSB_minus_strand_right, numOfreads_CSB_CPD, "CSB", "hg38_enhancers_plot_CSB_eTSS_right_750.pdf", c(1.5e-8,5e-8), c(0,75), c("eTSS","750b"))
make_profile_ratio(enhancers_left, CPD_plus_strand_left, CPD_minus_strand_left, CSB_plus_strand_left, CSB_minus_strand_left, numOfreads_CSB_CPD, "CSB", "hg38_enhancers_plot_CSB_eTSS_left_750.pdf", c(1.5e-8,5e-8), c(0,75), c("-750b","eTSS"))

# XPC
make_profile_ratio(enhancers_right, CPD_plus_strand_right, CPD_minus_strand_right, XPC_plus_strand_right, XPC_minus_strand_right, numOfreads_XPC_CPD, "XPC", "hg38_enhancers_plot_XPC_eTSS_right_750.pdf", c(5e-9,3.5e-8), c(0,75), c("eTSS","750b"))
make_profile_ratio(enhancers_left, CPD_plus_strand_left, CPD_minus_strand_left, XPC_plus_strand_left, XPC_minus_strand_left, numOfreads_XPC_CPD, "XPC", "hg38_enhancers_plot_XPC_eTSS_left_750.pdf", c(5e-9,3.5e-8), c(0,75), c("-750b","eTSS"))

 
#######################################################
#------Average repair normalized to TT profiles-------#
#######################################################

# TT bed files.
TT_plus_strand_right <- "/vol/sci/bio/data/sheera.adar/elisheva.h/enhancers/eTSS_analysis_data/intersectedData/TT/plus/hg38/hg38_enhancers_750_right_TT.bed"
TT_minus_strand_right <- "/vol/sci/bio/data/sheera.adar/elisheva.h/enhancers/eTSS_analysis_data/intersectedData/TT/minus/hg38/hg38_enhancers_750_right_TT.bed"

TT_plus_strand_left <- "/vol/sci/bio/data/sheera.adar/elisheva.h/enhancers/eTSS_analysis_data/intersectedData/TT/plus/hg38/hg38_enhancers_750_left_TT.bed"
TT_minus_strand_left <- "/vol/sci/bio/data/sheera.adar/elisheva.h/enhancers/eTSS_analysis_data/intersectedData/TT/minus/hg38/hg38_enhancers_750_left_TT.bed"

# NHF1

make_profile_ratio_dinuc(enhancers_right, TT_plus_strand_right, TT_minus_strand_right, NHF1_plus_strand_right, NHF1_minus_strand_right, numOfreads_NHF1_CPD, "NHF1", "hg38_enhancers_plot_NHF1_eTSS_right_750.pdf", c(3e-9,1.4e-8), c(0,75), c("eTSS","750b"))
make_profile_ratio_dinuc(enhancers_left, TT_plus_strand_left, TT_minus_strand_left, NHF1_plus_strand_left, NHF1_minus_strand_left, numOfreads_NHF1_CPD, "NHF1", "hg38_enhancers_plot_NHF1_eTSS_left_750.pdf", c(3e-9,1.4e-8), c(0,75), c("-750b","eTSS"))

# CSB
make_profile_ratio_dinuc(enhancers_right, TT_plus_strand_right, TT_minus_strand_right, CSB_plus_strand_right, CSB_minus_strand_right, numOfreads_CSB_CPD, "CSB", "hg38_enhancers_plot_CSB_eTSS_right_750.pdf", c(3.6e-9,1.3e-8), c(0,75), c("eTSS","750b"))
make_profile_ratio_dinuc(enhancers_left, TT_plus_strand_left, TT_minus_strand_left, CSB_plus_strand_left, CSB_minus_strand_left, numOfreads_CSB_CPD, "CSB", "hg38_enhancers_plot_CSB_eTSS_left_750.pdf", c(3.6e-9,1.3e-8), c(0,75), c("-750b","eTSS"))

# XPC
make_profile_ratio_dinuc(enhancers_right, TT_plus_strand_right, TT_minus_strand_right, XPC_plus_strand_right, XPC_minus_strand_right, numOfreads_XPC_CPD, "XPC", "hg38_enhancers_plot_XPC_eTSS_right_750.pdf", c(1.2e-9,8.4e-9), c(0,75), c("eTSS","750b"))
make_profile_ratio_dinuc(enhancers_left, TT_plus_strand_left, TT_minus_strand_left, XPC_plus_strand_left, XPC_minus_strand_left, numOfreads_XPC_CPD, "XPC", "hg38_enhancers_plot_XPC_eTSS_left_750.pdf", c(1.2e-9,8.4e-9), c(0,75), c("-750b","eTSS"))
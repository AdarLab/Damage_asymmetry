#!/usr/bin/Rscript 
# R version 3.5.2
library(genomation)
library(BSgenome)
library(Biostrings)

##############################################################################
# This script gets genomic coordinates and plot:
# 1) Average profiles of XR-seq reads.
# 2) Average profiles of repair normalized to damage (CPD Damage-seq reads).
# 3) Average profiles of repair normalized to TT counts.
##############################################################################

########################
#-----Num of reads-----#
########################

# Dinucs
numOfreadsTT <- 576715193 # TT in the genome
numOfreadsTC <- 350771128 # TC in the genome
# Damage
numOfreadsNCPD <- 68855889 # CPD dypirimidines
numOfreadsN64 <- 60088813 # (6-4)PP dypirimidines

# Repair
numOfreads_NHF1_CPD <- 11100235 # NHF1 CPD.
numOfreads_NHF1_64PP <- 42121243 # NHF1 (6-4)PP
numOfreads_CSB_CPD <- 19152949 # CSB CPD
numOfreads_CSB_64PP <- 32012948 # CSB (6-4)PP
numOfreads_XPC_CPD <- 14856828 # XPC CPD
numOfreads_XPC_64PP <- 15784721# XPC (6-4)PP

########################
#--------Pathes--------#
########################

dinucPath <- "averageProfilesData/patterns"
damagePath <-  "averageProfilesData/damageIntersect"
repairPath <- "XR_Seq_data/intrsectedDataPlot_4nt"
elementPath <- "averageProfilesData/annotation"
figurePath <- "figures"

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

TC_genes <- "hg38_RefSeq_genes_3Kb_extention_TC.bed"
TC_exons <- "hg38_RefSeq_uniq_chopped_exons_TC.bed"
TC_introns <- "hg38_RefSeq_uniq_chopped_introns_TC.bed"

# Damage files

NCPD_genes_up <- "hg38_RefSeq_uniq_coding_genes_3Kb_up_plot_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
NCPD_genes_down <- "hg38_RefSeq_uniq_coding_genes_3Kb_down_plot_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
NCPD_exons_start <- "hg38_RefSeq_uniq_chopped_exons_start_plot_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
NCPD_exons_end <- "hg38_RefSeq_uniq_chopped_exons_end_plot_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
NCPD_introns_start <- "hg38_RefSeq_uniq_chopped_introns_start_plot_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
NCPD_introns_end<- "hg38_RefSeq_uniq_chopped_introns_end_plot_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"

N64_genes <- "genes_3Kb_extention_N64_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
N64_exons <- "chopped_exons_N64_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"
N64_introns <- "chopped_introns_N64_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed"

# Repair file Names have the same name as the element files.

make_profile <- function(element, cell, type, fileName, elementFile, bins, xLab, coordinates, labels, counts, title, figureName, colors, yLim)
{
  # This function generates average profiles of damage/repair on the transcribed and non-transcribed strand 3 Kb upstream
  # and 10 Kb downstream of the TSS of human genes.
  
  transcribedFile <- readGeneric(paste(repairPath, cell, type, "transcribedStrand", paste(cell, type, fileName, sep = "_"), sep = "/"), strand = 6)
  nonTranscribedFile <- readGeneric(paste(repairPath, cell, type,"nonTranscribedStrand", paste(cell, type, fileName, sep = "_"), sep = "/"), strand = 6)
  elementFile <- readGeneric(paste(elementPath, elementFile, sep = "/"), strand = 6)
  
  sm1 = ScoreMatrixBin(target = transcribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  sm2 = ScoreMatrixBin(target = nonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  
  #Normalize the data
  sm1 <- (sm1/counts)
  sm2 <- (sm2/counts)
  
  plot.data = new("ScoreMatrixList",list(sm1,sm2))
  
  sm1_min <- min(colMeans(sm1))
  sm2_min <- min(colMeans(sm2))
  
  y_min <- min(sm1_min, sm2_min)
  
  sm1_max <- max(colMeans(sm1))
  sm2_max <- max(colMeans(sm2))
  
  y_max <- max(sm1_max, sm2_max)
  
  y_coordinates <- seq(y_min, y_max, length.out=3)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.
  
  pdf(paste0(figurePath, "/repair/", element, "/", figureName, ".pdf"),height=14,width=20)
  op <- par(mar = c(10,15,10,15) + 0.1)
  plotMeta(mat = plot.data, overlay = TRUE, line.col = colors, lwd = 8, xaxt = 'n', yaxt = 'n',ylim = c(y_min, y_max),
           ylab = "", xlab = "", cex.lab = 4, cex.axis = 2.5, font.main = 1)
  axis(1,at = coordinates, labels = labels, cex.axis = 4)
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(main = title, cex.main=4, font.main = 1, ylab = "Average normalized CPD XR-seq read counts\n", cex.lab=4, line = 6, xlab = paste("Relative distance from", xLab))
  dev.off()  
}


make_profile_ratio <- function(element, cell, damagePath, damageFileName, repairPath, repairFileName, elementFile, damage_counts, repair_counts, bins, xLab, coordinates, labels, title, figureName, yLim)
{
  # This function generates average ratio profiles of damage & repair on the transcribed and non-transcribed strands.
  
  damageTranscribedFile <- readGeneric(paste(damagePath, "transcribed", damageFileName, sep = "/"), strand = 6)
  damageNonTranscribedFile <- readGeneric(paste(damagePath, "nonTranscribed", damageFileName, sep = "/"), strand = 6)
  repairTranscribedFile <- readGeneric(paste(repairPath, "transcribedStrand", paste(cell, "CPD", repairFileName, sep = "_"), sep = "/"), strand = 6)
  repairNonTranscribedFile <- readGeneric(paste(repairPath, "nonTranscribedStrand", paste(cell, "CPD", repairFileName, sep = "_"), sep = "/"), strand = 6)
  elementFile <- readGeneric(paste(elementPath, elementFile, sep = "/"), strand = 6)
  
  sm1 = ScoreMatrixBin(target = damageTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  sm2 = ScoreMatrixBin(target = damageNonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  sm3 = ScoreMatrixBin(target = repairTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  sm4 = ScoreMatrixBin(target = repairNonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  
  #Normalize the data
  sm3 <- (sm3/repair_counts)
  sm4 <- (sm4/repair_counts)
  
  # Divide the mean value of damage in the mean value of repair.
  ratio_transcribed <- colMeans(sm3)/colMeans(sm1) 
  ratio_non_transcribed <- colMeans(sm4)/colMeans(sm2)
  
  sm1_min <- min(ratio_transcribed)
  sm2_min <- min(ratio_non_transcribed)
  
  y_min <- min(sm1_min, sm2_min)
  
  sm1_max <- max(ratio_transcribed)
  sm2_max <- max(ratio_non_transcribed)
  
  y_max <- max(sm1_max, sm2_max)
  
  y_coordinates <- seq(y_min, y_max, length.out=3)
  y_coordinates <- signif(y_coordinates, digits = 2) # Round the numbers to get only 2 digits after decimal point.
  
  pdf(paste0(figurePath, "/damage_repair/", element, "/", figureName, ".pdf"),height=14, width=20)
  op <- par(mar = c(10,15,10,10) + 0.1)

  matplot(cbind(ratio_transcribed, ratio_non_transcribed),type="l",col = c("darkturquoise","darkorchid4"),lty = c(1,1), lwd = 8,
          xaxt = 'n', yaxt = "n", ylab = "", ylim = c(y_min, y_max),
          cex.lab = 4, cex.axis = 2.5, main = "", font.main = 1, cex.main = 5)
  axis(1,at = coordinates, labels = labels, cex.axis = 4)
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(main = title, cex.main=4, font.main = 1, ylab = "XR-Seq reads/Damage-Seq reads", cex.lab=4, line = 6, xlab = paste("Relative distance from", xLab))
  dev.off()
}

make_profile_ratio_dinuc <- function(element, cell, dinucPath, dinucFileName, repairPath, repairFileName, elementFile, dinuc_counts, repair_counts, bins, xLab, coordinates, labels, title, figureName, yLim)
{
  # This function generates average ratio profiles of dinuc & repair on the transcribed and non-transcribed strands.
  
  dinucTranscribedFile <- readGeneric(paste(dinucPath, "transcribedStrand", dinucFileName, sep = "/"), strand = 6)
  dinucNonTranscribedFile <- readGeneric(paste(dinucPath, "nonTranscribedStrand", dinucFileName, sep = "/"), strand = 6)
  repairTranscribedFile <- readGeneric(paste(repairPath, "transcribedStrand", paste(cell, "CPD", repairFileName, sep = "_"), sep = "/"), strand = 6)
  repairNonTranscribedFile <- readGeneric(paste(repairPath, "nonTranscribedStrand", paste(cell, "CPD", repairFileName, sep = "_"), sep = "/"), strand = 6)
  elementFile <- readGeneric(paste(elementPath, elementFile, sep = "/"), strand = 6)
  
  sm1 = ScoreMatrixBin(target = dinucTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  sm2 = ScoreMatrixBin(target = dinucNonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  sm3 = ScoreMatrixBin(target = repairTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  sm4 = ScoreMatrixBin(target = repairNonTranscribedFile, windows = elementFile, strand.aware = TRUE, bin.num = bins)
  
#  #Normalize the data
#  sm1 <- (sm1/dinuc_counts)*10^9
#  sm2 <- (sm2/dinuc_counts)*10^9
  sm3 <- (sm3/repair_counts)
  sm4 <- (sm4/repair_counts)
  
  # Divide the mean value of dinuc in the mean value of repair.
  ratio_transcribed <- colMeans(sm3)/colMeans(sm1) 
  ratio_non_transcribed <- colMeans(sm4)/colMeans(sm2)
  
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
  
  pdf(paste0(figurePath, "/dinuc_repair/", element, "/", figureName, ".pdf"),height=14, width=20)
  op <- par(mar = c(10,15,10,10) + 0.1)
  
  matplot(cbind(ratio_transcribed, ratio_non_transcribed),type="l",col = c("springgreen4","black"),lty = c(1,1), lwd = 8,
          xaxt = 'n', yaxt = "n", ylab = "", ylim = c(y_min, y_max),
          cex.lab = 4, cex.axis = 2.5, main = "", font.main = 1, cex.main = 5)
  axis(1,at = coordinates, labels = labels, cex.axis = 4)
  axis(2, at = y_coordinates, labels = y_coordinates, cex.axis = 4)
  title(main = title, cex.main=4, font.main = 1, ylab = "XR-Seq reads/TT counts", cex.lab=4, line = 6, xlab = paste("Relative distance from", xLab))
  dev.off()
}

###################################
#----Average repair profiles------#
###################################

#-------------------------------------------#
#-------------------Genes-------------------#
#-------------------------------------------#
# NHF1

make_profile("genes", "NHF1", "CPD", genes_3Kb_up, genes_3Kb_up, 325, "TSS", c(0,75,325), c("-3Kb", "TSS", "10Kb"), numOfreads_NHF1_CPD, "NHF1 CPD XR-Seq reads", "CPD_repair_NHF1_genes_up_profile", c("darkred","darkorange1"), c(0.5e-09,6.5e-09))

make_profile("genes", "NHF1", "CPD", genes_3Kb_down, genes_3Kb_down, 325,  "TSS",  c(0,250,325), c("-10Kb", "TES", "3Kb"), numOfreads_NHF1_CPD, "NHF1 CPD XR-Seq reads", "CPD_repair_NHF1_genes_down_profile", c("darkred","darkorange1"), c(0.5e-09,6.5e-09))

#CSB
make_profile("genes", "CSB", "CPD", genes_3Kb_up, genes_3Kb_up, 325,  "TSS",  c(0,75,325), c("-3Kb", "TSS", "10Kb"), numOfreads_CSB_CPD, "CSB CPD XR-Seq reads", "CPD_repair_CSB_genes_up_profile", c("darkred","darkorange1"), c(1.2e-09,6.5e-09))

make_profile("genes", "CSB", "CPD", genes_3Kb_down, genes_3Kb_down, 325, "TSS",  c(0,250,325), c("-10Kb", "TES", "3Kb"), numOfreads_CSB_CPD, "CSB CPD XR-Seq reads", "CPD_repair_CSB_genes_down_profile", c("darkred","darkorange1"), c(1.2e-09,6.5e-09))
# XPC
make_profile("genes", "XPC", "CPD", genes_3Kb_up, genes_3Kb_up, 325,  "TSS", c(0,75,325), c("-3Kb", "TSS", "10Kb"), numOfreads_XPC_CPD, "XPC CPD XR-Seq reads", "CPD_repair_XPC_genes_up_profile", c("darkred","darkorange1"), genesLimitis)

make_profile("genes", "XPC", "CPD", genes_3Kb_down, genes_3Kb_down, 325,  "TSS", c(0,250,325), c("-10Kb", "TES", "3Kb"), numOfreads_XPC_CPD, "XPC CPD XR-Seq reads", "CPD_repair_XPC_genes_down_profile", c("darkred","darkorange1"), genesLimitis)

#-------------------------------------------#
#-------------------Exons-------------------#
#-------------------------------------------#
genesLimitis <- c(0,10.5)
# NHF1
make_profile("exons", "NHF1", "CPD", exons_start, exons_start, 120, "Exon start", c(0,120), c("Exon start", "120b"), numOfreads_NHF1_CPD, "NHF1 CPD XR-Seq reads", "CPD_repair_NHF1_exons_start_profile", c("darkred","darkorange1"), genesLimitis)

make_profile("exons", "NHF1", "CPD", exons_end, exons_end, 120,  "Exon end",  c(0,120), c("-120b", "Exon end"), numOfreads_NHF1_CPD, "NHF1 CPD XR-Seq reads", "CPD_repair_NHF1_exons_end_profile", c("darkred","darkorange1"), genesLimitis)

# CSB
make_profile("exons", "CSB", "CPD", exons_start, exons_start, 120,  "Exon start", c(0,120), c("Exon start", "120b"), numOfreads_NHF1_CPD, "CSB CPD XR-Seq reads", "CPD_repair_CSB_exons_start_profile", c("darkred","darkorange1"), genesLimitis)

make_profile("exons", "CSB", "CPD", exons_end, exons_end, 120, "Exon end",  c(0,120), c("-120b", "Exon end"), numOfreads_NHF1_CPD, "CSB CPD XR-Seq reads", "CPD_repair_CSB_exons_end_profile", c("darkred","darkorange1"), genesLimitis)
# XPC
make_profile("exons", "XPC", "CPD", exons_start, exons_start, 120,  "Exon start", c(0,120), c("Exon start", "120b"), numOfreads_NHF1_CPD, "XPC CPD XR-Seq reads", "CPD_repair_XPC_exons_start_profile", c("darkred","darkorange1"), genesLimitis)

make_profile("exons", "XPC", "CPD", exons_end, exons_end, 120,  "Exon end", c(0,120), c("-120b", "Exon end"), numOfreads_NHF1_CPD, "XPC CPD XR-Seq reads", "CPD_repair_XPC_exons_end_profile", c("darkred","darkorange1"), genesLimitis)

#-------------------------------------------#
#-------------------Introns-----------------#
#-------------------------------------------#
# NHF1
make_profile("introns", "NHF1", "CPD", introns_start, introns_start, 180, "Intron start", c(0,180), c("Intron start", "1.86Kb"), numOfreads_NHF1_CPD, "NHF1 CPD XR-Seq reads", "CPD_repair_NHF1_introns_start_profile", c("darkred","darkorange1"), genesLimitis)

make_profile("introns", "NHF1", "CPD", introns_end, introns_end, 180,  "Intron end",  c(0,180), c("-1.86Kb", "Intron end"), numOfreads_NHF1_CPD, "NHF1 CPD XR-Seq reads", "CPD_repair_NHF1_introns_end_profile", c("darkred","darkorange1"), genesLimitis)

# CSB
make_profile("introns", "CSB", "CPD", introns_start, introns_start, 180,  "Intron start", c(0,180), c("Intron start", "1.86Kb"), numOfreads_NHF1_CPD, "CSB CPD XR-Seq reads", "CPD_repair_CSB_introns_start_profile", c("darkred","darkorange1"), genesLimitis)

make_profile("introns", "CSB", "CPD", introns_end, introns_end, 180, "Intron end",  c(0,180), c("-1.86Kb", "Intron end"), numOfreads_NHF1_CPD, "CSB CPD XR-Seq reads", "CPD_repair_CSB_introns_end_profile", c("darkred","darkorange1"), genesLimitis)
# XPC
make_profile("introns", "XPC", "CPD", introns_start, introns_start, 180,  "Intron start", c(0,180), c("Intron start", "1.86Kb"), numOfreads_NHF1_CPD, "XPC CPD XR-Seq reads", "CPD_repair_XPC_introns_start_profile", c("darkred","darkorange1"), genesLimitis)

make_profile("introns", "XPC", "CPD", introns_end, introns_end, 180,  "Intron end", c(0,180), c("-1.86Kb", "Intron end"), numOfreads_NHF1_CPD, "XPC CPD XR-Seq reads", "CPD_repair_XPC_introns_end_profile", c("darkred","darkorange1"), genesLimitis)


#######################################################
#----Average repair normalized to damage profiles-----#
#######################################################
#-------------------------------------------#
#-------------------Genes-------------------#
#-------------------------------------------#

yLim <- c(0,25)

# NHF1
make_profile_ratio("genes", "NHF1", damagePath, paste0("NCPD/", NCPD_genes_up), paste0(repairPath, "/NHF1/CPD"), genes_3Kb_up, genes_3Kb_up, numOfreadsNCPD, numOfreads_NHF1_CPD, 325, "TSS", c(0,75,325), c("-3Kb", "TSS", "10Kb"), "NHF1 CPD repair/damage ratio", "CPD_damage_repair_ratio_NHF1_genes_up_profile", yLim)
make_profile_ratio("genes", "NHF1", damagePath, paste0("NCPD/", NCPD_genes_down), paste0(repairPath, "/NHF1/CPD"), genes_3Kb_down, genes_3Kb_down, numOfreadsNCPD, numOfreads_NHF1_CPD, 325, "TES", c(0,250,325), c("-10Kb", "TES", "3Kb"), "NHF1 CPD repair/damage ratio", "CPD_damage_repair_ratio_NHF1_genes_down_profile", yLim)

# CSB
make_profile_ratio("genes", "CSB", damagePath, paste0("NCPD/", NCPD_genes_up), paste0(repairPath, "/CSB/CPD"), genes_3Kb_up, genes_3Kb_up, numOfreadsNCPD, numOfreads_CSB_CPD, 325, "TSS", c(0,75,325), c("-3Kb", "TSS", "10Kb"), "CSB CPD repair/damage ratio", "CPD_damage_repair_ratio_CSB_genes_up_profile", yLim)
make_profile_ratio("genes", "CSB", damagePath, paste0("NCPD/", NCPD_genes_down), paste0(repairPath, "/CSB/CPD"), genes_3Kb_down, genes_3Kb_down, numOfreadsNCPD, numOfreads_CSB_CPD, 325, "TES", c(0,250,325), c("-10Kb", "TES", "3Kb"), "CSB CPD repair/damage ratio", "CPD_damage_repair_ratio_CSB_genes_down_profile", yLim)

# XPC
make_profile_ratio("genes", "XPC", damagePath, paste0("NCPD/", NCPD_genes_up), paste0(repairPath, "/XPC/CPD"), genes_3Kb_up, genes_3Kb_up, numOfreadsNCPD, numOfreads_XPC_CPD, 325, "TSS", c(0,75,325), c("-3Kb", "TSS", "10Kb"), "XPC CPD repair/damage ratio", "CPD_damage_repair_ratio_XPC_genes_up_profile", yLim)
make_profile_ratio("genes", "XPC", damagePath, paste0("NCPD/", NCPD_genes_down), paste0(repairPath, "/XPC/CPD"), genes_3Kb_down, genes_3Kb_down, numOfreadsNCPD, numOfreads_XPC_CPD, 325, "TES", c(0,250,325), c("-10Kb", "TES", "3Kb"), "XPC CPD repair/damage ratio", "CPD_damage_repair_ratio_XPC_genes_down_profile", yLim)

#-------------------------------------------#
#-------------------Exons-------------------#
#-------------------------------------------#
yLim <- c(0,12)
# NHF1
make_profile_ratio("exons", "NHF1", damagePath, paste0("NCPD/", NCPD_exons_start), paste0(repairPath, "/NHF1/CPD"), exons_start, exons_start, numOfreadsNCPD, numOfreads_NHF1_CPD, 120, "Exons start", c(0,120), c("Exons start", "120b"), "NHF1 CPD repair/damage ratio", "CPD_damage_repair_ratio_NHF1_exons_start_profile", yLim)
make_profile_ratio("exons", "NHF1", damagePath, paste0("NCPD/", NCPD_exons_end), paste0(repairPath, "/NHF1/CPD"), exons_end, exons_end, numOfreadsNCPD, numOfreads_NHF1_CPD, 120, "Exon end", c(0,120), c("-120b", "Exon end"), "NHF1 CPD repair/damage ratio", "CPD_damage_repair_ratio_NHF1_exons_end_profile", yLim)

# CSB
make_profile_ratio("exons", "CSB", damagePath, paste0("NCPD/", NCPD_exons_start), paste0(repairPath, "/CSB/CPD"), exons_start, exons_start, numOfreadsNCPD, numOfreads_CSB_CPD, 120, "Exons start", c(0,120), c("Exons start", "120b"), "CSB CPD repair/damage ratio", "CPD_damage_repair_ratio_CSB_exons_start_profile", yLim)
make_profile_ratio("exons", "CSB", damagePath, paste0("NCPD/", NCPD_exons_end), paste0(repairPath, "/CSB/CPD"), exons_end, exons_end, numOfreadsNCPD, numOfreads_CSB_CPD, 120, "Exon end", c(0,120), c("-120b", "Exon end"), "CSB CPD repair/damage ratio", "CPD_damage_repair_ratio_CSB_exons_end_profile", yLim)

# XPC
make_profile_ratio("exons", "XPC", damagePath, paste0("NCPD/", NCPD_exons_start), paste0(repairPath, "/XPC/CPD"), exons_start, exons_start, numOfreadsNCPD, numOfreads_XPC_CPD, 120, "Exons start", c(0,120), c("Exons start", "120b"), "XPC CPD repair/damage ratio", "CPD_damage_repair_ratio_XPC_exons_start_profile", yLim)
make_profile_ratio("exons", "XPC", damagePath, paste0("NCPD/", NCPD_exons_end), paste0(repairPath, "/XPC/CPD"), exons_end, exons_end, numOfreadsNCPD, numOfreads_XPC_CPD, 120, "Exon end", c(0,120), c("-120b", "Exon end"), "XPC CPD repair/damage ratio", "CPD_damage_repair_ratio_XPC_exons_end_profile", yLim)

#-------------------------------------------#
#------------------Introns------------------#
#-------------------------------------------#
# NHF1
make_profile_ratio("introns", "NHF1", damagePath, paste0("NCPD/", NCPD_introns_start), paste0(repairPath, "/NHF1/CPD"), introns_start, introns_start, numOfreadsNCPD, numOfreads_NHF1_CPD, 180, "Intron start", c(0,180), c("Intron start", "1.86Kb"), "NHF1 CPD repair/damage ratio", "CPD_damage_repair_ratio_NHF1_introns_start_profile", yLim)
make_profile_ratio("introns", "NHF1", damagePath, paste0("NCPD/", NCPD_introns_end), paste0(repairPath, "/NHF1/CPD"), introns_end, introns_end, numOfreadsNCPD, numOfreads_NHF1_CPD, 180, "Intron end", c(0,180), c("-1.86Kb", "Intron end"), "NHF1 CPD repair/damage ratio", "CPD_damage_repair_ratio_NHF1_introns_end_profile", yLim)

# CSB
make_profile_ratio("introns", "CSB", damagePath, paste0("NCPD/", NCPD_introns_start), paste0(repairPath, "/CSB/CPD"), introns_start, introns_start, numOfreadsNCPD, numOfreads_CSB_CPD, 180, "Intron start", c(0,180), c("Intron start", "1.86Kb"), "CSB CPD repair/damage ratio", "CPD_damage_repair_ratio_CSB_introns_start_profile", yLim)
make_profile_ratio("introns", "CSB", damagePath, paste0("NCPD/", NCPD_introns_end), paste0(repairPath, "/CSB/CPD"), introns_end, introns_end, numOfreadsNCPD, numOfreads_CSB_CPD, 180, "Intron end", c(0,180), c("-1.86Kb", "Intron end"), "CSB CPD repair/damage ratio", "CPD_damage_repair_ratio_CSB_introns_end_profile", yLim)

# XPC
make_profile_ratio("introns", "XPC", damagePath, paste0("NCPD/", NCPD_introns_start), paste0(repairPath, "/XPC/CPD"), introns_start, introns_start, numOfreadsNCPD, numOfreads_XPC_CPD, 180, "Intron start", c(0,180), c("Intron start", "1.86Kb"), "XPC CPD repair/damage ratio", "CPD_damage_repair_ratio_XPC_introns_start_profile", yLim)
make_profile_ratio("introns", "XPC", damagePath, paste0("NCPD/", NCPD_introns_end), paste0(repairPath, "/XPC/CPD"), introns_end, introns_end, numOfreadsNCPD, numOfreads_XPC_CPD, 180, "Intron end", c(0,180), c("-1.86Kb", "Intron end"), "XPC CPD repair/damage ratio", "CPD_damage_repair_ratio_XPC_introns_end_profile", yLim)


#######################################################
#------Average repair normalized to TT profiles-------#
#######################################################
#-------------------------------------------#
#-------------------Genes-------------------#
#-------------------------------------------#
yLim <- c(0,50)
# NHF1
make_profile_ratio_dinuc("genes", "NHF1", dinucPath, TT_genes_up, paste0(repairPath, "/NHF1/CPD"), genes_3Kb_up, genes_3Kb_up, numOfreadsTT, numOfreads_NHF1_CPD, 325, "TSS", c(0,75,325), c("-3Kb", "TSS", "10Kb"), "NHF1 CPD repair/TT ratio", "CPD_dinuc_repair_ratio_NHF1_genes_up_profile", yLim)
make_profile_ratio_dinuc("genes", "NHF1", dinucPath, TT_genes_down, paste0(repairPath, "/NHF1/CPD"), genes_3Kb_down, genes_3Kb_down, numOfreadsTT, numOfreads_NHF1_CPD, 325, "TES", c(0,250,325), c("-10Kb", "TES", "3Kb"), "NHF1 CPD repair/TT ratio", "CPD_dinuc_repair_ratio_NHF1_genes_down_profile", yLim)
#
 CSB
make_profile_ratio_dinuc("genes", "CSB", dinucPath, TT_genes_up, paste0(repairPath, "/CSB/CPD"), genes_3Kb_up, genes_3Kb_up, numOfreadsTT, numOfreads_CSB_CPD, 325, "TSS", c(0,75,325), c("-3Kb", "TSS", "10Kb"), "CSB CPD repair/TT ratio", "CPD_dinuc_repair_ratio_CSB_genes_up_profile", yLim)
make_profile_ratio_dinuc("genes", "CSB", dinucPath, TT_genes_down, paste0(repairPath, "/CSB/CPD"), genes_3Kb_down, genes_3Kb_down, numOfreadsTT, numOfreads_CSB_CPD, 325, "TES", c(0,250,325), c("-10Kb", "TES", "3Kb"), "CSB CPD repair/TT ratio", "CPD_dinuc_repair_ratio_CSB_genes_down_profile", yLim)

# XPC
make_profile_ratio_dinuc("genes", "XPC", dinucPath, TT_genes_up, paste0(repairPath, "/XPC/CPD"), genes_3Kb_up, genes_3Kb_up, numOfreadsTT, numOfreads_XPC_CPD, 325, "TSS", c(0,75,325), c("-3Kb", "TSS", "10Kb"), "XPC CPD repair/TT ratio", "CPD_dinuc_repair_ratio_XPC_genes_up_profile", yLim)
make_profile_ratio_dinuc("genes", "XPC", dinucPath, TT_genes_down, paste0(repairPath, "/XPC/CPD"), genes_3Kb_down, genes_3Kb_down, numOfreadsTT, numOfreads_XPC_CPD, 325, "TES", c(0,250,325), c("-10Kb", "TES", "3Kb"), "XPC CPD repair/TT ratio", "CPD_dinuc_repair_ratio_XPC_genes_down_profile", yLim)

#-------------------------------------------#
#-------------------Exons-------------------#
#-------------------------------------------#
yLim <- c(2.6e-09,1.7e-08)
# NHF1
make_profile_ratio_dinuc("exons", "NHF1", dinucPath, TT_exons_start, paste0(repairPath, "/NHF1/CPD"), exons_start, exons_start, numOfreadsTT, numOfreads_NHF1_CPD, 120, "Exons start", c(0,120), c("Exons start", "120b"), "NHF1 CPD repair/TT ratio", "CPD_dinuc_repair_ratio_NHF1_exons_start_profile", yLim)
make_profile_ratio_dinuc("exons", "NHF1", dinucPath, TT_exons_end, paste0(repairPath, "/NHF1/CPD"), exons_end, exons_end, numOfreadsTT, numOfreads_NHF1_CPD, 120, "Exon end", c(0,120), c("-120b", "Exon end"), "NHF1 CPD repair/TT ratio", "CPD_dinuc_repair_ratio_NHF1_exons_end_profile", yLim)

# CSB
make_profile_ratio_dinuc("exons", "CSB", dinucPath, TT_exons_start, paste0(repairPath, "/CSB/CPD"), exons_start, exons_start, numOfreadsTT, numOfreads_CSB_CPD, 120, "Exons start", c(0,120), c("Exons start", "120b"), "CSB CPD repair/TT ratio", "CPD_dinuc_repair_ratio_CSB_exons_start_profile", c(3.2e-9, 1e-8))
make_profile_ratio_dinuc("exons", "CSB", dinucPath, TT_exons_end, paste0(repairPath, "/CSB/CPD"), exons_end, exons_end, numOfreadsTT, numOfreads_CSB_CPD, 120, "Exon end", c(0,120), c("-120b", "Exon end"), "CSB CPD repair/TT ratio", "CPD_dinuc_repair_ratio_CSB_exons_end_profile", yLim)

# XPC
make_profile_ratio_dinuc("exons", "XPC", dinucPath, TT_exons_start, paste0(repairPath, "/XPC/CPD"), exons_start, exons_start, numOfreadsTT, numOfreads_XPC_CPD, 120, "Exons start", c(0,120), c("Exons start", "120b"), "XPC CPD repair/TT ratio", "CPD_dinuc_repair_ratio_XPC_exons_start_profile", c(6.8e-10, 3.5e-8))
make_profile_ratio_dinuc("exons", "XPC", dinucPath, TT_exons_end, paste0(repairPath, "/XPC/CPD"), exons_end, exons_end, numOfreadsTT, numOfreads_XPC_CPD, 120, "Exon end", c(0,120), c("-120b", "Exon end"), "XPC CPD repair/TT ratio", "CPD_dinuc_repair_ratio_XPC_exons_end_profile", yLim)

#-------------------------------------------#
#------------------Introns------------------#
#-------------------------------------------#
# NHF1
make_profile_ratio_dinuc("introns", "NHF1", dinucPath, TT_introns_start, paste0(repairPath, "/NHF1/CPD"), introns_start, introns_start, numOfreadsTT, numOfreads_NHF1_CPD, 180, "Intron start", c(0,180), c("Intron start", "1.86Kb"), "NHF1 CPD repair/TT ratio", "CPD_dinuc_repair_ratio_NHF1_introns_start_profile", yLim)
make_profile_ratio_dinuc("introns", "NHF1", dinucPath, TT_introns_end, paste0(repairPath, "/NHF1/CPD"), introns_end, introns_end, numOfreadsTT, numOfreads_NHF1_CPD, 180, "Intron end", c(0,180), c("-1.86Kb", "Intron end"), "NHF1 CPD repair/TT ratio", "CPD_dinuc_repair_ratio_NHF1_introns_end_profile", yLim)

# CSB
make_profile_ratio_dinuc("introns", "CSB", dinucPath, TT_introns_start, paste0(repairPath, "/CSB/CPD"), introns_start, introns_start, numOfreadsTT, numOfreads_CSB_CPD, 180, "Intron start", c(0,180), c("Intron start", "1.84Kb"), "CSB CPD repair/TT ratio", "CPD_dinuc_repair_ratio_CSB_introns_start_profile", c(3.2e-9, 1e-8))
make_profile_ratio_dinuc("introns", "CSB", dinucPath, TT_introns_end, paste0(repairPath, "/CSB/CPD"), introns_end, introns_end, numOfreadsTT, numOfreads_CSB_CPD, 180, "Intron end", c(0,180), c("-1.86Kb", "Intron end"), "CSB CPD repair/TT ratio", "CPD_dinuc_repair_ratio_CSB_introns_end_profile", yLim)

# XPC
make_profile_ratio_dinuc("introns", "XPC", dinucPath, TT_introns_start, paste0(repairPath, "/XPC/CPD"), introns_start, introns_start, numOfreadsTT, numOfreads_XPC_CPD, 180, "Intron start", c(0,180), c("Intron start", "1.84Kb"), "XPC CPD repair/TT ratio", "CPD_dinuc_repair_ratio_XPC_introns_start_profile", c(6.8e-10, 3.5e-8))
make_profile_ratio_dinuc("introns", "XPC", dinucPath, TT_introns_end, paste0(repairPath, "/XPC/CPD"), introns_end, introns_end, numOfreadsTT, numOfreads_XPC_CPD, 180, "Intron end", c(0,180), c("-1.86Kb", "Intron end"), "XPC CPD repair/TT ratio", "CPD_dinuc_repair_ratio_XPC_introns_end_profile", yLim)

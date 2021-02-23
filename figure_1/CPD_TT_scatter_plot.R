library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(ggpointdensity)
warnings()

rm(list = ls())
# Paths
current_path <<- "your_path"
genes_counting_path <<- paste(current_path, "tables/counting/genes_seq_counting.csv", sep = "/")
damages_frequencies_path <<- paste(current_path, "damageCoverage", sep = "/")
figure_path <<- "your_figure_path"

creatData <- function(damageFile, dinuc_trans, dinuc_non_trans) {
  # This function gets dinucleotide counting table and CPD counts over genes and combines the given nucleotide counts with CPD frequencies.
  fileName <- gsub(".bed", "", damageFile)
  # Extract the dinucleotide counts and normalize them to 1 Kb.
  counting_table <- read.csv(genes_counting_path, header = T)
  counting_table_dinuc <- cbind(as.data.frame(counting_table$Id), counting_table$length, ((counting_table[dinuc_non_trans] + counting_table[dinuc_trans]) / counting_table$length) * 1000, counting_table[dinuc_non_trans] / counting_table$length, counting_table[dinuc_trans] / counting_table$length)
  names(counting_table_dinuc) <- c("Id", "length", dinuc_trans, paste(dinuc_trans, "transcribed", sep = "_"), paste(dinuc_trans, "non_transcribed", sep = "_")) # Join the freuencies of the given dinucleotide on the transcribed strand and the non-transcribed strand.
  counting_table_dinuc$Id <- gsub("::.*", "", counting_table_dinuc$Id)

  # Damage frequencies
  damage_file_name <- damageFile
  damage_table_trans <- read.table(paste(damages_frequencies_path, "transcribed", damageFile, sep = "/"), sep = "\t")
  damage_table_non_trans <- read.table(paste(damages_frequencies_path, "nonTranscribed", damageFile, sep = "/"), sep = "\t")
  # Extract the CPD counts and normalize them to 1 Kb.
  # $V4 - The Id, $V7 - the CPD counts, $V3-$V2 - the element length.
  damage_frequency_table <- cbind(as.data.frame(damage_table_trans$V4), ((damage_table_trans$V7 + damage_table_non_trans$V7) / (damage_table_trans$V3 - damage_table_trans$V2)) * 1000, damage_table_trans$V8, damage_table_non_trans$V8 # Join the CPD frequencies on the transcribed strand and the non-transcribed strand.
  names(damage_frequency_table) <- c("Id", "damage", "damage_transcribed", "damage_non_transcribed") 
  damage_frequency_table$Id <- gsub("::.*", "", damage_frequency_table$Id)
  
  # Combined table
  dinuc_damage_table <- merge(counting_table_dinuc, damage_frequency_table, by = "Id") # Merge the dinucs and the damage tables by the Id.

  return(dinuc_damage_table)
}

NCPD_TT_table <- creatData("uniq_coding_genes_NCPD_T0_RepA_RepB_sorted_uniq_dipyrimidine.bed" "TT", "AA")

pdf(paste(figure_path, "CPD_TT_genes_scatter_plot.pdf", sep = "/"), height = 6, width = 11)
ggscatter(NCPD_TT_table, x = "TT", y = "damage", color = "black", shape = 21, size = 3, add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "spearman", xlab = "TT frequency (per Kb)", ylab = "CPD levels (per Kb)",  font.label = list(size = 7)  ,add.params = list(color = "#D11515", fill = "lightgray")) +  geom_pointdensity() + scale_color_viridis_c()
dev.off()
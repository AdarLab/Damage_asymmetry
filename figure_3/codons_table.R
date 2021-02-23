#!/usr/bin/Rscript 
# R version 3.5.2
library(seqinr)
library(ggplot2)
############################################################################################
# This script gets CDS sequence and creates tables of codon counts and adjacent codon counts
############################################################################################
codons <- c(
  "ATT", "ATC", "ATA", "CTT", "CTC", "CTA", "CTG", "TTA", "TTG", "GTT", "GTC", "GTA", "GTG", "TTT", "TTC", "ATG", "TGT", "TGC",
  "GCT", "GCC", "GCA", "GCG", "GGT", "GGC", "GGA", "GGG", "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC",
  "TAT", "TAC", "TGG", "CAA", "CAG", "AAT", "AAC", "CAT", "CAC", "GAA", "GAG", "GAT", "GAC", "AAA", "AAG", "CGT", "CGC", "CGA", "CGG", "AGA", "AGG", "TAA", "TAG", "TGA"
)

stopCodons <- c("TAA", "TAG", "TGA")
nucs <<- c("A", "T", "C", "G")

codonsTable <- as.data.frame(matrix(0L, nrow = length(codons), ncol = 1))
colnames(codonsTable) <- "frequency"
rownames(codonsTable) <- codons

adjacentCodonsTable <- as.data.frame(matrix(0L, nrow = length(codons), ncol = length(codons)))
colnames(adjacentCodonsTable) <- codons
rownames(adjacentCodonsTable) <- codons
numOfCodons <- 0
extract_codons <- function(sequence, adjacentCodonsTable, numOfCodons, i) {
  "This function builds a table of 2 adjoining codons frequencies,
  the rows of the table represent the frequency of a the first and the columns represent the frequency of the second codon"
  flag <- F
  startCodon <- paste(sequence[1:3], collapse = "") # Extract the first codon of the sequence.
  stopCodon <- paste(sequence[length(sequence) - 2:0], collapse = "") # Extract the last codon of the sequence.
  print(i)
  # Check if:
  # 1) The first codon is Methionine.
  # 2) The last codon is stop codon.
  # 3) The sequence length is divided by 3 without remainder.

  if (startCodon == "ATG" & (any(stopCodon == stopCodons)) & (length(sequence) %% 3 == 0) & all(sequence %in% nucs)) {
    flag <- T

    # Scan all the codons of the sequence and check if there is a stop codon in the middle of the sequence.
    for (i in seq(1, length(sequence) - 3, 3)) {
      # In case there is stop a codon in the middle of the sequence, set the bolean flag to False.
      if (any(paste(sequence[i:(i + 2)], collapse = "") == stopCodons)) {
        flag <- F        
      }
  
    }
  }

    if (flag) {
      numOfCodons <- numOfCodons + (length(sequence) / 3) # Add the amount of codons in the current sequence.
       for (i in seq(1, length(sequence), 3)) {
         currentCodon <- paste(sequence[i:(i + 2)], collapse = "")
         codonsTable[currentCodon, "frequency"] <- codonsTable[currentCodon, "frequency"] + 1 # Add the occurrence of the current codon to the table.
         assign("codonsTable", codonsTable, envir = .GlobalEnv) # Update the table.
       }
      for (i in seq(1, length(sequence) - 3, 3)) {
        firstCodon <- paste(sequence[i:(i + 2)], collapse = "") # Extract the first of the 2 adjoining codons.
        secondCodon <- paste(sequence[(i + 3):(i + 5)], collapse = "") # Extract the second codon.
        # Check if there is any stop codon inside the sequence.
        print(firstCodon)
        print(secondCodon)
        adjacentCodonsTable[firstCodon, secondCodon] <- adjacentCodonsTable[firstCodon, secondCodon] + 1 # Add the occurrence of the 2 adjoining codons to the table.
        assign("adjacentCodonsTable", adjacentCodonsTable, envir = .GlobalEnv) # Update the table.
      }
    }


  return(numOfCodons)
}

# #####################
# #------ Human ------#
# #####################
 path <- "your_path"
 codonsFile <- "codonsTable.csv"
 adjacentCodonsFile <- "adjacentCodonsTable.csv"
 numOfCodons <- 0
 cdsFile <- read.fasta(paste(path, "RefSeq_hg38_CDS_13_01_20.fa", sep = "/"), forceDNAtolower = F, seqonly = T)
 for (i in 1:length(cdsFile)) {
   # Send evey CDS to the extractCodons function
   numOfCodons <- extract_codons(strsplit(cdsFile[[i]], "")[[1]], adjacentCodonsTable, numOfCodons)
 }
 
 write.csv(codonsTable, file = codonsFile)
 write.csv(adjacentCodonsTable, file = adjacentCodonsFile)

#####################
#--- Arabidopsis ---#
#####################

path <- "Arabidopsis_path"
codonsFile <- "Arabidopsis/codonsTable.csv"
adjacentCodonsFile <- "Arabidopsis/adjacentCodonsTable.csv"
numOfCodons <- 0
cdsFile <- read.fasta(paste(path, "Arabidopsis_thaliana_cds_ensembl_release_47_01_06_2020.fa", sep = "/"), forceDNAtolower = F, seqonly = T)
for (i in 1:length(cdsFile)) {
  # Send evey CDS to the extractCodons function
  numOfCodons <- extract_codons(strsplit(cdsFile[[i]], "")[[1]], adjacentCodonsTable, numOfCodons, i)
}

print(numOfCodons)
write.csv(codonsTable, file = codonsFile)
write.csv(adjacentCodonsTable, file = adjacentCodonsFile)

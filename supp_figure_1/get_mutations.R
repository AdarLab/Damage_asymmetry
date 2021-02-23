

# Directories -----------------------------------------------------------

setwd(file.path("~","projects","paperEnhancers"))
projDir <- getwd()
resDir <-  file.path(projDir,"results")
dataDir <- file.path(projDir,"data")



# Load packages -----------------------------------------------------------

library("dplyr")
library("readr")
library("VariantAnnotation")
library("MutationalPatterns")
library("data.table")
library("stringr")


# Data --------------------------------------------------------------------

#~~~ VCF files: final_consensus_snv_indel_icgc.controlled.tgz from https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel 
vcfDir <- file.path(dataDir,"icgc_vcf") ## only SNV files
vcfFiles <- list.files(vcfDir,recursive=T,full.names = T)
vcfNames <- list.files(vcfDir,recursive=T,full.names = F)

#~~~ pcawg sample sheet pcawg_sample_sheet.tsv from https://dcc.icgc.org/releases/PCAWG/donors_and_biospecimens
pcawg_sample_sheet <- read.delim(file = file.path(dataDir,"pcawg_sample_sheet.txt"))
pcawg_sample_sheet <- pcawg_sample_sheet[,c("aliquot_id","dcc_project_code","icgc_donor_id")] 

pcawg_sample_sheet[["aliquot_id"]] <- as.character(pcawg_sample_sheet[["aliquot_id"]])
pcawg_sample_sheet[["aliquot_id"]] <- paste(pcawg_sample_sheet[["aliquot_id"]], "vcf", sep = ".") 
d_pcawg_sample_sheet <- distinct(pcawg_sample_sheet)

#~~~ VCF files information 
vcfTbl <- as.data.frame(vcfNames) 
vcfTbl$vcfNames <- as.character(vcfTbl$vcfNames)
vcfTbl$project <- d_pcawg_sample_sheet[match(vcfTbl[["vcfNames"]], d_pcawg_sample_sheet[['aliquot_id']] ) , 'dcc_project_code'] 
vcfTbl$donor <- d_pcawg_sample_sheet[match(vcfTbl[["vcfNames"]], d_pcawg_sample_sheet[['aliquot_id']] ) , 'icgc_donor_id'] 

vcfTbl$project <- as.character(vcfTbl$project)
vcfTbl$donor <- as.character(vcfTbl$donor)



# Mutations --------------------------------------------------------------------

cancerTypes <- unique(vcfTbl$project)

for (i in 1:length(cancerTypes))
  
{
  project_id <- as.character(cancerTypes[[i]]) 
  saveName <- str_replace_all(project_id,"-","_")
  vcfTbl_i <- subset(vcfTbl, project == project_id)
  vcfNames_i <- vcfTbl_i[["vcfNames"]]
  vcfFiles_i <- file.path(vcfDir,vcfNames_i)

  final_CT <- c()
  final_GA <- c()
  
  for (j in 1:length(vcfFiles_i))
    
  {
    vcfs_j<- readVcf(vcfFiles_i[[j]], "hg19")
    mut_j <- rownames(vcfs_j)
    CT_j <- mut_j[grepl("C/T", mut_j)]
    GA_j <- mut_j[grepl("G/A", mut_j)]
    
    final_CT <- c(final_CT, CT_j)
    final_GA <- c(final_GA, GA_j)
    
  }
  
  
  ####CT table 
  chr_startMutation <- tstrsplit(final_CT, ":", names = c("chr", "beta")) 
  start_mutation <- tstrsplit(chr_startMutation$beta, "_", names = c("start", "mutation"))
  
  chr <- chr_startMutation$chr
  start <- start_mutation$start
  end <- start_mutation$start
  
  CT_final_project_table <- data.table(chr,start,end)
  write.table(CT_final_project_table, file = file.path(resDir,paste0("CtoT_",saveName,".tsv")), quote = FALSE, row.names=FALSE )
  
  
  ####GA table 
  chr_startMutation <- tstrsplit(final_GA, ":", names = c("chr", "beta")) 
  start_mutation <- tstrsplit(chr_startMutation$beta, "_", names = c("start", "mutation"))
  
  chr <- chr_startMutation$chr
  start <- start_mutation$start
  end <- start_mutation$start
  
  GA_final_project_table <- data.table(chr,start,end)
  write.table(GA_final_project_table, file = file.path(resDir,paste0("GtoA_",saveName,".tsv")), quote = FALSE, row.names=FALSE )
  
  
  

 } 





















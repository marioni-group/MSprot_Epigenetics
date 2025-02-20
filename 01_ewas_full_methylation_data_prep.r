###--------------------EWAS data prep for full GS data set-------------------------
#Methylation data: all GS – pre-regress age, sex, batch -> residuals --> scale

#Load required packages
library(tidyverse)
library(limma)
library(data.table)
library(QCEWAS)
library(here)

library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")

#Load required data
ext_data <- "/data_path/"
probes <- read.table(here(ext_data, "gs_methylation/cpgs_tokeep.txt"), header=F) 
protein_df <- readRDS(here("GS_ProteinGroups_RankTransformed.rds")) # rank transformed protein data
pheno <- readRDS(here(ext_data, "gs_methylation/GS20k/GS20k_Targets_18869.rds")) #phenotype and ids to keep

#subset to match IDs
pheno <- pheno %>%
filter(Sample_Name %in% protein_df$id) %>% #18869 -> 14671
select(Sample_Name, Sample_Sentrix_ID, age, sex, Batch) 

#Check samples also in DNAm_QC file
table(protein_df$id %in% DNAm_QC$Sample_Name) #N = 14671
table(protein_df$id %in% pheno$Sample_Name) #N = 14671
#phenotype object now contains the IDs passing QC and who have protein data

#set up locations for input & output 
datadir <- "/data_dir/input/"
methdir <- "methylation_data/"
localdir <- "/Local_dir/"

# Set up paralellizing
cores <- detectCores()  
cl <- makeCluster(6, outfile=paste0(datadir, "parallel_test.txt")) 
registerDoParallel(cl)

# Iterate per chromosome
foreach(i=1:22, .packages = "data.table") %dopar% { 
  print(paste0("Working on chromosome ",i))

  # Import 
  meth <- readRDS(paste0(localdir, "GS20k_chr", i, "_mvals.rds")) #Each chromosome saved out separately

  # Subset 
  meth <- meth[,which(colnames(meth) %in% pheno$Sample_Sentrix_ID)] # Subset to those with phenotype in question 
  meth <- meth[which(rownames(meth) %in% probes$V1),] # Subset to probes passing QC 

  # Match order of IDs in phenotype and methylation file 
  meth <- meth[,match(pheno$Sample_Sentrix_ID, colnames(meth))]

  # Mean impute - cannot have NAs in final file 
  meth <- apply(meth, 1, meanimpute)
  meth <- t(meth) 


  # Regression step - residualise for age, sex and batch 
  design.resid <- model.matrix(~age + as.factor(sex) +  as.factor(Batch), data=pheno)
  fit.resid <- limma::lmFit(meth, design.resid)
  gc()
  meth <- limma::residuals.MArrayLM(fit.resid, meth) 
  meth <- meth[!is.infinite(rowSums(meth)),] 
  rm(fit.resid)
  gc()

  # Write out CpGs & SSID
  cpgs <- as.data.frame(row.names(meth))
  SSID <- as.data.frame(colnames(meth))
  names(cpgs)[1] <- "CpG"
  names(SSID)[1] <- "SSID"
  fwrite(cpgs, paste0(datadir, methdir, "/GS_chr", i, "_cpgs.txt"),row.names=F)

  # Scale across each CpG (row) & transpose back 
  meth <- t(apply(meth,1,scale)) #this step loses SSID column names, keeps cpgs as another column
  colnames(meth) <- SSID$SSID  

  # Save out residualised file scaled
  meth <- data.table(meth, keep.rownames = TRUE)
  fwrite(meth, paste0(datadir, methdir, "/GS_chr", i, "_resid_mvals_scaled.txt"),row.names=F)  

 # Remove methylation object and clean up environment 
  rm(meth)
  gc()
}   

# End parallel
stopCluster(cl)

# Fuse methylation files - scaled 
files <- list.files(path = paste0(datadir, methdir, "/"), pattern = "_resid_mvals_scaled.txt") # Get files
files <- files[order(files)]
data <- rbindlist(lapply(paste0(datadir, methdir, "/", files),fread))
gc()

# Export fused methylation file - scaled 
data <- data.table(data, keep.rownames = TRUE)
data[1:5, 1:5]
data_t <- transpose(data, keep.names = "rn", make.names = "rn") 
data_t[1:5, 1:5]

#----------------------------------------------------------------------------------
#Check 
identical(data_t$rn, SSID$SSID) #TRUE

order <- data_t$rn
order <- as.data.frame(order)

#Save out order of methylation IDs
fwrite(x = order, paste0(datadir, methdir, "order_resid_scaled_methylation_data.csv"), sep = ",", row.names = F, col.names = TRUE, quote = F)  

#save out full methylation data in format which can be prepared for osca
write_rds(data_t, paste0(datadir, methdir, "GS_allchrom_resid_mvals_scaled.rds"), compress = 'none')

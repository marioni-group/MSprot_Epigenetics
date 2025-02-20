#-------------------EWAS data prep for full GS data set-------------------------
#Protein data: rank transformed, pre-regress age, sex, kinship matrix -> residuals,scale

##Load required packages
library(tidyverse)
library(data.table)
library(reshape)
library(finalfit)
library(forcats)
library(here)
library(Matrix)
library(coxme)

#--------------------------------------------------------------------------------
#load data
here()
ext_data <- "/data_dir/"
protein_df <- readRDS(here("data/GS_ProteinGroups_RankTransformed.rds")) #Protein data: rank transformed 
target <- readRDS(here(ext_data, "GS_methylation/GS20k_Targets_18869.rds")) # IDs for those passing QC 
kinship <- readRDS(here("ewas/input/kinship_matrix.rds")) #kinship matrix
output <- here("ewas/input/protein_data/") #folder location for data

#--------------------------------------------------------------------------
# Prepare data

protein_df$id <- as.character(protein_df$id) #change class for merging
table(protein_df$id %in% target$Sample_Name) #14671 individuals have protein data & methylation data passing QC

# Filter target file to contain IDs with protein data
target <- target %>%
  filter(target$Sample_Name %in% protein_df$id) #14671 remaining
protein_df <- protein_df %>%
  filter(protein_df$id %in% target$Sample_Name) #14671


# Match orders
names(target)[names(target) == "Sample_Name"] <- "id" #rename Sample_Name to id
target2 <- target[match(protein_df$id, target$id),]
identical(target2$id, protein_df$id) #check orders match
target <- target2 #reassign 


#--------------------------------------------------------------------------
#Linear model 

proteins <- colnames(protein_df) #extracts protein names 
proteins <- proteins[2:440]

residualised = protein_df
for (protein in proteins) {
  print(protein)
  residualised[[protein]] <- paste(resid(lmekin(residualised[[protein]] ~ target$age + as.factor(target$sex) + (1|residualised$id), 
  varlist = kinship*2,
   data = residualised, 
   na.action = na.exclude)))
} #Runs through each protein in data frame, 
#runs linear model, 
#takes the residuals and pastes them back into the data-frame

#write results out 
fwrite(residualised, file = here(output, "resid_proteins_basic_notscaled.csv", sep="\t", quote=F, row.names=F, col.names=T))

#________________________________________________________________________
##Scale results
proteins <- colnames(residualised) #extracts protein names 
proteins <- proteins[2:440]

#residualised data frame contains columns as type character -> need to be converted
residualised_2 <- residualised
residualised_2[2:440] <- lapply(residualised[2:440], as.numeric)

residualised_scaled_2 = residualised_2 
for (protein in proteins) {
    print(protein)
    residualised_scaled_2[[protein]] <- paste(scale(residualised_2[[protein]]))
}
residualised_scaled_2[1:5, 1:5]

#write out results
fwrite(residualised_scaled_2, file = here(output,"resid_proteins_basic_scaled.csv"),row.names=F)


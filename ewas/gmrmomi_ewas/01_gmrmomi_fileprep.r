#BayesR file preparation

library(data.table)
library(tidyverse)
library(missMethods)
library(tools)
library(bestNormalize)


# #set up locations for input & output
datadir <- "/input_dir/"


#FULL DATA - using methylation data prepped for EWAS________________________________________________________________________________________________________________________________
#This methylation data is pre-regressed for age, sex, batch, residuals taken and scaled (as below steps). 
# #1. Subset to IDs in phenotype prepped file (IDs in target file & with protein data)
# #2. Subset to CpGs that have passed QC
# #3. Regress out covariates of interest and take residuals (age, sex, batch in this case)
# #3. Scale the residuals
# #4. Ensure order of IDs matches that of the prepared phenotype data & save out order 
# #5. Save out in bin format. Go from format of rows(indiviuals) x columns(markers) -> vector where full data for CpG is printed consecutively
# #6. Also save out order of CpGs (will need this later)

meth <- readRDS("/data_dir/ewas/input_files/methylation_data/GS_allchrom_resid_mvals_scaled.rds")
meth[1:5, 1:5]
dim(meth)

#check the order of SSIDsdata
order <- fread("/data_dir/ewas/input_files/methylation_data/order_resid_scaled_methylation_data.csv") 
identical(order$order, meth$rn) #Check orders match
#Also save out SSID orders in bayesR folder
fwrite(order, file = "/data_dir/ewas/bayesr/order_files/SSID_order.csv") 

#write out cpg order
CpGs <- colnames(meth)
CpGs <- CpGs[-1]
fwrite(as.data.frame(CpGs), file = "/data_dir/ewas/bayesr/order_files/CpG_order_bayesR.csv")

#create .bin file
meth2 <- meth %>% select(-rn) #meth should be in the format of IDs(rows) x CpGs(columns)
data <- unlist(meth2, use.names = FALSE)
data <- as.vector(data) 
head(data) #data should start with the first column as a vector

new = file(paste0(datadir, "full_scaled_methylation_2602.bin"), "wb")
writeBin(data, new)
close(new)

#check bin file with below
help(readBin)
check = file(paste0(datadir, "full_scaled_methylation_2602.bin"), "rb")
readBin(check, numeric(), n = 6)
close(check)

#check with protein data (see below)


#FULL DATA________________________________________________________________________________________________________________________________
# Prepare phenotype files
# Aim is file with 3 columns: FID, IID, Phenotype, saved as space delimited
# NB - bayesR will not work if NA values present
# Remove column names and row names & save out (will need to ensure ID orders are matched with methylation data & covariates)
# Data should be scaled 
#Covariates data prep - to be available for residualising from the protein data
#Covariates to include: age, sex, logbmi, logpackyears to match model4 from osca model
#Cannot have NAs, so need to impute any missing data (KNN vs mean-impute)
#match to order

target <- readRDS("/ext_dir/GS_methylation/GS20k/GS20k_Targets_18869.rds")
cov <- readRDS("/data_dir/data/GS_phenos_internal_23Aug2023_REM.rds")
protein_df <- readRDS(file = "/data_dir/data/GS_ProteinGroups_RankTransformed3.rds")
order <- fread("/data_dir/ewas/input_files/methylation_data/order_resid_scaled_methylation_data.csv")

target_id <- target %>%
  select(Sample_Name, Sample_Sentrix_ID) 
target_id <- target_id %>% rename(id = Sample_Name)
covariates <- cov %>% select(id, age, sex, bmi, pack_years)
target_id$id <- as.character(target_id$id)
covariates$id <- as.character(covariates$id)


#add in SSID for matching
covariates <- left_join(covariates, target_id, by = "id")
head(covariates)

#match to ID order & subset
covariates <- covariates %>%
  filter(., Sample_Sentrix_ID %in% order$order)

covariates <- covariates[match(order$order, covariates$Sample_Sentrix_ID),]
identical(covariates$Sample_Sentrix_ID, order$order) #check orders match

#identify NA values
sum(is.na(covariates))
sapply(covariates, function(x) sum(is.na(x)))

#mean impute missing covariates (bmi n = 92, smoking n = 275)
covariates <- impute_mean(covariates, type = "columnwise", convert_tibble = TRUE)
sapply(covariates, function(x) sum(is.na(x)))

#create log covariates
covariates <- covariates %>%
  mutate(logbmi = log10(bmi),
            logpckyrs = log10(pack_years + 1))

#Use one protein for phenotype file
protein_df <- readRDS(file = "/data_dir/data/GS_ProteinGroups_RankTransformed.rds")
pheno <- protein_df %>% 
  select(!matches("\\.")) #filter to only single proteins
sum(is.na(pheno)) #Check for NA values -> need to impute if present

#add SSID and subset to those passed QC for individuals
target <- readRDS("/ext_dir/GS_methylation/GS20k/GS20k_Targets_18869.rds") ### ID for those passing QC
SSIDs <- target %>% select(Sample_Name, Sample_Sentrix_ID)
rownames(SSIDs) <- NULL
SSIDs <- SSIDs %>% rename(id = Sample_Name)
pheno$id <- as.character(pheno$id)


pheno <- left_join(pheno, SSIDs, by = "id") #add SSID for matching to methylation data
pheno <- pheno %>%
  filter(Sample_Sentrix_ID %in% target$Sample_Sentrix_ID) #filter to individuals in target file
sum(is.na(pheno)) #second check for NA after manipulation
pheno <- pheno %>% dplyr::select(c(Sample_Sentrix_ID, everything()))
pheno <- pheno %>% dplyr::select(-id)

#match order to covariates order
pheno <- pheno[match(order$order, pheno$Sample_Sentrix_ID),]
identical(pheno$Sample_Sentrix_ID, order$order)
identical(pheno$Sample_Sentrix_ID, covariates$Sample_Sentrix_ID)

#Residualise
#create list of protein names
proteins <- unlist(colnames(pheno)[-1])
proteins <- as.list(proteins)

residualised <- pheno %>% select(-Sample_Sentrix_ID)
for (protein in proteins) {
  print(protein) 
 model <- lm(residualised[[protein]] ~ covariates$age + as.factor(covariates$sex) + covariates$logbmi + covariates$logpckyrs,
   na.action = na.exclude)
  residualised[[protein]] <- resid(model) 
}

residualised[1:5, 1:5]
pheno[1:5, 1:5]

#Scale 
scaled_pheno <- residualised #create duplicate df for scaling

for (protein in proteins) {
    print(protein)
    scaled_pheno[[protein]] <- paste(scale(residualised[[protein]]))
}

#check
scaled_pheno[1:5, 1:5]
residualised[1:5, 1:5]

#Create two ID columns as required for bayesR input
scaled_pheno <- scaled_pheno %>% 
  mutate(FID = pheno$Sample_Sentrix_ID, IID = pheno$Sample_Sentrix_ID)

scaled_pheno <- scaled_pheno %>% 
  dplyr::select(c(FID, IID, everything()))
scaled_pheno[1:5, 1:5]

ID_order <- order
scaled_pheno <- scaled_pheno[match(ID_order$order, scaled_pheno$V1),]
identical(scaled_pheno$FID, order$order) #TRUE

#Save out, space or tab delimited, no column names
fwrite(scaled_pheno, file = "/data_dir/ewas/bayesr/scaled_proteins.phen", sep = " ", col.names = FALSE)
#this file also has the ID order

#-----------------------------------------------------------------------------------
#Write out protein data in separate files for input to bayesr

location <- "/home_dir/" #set location for output

for(i in 3:135) { 
name <- as.character(names(scaled_pheno)[i])
name_data <- scaled_pheno %>%
  select(all_of(c("FID", "IID", name)))
fwrite(name_data, file = paste0(location, name,".phen"), sep=" ", col.names = FALSE)
}

#________________________________________________________________________________________________________________________________
###Dim file
#x and y dimensions of original methylation matrix
#save with no column names, space separated

#example
example_dim <- fread("/data_dir/ewas/bayesr/gmrm-omi/example/example.dim")
head(example_dim)

test_dim <- example_dim %>% mutate(V1 = 14671, V2 = 752722)
fwrite(test_dim, "/home_dir/gmrm/full_dim.dim", sep = " ", col.names = FALSE)

#________________________________________________________________________________________________________________________________
###Preparation of groups file
# 2 columns, first are probes, second is group assignation (cpgs vs. SNPs vs....)
# column length should = length of no. probes

example_group <- fread("/data_dir/ewas/bayesr/gmrm-omi/example/example.gri")

#take CpG labels from prepped methylation data
Cpgs <- colnames(meth)
Cpgs <- Cpgs[-1]
Cpgs <- as.data.frame(Cpgs)
Cpgs$V2 <- c(0)


fwrite(Cpgs, "/home_dir/gmrm/full_gri.gri", sep = " ", col.names = FALSE)

#________________________________________________________________________________________________________________________________
###Preparation of group mixture file

#Look at example
example_grm <- fread("/data_dir/ewas/bayesr/gmrm-omi/example/example.grm")
head(example_grm)

#edit this to save out for test use
grm <- example_grm %>% select(V1,V4,V5)
grm <- grm %>% mutate(V6 = 0.1)
head(grm)
grm <- grm %>% rename(V2 = V4, V3 = V5, V4 = V6)
head(grm)
grm <- as.matrix(grm)
names(grm) <- NULL
#Save out
fwrite(grm, file = "/data_dir/ewas/bayesr/test_grm.grm", sep="\t", col.names = FALSE)




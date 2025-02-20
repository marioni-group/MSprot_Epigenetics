#Step 1 for Episcores:: Data prep

#Load required packages
library(tidyverse)
library(limma)
library(data.table)
library(QCEWAS)
library(missMethods)

library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")

#Split into test and train groups
pheno <- readRDS("/data_dir/data/GS_phenos_internal_23Aug2023_REM.rds")
target <- readRDS("/ext_dir/GS_methylation/GS20k/GS20k_Targets_18869.rds")
pedigree <- fread("/data_dir/data/pedigree.txt")
protein_df <- readRDS("/data_dir/data/GS_ProteinGroups_RankTransformed_23Aug2023.rds")
probes <- read.table("/ext_dir/gs_osca/data/cpgs_tokeep.txt", header=F) #n = 752722

#Pre-prep: create a file to help with subsetting to test set
#1. Add family ID to target file
#2. Create training and testing subsets
#3. Subset to IDs in protein data 
#4. Remove IDs from training set that have a family ID matching wv3 family ID 


#1. Add family ID and wave to target file
pedigree <- pedigree %>% dplyr::select(famid, volid) %>%
  dplyr::rename(id = volid)
pedigree$id <- as.character(pedigree$id)  

#Create training_target set
training_target <- target %>% 
  filter(., !grepl("W3", Batch)) #remove wave 3
training_target <- training_target %>% dplyr::rename(id = Sample_Name)
training_target <- left_join(training_target, pedigree, by = "id")
training_target <- training_target %>%
  filter(., id %in% protein_df$id) #11208

test_target <- target %>%
  filter(., grepl("W3", Batch)) #select wave 3
test_target <- test_target %>% dplyr::rename(id = Sample_Name)
test_target <- left_join(test_target, pedigree, by = "id")
test_target <- test_target %>%
  filter(., id %in% protein_df$id) #3463

training_target <- training_target %>% #remove anyone in training set with same famid as anyone in w3
  filter(., !famid %in% test_target$famid) #6816

#save out separate target files for later use
fwrite(training_target, file = "/data_dir/ewas/bayesr/episcore/training_target.csv")
fwrite(test_target, file = "/data_dir/ewas/bayesr/episcore/test_target.csv")  

#Prepare methylation Data
# Set up input & output locations
datadir <- "/data_dir/"
methdir <- "methylation_data/"
localdir <- "/home_dir/GS_20k/"

#___________________________________________________________________________________________________________
#1. Load Data by Chromosome
#2. Subset to training set
#3. Subset to those with protein data
#4. Subset probes to those passing QC
#5. Order
#6. Impute if any NA values
#7. Regress for age, sex, batch
#8. Write out CpGs
#9. Scale across each CpG 
#10. Save out & exit 
#11. Fuse methylation files
#12. Make sure order is saved
#13. ?Save out full object

# Set up paralellizing
cores <- detectCores()  
cl <- makeCluster(3, outfile=paste0(datadir, "parallel_test.txt")) 
registerDoParallel(cl)

# Iterate per chromosome
foreach(i=1:22, .packages = "data.table") %dopar% { 

print(paste0("Working on chromosome ",i))

# Import 
  meth <- readRDS(paste0(localdir, "GS20k_chr", i, "_mvals.rds")) #Each chromosome saved out separately

  # Subset 
  meth <- meth[,which(colnames(meth) %in% training_target$Sample_Sentrix_ID)] # Subset to those with phenotype in question 
  meth <- meth[which(rownames(meth) %in% probes$V1),] # Subset to probes passing QC     

  # Match order of IDs in phenotype and methylation file 
  meth <- meth[,match(training_target$Sample_Sentrix_ID, colnames(meth))]

   # Regression step - residualise for age, sex and batch 
  design.resid <- model.matrix(~age + as.factor(sex) +  as.factor(Batch), data=training_target)
  fit.resid <- limma::lmFit(meth, design.resid)
  gc()
  meth <- limma::residuals.MArrayLM(fit.resid, meth) 
  meth <- meth[!is.infinite(rowSums(meth)),]
  rm(fit.resid)
  gc()

 # Write out CpGs 
  cpgs <- as.data.frame(row.names(meth))
  SSID <- as.data.frame(colnames(meth))
  names(cpgs)[1] <- "CpG"
  names(SSID)[1] <- "SSID"
  fwrite(cpgs, paste0(datadir, methdir, "/GS_chr", i, "_cpgs.txt"),row.names=F)

 # Scale each CpG (column on transposed df) & transpose back    
  meth <- scale(t(meth)) 
  meth <- t(meth)
  gc()

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

dim(data)

# Export fused methylation file - scaled 
data <- data.table(data, keep.rownames = TRUE)
data[1:5, 1:5]
data_t <- transpose(data, keep.names = "rn", make.names = "rn") 
data_t[1:5, 1:5]

#Check 
all.equal(data_t$rn, training_target$Sample_Sentrix_ID) #TRUE
identical(data_t$rn, training_target$Sample_Sentrix_ID) #TRUE

order <- data_t$rn
order <- as.data.frame(order)

cpgs <- colnames(data_t)
cpgs <- cpgs[-1]
cpgs <- as.data.frame(cpgs)

fwrite(x = cpgs, paste0(datadir, methdir, "cpg_order_training.csv"), sep = ",", row.names = F, col.names = TRUE, quote = F)

#Save out order of methylation IDs
fwrite(x = order, paste0(datadir, methdir, "training_order_resid_scaled_methylation_data.csv"), sep = ",", row.names = F, col.names = TRUE, quote = F)

#create .bin file
meth2 <- data_t %>% select(-rn) #meth should be in the format of IDs(rows) x CpGs(columns)
data <- unlist(meth2, use.names = FALSE)
data <- as.vector(data) 
head(data) #data should start with the first column as a vector

new = file("/home_dir/gmrm_scores/input/scores_methylation.bin", "wb")
writeBin(data, new)
close(new)


#check bin file with below
help(readBin)
check = file("/home_dir/gmrm_scores/input/scores_methylation.bin", "rb")
readBin(check, numeric(), n = 6)
close(check)


#___________________________________________________________________________________________________________
#Phenotype data prep
#1. Subset protein data to IDs in training_target
#2. Match order to methylation data order
#3. subset phenotype data & impute any NA values if present - check covariates
#4. regress for age, sex, logbmi, logpackyears
#5. take residuals
#6. scale
#7. check order
#8. Save out as 3 columns: FID, IID, Phenotype, saved as space delimited

training_target <- fread("/data_dir/ewas/bayesr/episcore/training_target.csv")

#1. & 2: subset & match order
protein_df <- protein_df %>% dplyr::filter(., id %in% training_target$id)
protein_df <- protein_df %>% 
  select(!matches("\\.")) #remove protein groups
training_target$id <- as.character(training_target$id)
protein_df$id <- as.character(protein_df$id)
identical(training_target$id, protein_df$id) #FALSE
protein_df <- protein_df[match(training_target$id, protein_df$id),]
identical(training_target$id, protein_df$id) #TRUE

#Check order with ID order from methylation file
meth_order <- fread(paste0(datadir, methdir, "training_order_resid_scaled_methylation_data.csv"))
identical(training_target$Sample_Sentrix_ID, meth_order$order) #TRUE


#3. 
covariates <- pheno %>% dplyr::filter(., id %in% training_target$id)
identical(training_target$id, covariates$id) #FALSE
covariates <- covariates[match(training_target$id, covariates$id),]
covariates$id <- as.character(covariates$id)
identical(training_target$id, covariates$id) #TRUE
covariates <- covariates %>% dplyr::select(id, age, sex, bmi, pack_years)

#identify NA values
sum(is.na(covariates))
sapply(covariates, function(x) sum(is.na(x))) #44 NA bmi, 155 NA pack_years

#mean impute
covariates <- impute_mean(covariates, type = "columnwise", convert_tibble = TRUE)
sapply(covariates, function(x) sum(is.na(x))) #0

#4. create log covariates
covariates <- covariates %>%
  mutate(logbmi = log10(bmi),
            logpckyrs = log10(pack_years + 1))

# check order again/confirm have made correct subset
identical(training_target$id, covariates$id) #TRUE
identical(training_target$Sample_Sentrix_ID, meth_order$order) #TRUE
identical(covariates$id, protein_df$id) #TRUE


##Residualise
proteins <- unlist(colnames(protein_df)[-1])
proteins <- as.list(proteins)

residualised <- protein_df %>% select(-id)
for (protein in proteins) {
  print(protein) 
 model <- lm(residualised[[protein]] ~ covariates$age + as.factor(covariates$sex) + covariates$logbmi + covariates$logpckyrs,
   na.action = na.exclude)
  residualised[[protein]] <- resid(model) 
}

residualised[1:5, 1:5]
protein_df[1:5, 1:5]

#Scale 
scaled_pheno <- scale(residualised) #create duplicate df for scaling

#compare dataframes
scaled_pheno[1:5, 1:5]
residualised[1:5, 1:5]

#Create two ID columns as required for bayesR input
identical(training_target$id, protein_df$id)
scaled_pheno <- as.data.frame(scaled_pheno) %>% 
  mutate(FID = training_target$Sample_Sentrix_ID, IID = training_target$Sample_Sentrix_ID)

scaled_pheno <- scaled_pheno %>% 
  dplyr::select(c(FID, IID, everything()))
scaled_pheno[1:5, 1:5]
identical(scaled_pheno$FID, training_target$Sample_Sentrix_ID) #TRUE
identical(scaled_pheno$FID, meth_order$order) #TRUE

#Save out, space or tab delimited, no column names
fwrite(scaled_pheno, file = "/data_dir/ewas/bayesr/episcore/protein_data/training_scaled_proteins.phen", sep = " ", col.names = TRUE)
#this file also has the ID order

#-----------------------------------------------------------------------------------
#Write out protein data in separate files for input to bayesr

location <- "/home_dir/gmrm_scores/input/" #set location for output

for(i in 3:135) { 
name <- as.character(names(scaled_pheno)[i])
name_data <- scaled_pheno %>%
  select(all_of(c("FID", "IID", name)))
fwrite(name_data, file = paste0(location, name,".phen"), sep=" ", col.names = FALSE)
}


###Dim file________________________________________________________________________________________________________________________________
#x and y dimensions of original methylation matrix
#save with no column names, space separated

#example
example_dim <- fread("/data_dir/ewas/bayesr/gmrm-omi/example/example.dim")
head(example_dim)
test_dim <- example_dim %>% mutate(V1 = 6816, V2 = 752722)
fwrite(test_dim, "/home_dir/input/scores_dim.dim", sep = " ", col.names = FALSE)


#### Preparation of groups file________________________________________________________________________________________________________________________________
# 2 columns, first are probes, second is group assignation (cpgs vs. SNPs vs....)
# column length should = length of no. probes

example_group <- fread("/data_dir/ewas/bayesr/gmrm-omi/example/example.gri")

#take CpG labels from prepped methylation data ****check this
probes <- fread(paste0(datadir, methdir, "cpg_order_training.csv"))
probes$V2 <- c(0)
fwrite(probes, "/home_dir/gmrm_scores/input/full_gri.gri", sep = " ", col.names = FALSE)


#grm file________________________________________________________________________________________________________________________________
grm <- fread("/data_dir/ewas/bayesr/input/test_grm.grm") #use same prior mix as full ewas
head(grm)
fwrite(grm, "/home_dir/gmrm_scores/input/grm.grm", sep = " ", col.names = FALSE)



##TEST DATA##______________________________________________________________________________________________________________________________

#Prepare methylation data for episcore test_____________________
test_target <- fread("/data_dir/ewas/bayesr/episcore/test_target.csv")
probes <- read.table("/ext_dir/gs_osca/data/cpgs_tokeep.txt", header=F) 

head(test_target)
dim(test_target)

datadir <- "/data_dir/"
methdir <- "/test_data/"
localdir <- "/home_dir/GS_20k/"

#___________________________________________________________________________________________________________
#1. Load Data by Chromosome
#2. Subset to test set
#3. Subset to those with protein data
#4. Subset probes to those passing QC
#5. Match Order
#6. Impute if any NA values
#7. Write out CpGs
#8. Scale across each CpG 
#9. Save out & exit 
#10. Fuse methylation files
#11. Make sure order is saved
#12. ?Save out full object

# Set up paralellizing
cores <- detectCores()  
cl <- makeCluster(3, outfile=paste0(datadir, "parallel_test.txt")) 
registerDoParallel(cl)

# Iterate per chromosome
foreach(i=1:22, .packages = "data.table") %dopar% { 

print(paste0("Working on chromosome ",i))

# Import 
  meth <- readRDS(paste0(localdir, "GS20k_chr", i, "_mvals.rds")) #Each chromosome saved out separately

  # Subset 
  meth <- meth[,which(colnames(meth) %in% test_target$Sample_Sentrix_ID)] # Subset to those with phenotype in question 
  meth <- meth[which(rownames(meth) %in% probes$V1),] # Subset to probes passing QC     

  # Match order of IDs in phenotype and methylation file 
  meth <- meth[,match(test_target$Sample_Sentrix_ID, colnames(meth))]

 # Write out CpGs 
  cpgs <- as.data.frame(row.names(meth))
  SSID <- as.data.frame(colnames(meth))
  names(cpgs)[1] <- "CpG"
  names(SSID)[1] <- "SSID"
  fwrite(cpgs, paste0(datadir, methdir, "/GS_chr", i, "_cpgs.txt"),row.names=F)

 # Scale each CpG (column on transposed df) & transpose back    
  meth <- scale(t(meth)) 
  meth <- t(meth)
  gc()

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

dim(data)

# Export fused methylation file - scaled 
data1 <- data.table(data, keep.rownames = TRUE)
data1[1:5, 1:5]
codata_t <- transpose(data1, keep.names = "rn", make.names = "rn") 
codata_t[1:5, 1:5]

#Check 
all.equal(codata_t$rn, test_target$Sample_Sentrix_ID) #TRUE
identical(codata_t$rn, test_target$Sample_Sentrix_ID) #TRUE

order <- data_t$rn
order <- as.data.frame(order)

cpgs <- colnames(data_t)
cpgs <- cpgs[-1]
cpgs <- as.data.frame(cpgs)

fwrite(x = cpgs, paste0(datadir, methdir, "/cpg_order_test.csv"), sep = ",", row.names = F, col.names = TRUE, quote = F)

#Save out order of methylation IDs
fwrite(x = order, paste0(datadir, methdir, "/test_order_scaled_methylation_data.csv"), sep = ",", row.names = F, col.names = TRUE, quote = F)
fwrite(x = codata_t, paste0(datadir, methdir, "/test_scaled_meth_data.csv"), sep = ",", row.names = F, col.names = TRUE, quote = F)


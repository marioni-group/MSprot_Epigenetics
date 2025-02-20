#OSCA fast-linear covariate files for use in models 2 - 4

library(tidyverse)
library(data.table)
library(here)
library(performance)
library(missMethods)

 
#-------------------------------------------------------------------------------------------------------------------------------------------------
#Model 2: add estimated white cells as covariates
here()
data_dir <- "/data_dir/" #for group data external to project folder
target <- readRDS(paste0(data_dir, "GS/GS_methylation/GS20k/GS20k_Targets_18869.rds")) #Covariates: white cell counts, age and sex
str(target) #7 white cell types

protein_df <- readRDS("data/GS_ProteinGroups_RankTransformed_23Aug2023.rds") #need to ensure IDs are correct
order <- read.csv("ewas/input/methylation_data/order_resid_scaled_methylation_data.csv") #order of IDs for methylation input: created in script 01

#filter target to covariates of interest
target <- target %>%
select(Sample_Name, Sample_Sentrix_ID, sex, age, CD8T, Mono, Bcell, Eos, NK, Neu, CD4T)
str(target)


#Filter by Ids in protein_df 
target <- target %>%
  filter(Sample_Name %in% protein_df$id) 
identical(order$order, target$Sample_Sentrix_ID)  #TRUE

#will need to drop one wc proportion as all add to 1 -> collinearity
#inspect correlation matrix to see which likely candidate -> neutrophils
cordf <- cor(target[4:11])

#Run test lm to identify which to remove
#Repeat using protein data: "protein"
target$id <- target$Sample_Name
protein_example <- protein_df %>% select(id,protein)
protein_example$id <- as.character(protein_example$id)
pheno <- left_join(target, protein_example, by = "id")

test.model<-lm(protein ~ age + sex + Bcell + CD4T + CD8T + Eos + Mono + NK + Neu, data=pheno)
summary(test.model)
confint(test.model)
anova(test.model)

#Test
check <- test.model %>% 
  check_model() #High collinearity with interaction terms

check_heteroscedasticity(test.model) #variance 
check_collinearity(test.model) #check collinearity
check_normality(test.model)
model_performance(test.model)

#Remove neutrophils are re-run
test.model<-lm(protein ~ age + sex + Bcell + CD4T + CD8T + Eos + Mono + NK, data=pheno)
summary(test.model)
confint(test.model)
anova(test.model)
test.model %>% 
  check_model()
check_heteroscedasticity(test.model) #variance 
check_collinearity(test.model) #check collinearity
check_normality(test.model)
model_performance(test.model)
#removing neutrophils has resolved model issues

# Make quantitative covariates file (cellcounts only) 
# First two columns are FID and IID as before
# Drop neutrophils as this resolved collinearity issues
#check orders again
identical(target$Sample_Sentrix_ID, order$order) #TRUE

quant_cov <- data.frame(FID = target$Sample_Sentrix_ID,
	                    IID = target$Sample_Sentrix_ID,
	                    CD8T = target$CD8T,
                        Mono = target$Mono,
                        Bcell = target$Bcell,
                        Eos = target$Eos,
                        NK = target$NK,
                        CD4T = target$CD4T )
sum(is.na(quant_cov)) #0                    
write.table(quant_cov, file="ewas/input/protein/batches/model/quant.cov", quote =F, row.names=F, sep=' ')

#-------------------------------------------------------------------------------------------------------------------------------------------------
#Model 3: wcc + bmi + smoking

here()
target <- readRDS(here(data_dir, "GS/GS_methylation/GS20k/GS20k_Targets_18869.rds"))
pheno <- readRDS("data/GS_phenos_internal_23Aug2023_REM.rds") #phenotypes
protein_df <- readRDS("data/GS_ProteinGroups_RankTransformed.rds") #protein data
order <- read.csv(here(datadir, "/order_resid_scaled_methylation_data.csv")) #written in 01_ewas_full_methylation_dataprep.r

#select id and SSID from target file
target <- target %>%
  select(Sample_Name, Sample_Sentrix_ID) %>%
  mutate(id = Sample_Name) 
head(target)

#join to phenotype file
pheno$id <- as.character(pheno$id)
pheno <- left_join(pheno, target, by = "id")

#filter to SSID in order file
pheno <- pheno %>%
  dplyr::filter(Sample_Sentrix_ID %in% order$order) #14671

#match order
pheno <- pheno[match(order$order, pheno$Sample_Sentrix_ID),]
identical(pheno$Sample_Sentrix_ID, order$order) #TRUE

#select covariates of interest: age, sex, bmi, smoking
pheno <- pheno %>%
  dplyr::select(id, Sample_Sentrix_ID, age, sex, bmi, pack_years)

#check for missing values
sapply(pheno, function(x) sum(is.na(x))) #both bmi and pack_years have missing values


#mean impute missing values
pheno <- impute_mean(pheno, type = "columnwise", convert_tibble = TRUE)
sapply(pheno, function(x) sum(is.na(x))) #no remaining missing values

pheno %>%
  dplyr::select(bmi) %>%
  ggplot(aes(x = bmi)) +
  geom_density() #right skewed distribution

 pheno %>%
  dplyr::select(pack_years) %>%
  ggplot(aes(x = log10(pack_years + 1))) +
  geom_density()  #some improvement 

#create log bmi and pack years columns
pheno <- pheno %>%
  mutate(logbmi = log10(bmi)) %>%
  mutate(logpackyears = log10(pack_years + 1))

#add to previous quant_cov file to create covariate file for next model
quant_cov <- fread("ewas/input/protein_data/batches/model/quant.cov")

#check orders
pheno <- pheno %>%
  select(Sample_Sentrix_ID, logbmi, logpackyears)
identical(pheno$Sample_Sentrix_ID, quant_cov$FID) #TRUE
identical(pheno$Sample_Sentrix_ID, order$order) #TRUE
pheno$FID <- pheno$Sample_Sentrix_ID

quant_smok <- left_join(pheno, quant_cov, by = "FID")
quant_smok <- quant_smok %>%
  dplyr::select(-c(Sample_Sentrix_ID))

#maintain order as per quant_cov
quant_smok <- quant_smok %>% 
  dplyr::select(FID, IID,CD8T, Mono, Bcell, Eos, NK, CD4T, logbmi, logpackyears)
identical(quant_smok$FID, order$order) #TRUE

#save out 
here()
write.table(quant_smok, file=here(data_dir, "model/quant_smok.cov", quote =F, row.names=F, sep=' '))

#-------------------------------------------------------------------------------------------------------------------------------------------------
#Model 4: wcc + bmi + smoking + 20 Methylation PCs

#load data
methPC <- fread(here(data_dir, "GS_methylation/GS20k/2023-01-04_gs20k_pcs.csv"))
cov <- fread(here(data_dir, "model/quant_smok.cov"))
order <- read.csv(here(data_dir, "/order_resid_scaled_methylation_data.csv"))

#Prep for osca
methPC$IID <- methPC$ID
cov4 <- left_join(cov, methPC, by = "IID")
cov4 <- cov4[match(order$order, cov4$IID),]
identical(cov4$IID, order$order) #TRUE

#Select 20 first PCs
cov_PCs <- cov4[,1:30]
colnames(cov_PCs)
write.table(cov_PCs, file=here(datadir, "/model/quant_20_PC.cov", quote =F, row.names=F, sep=' '))


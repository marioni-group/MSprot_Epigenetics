#Episcores and CVD

#1. Comparing EpiScores and Proteins as predictors of cardiovascular disease
# in testing subset
# load predicted protein levels for all test individuals
# load actual protein levels for all test individuals
# calculate numbers for CVD phenotypes
# set up cox model data
# run cox models for prediction
# compare to actual protein levels as predictors

#load libraries
library(tidyverse)
library(data.table)
library(stringr)
library(here)
library(magrittr)
library(lubridate)

#Cox model libraries
library(survival)
library(survminer)

#for plot ordering
library(forcats)

#_____________________________________________________________________________________
#### DATA PREP ###

#Load files
here()
data_path <- here::here("") #set data path
test_target <- fread(here::here(data_path, "test_target.csv"))  #test group target file
episcore_df <- fread(here(data_path, "episcore_protein_predictions.csv")) #load predicted protein levels
protein_df <- readRDS(here::here("data/GS_ProteinGroups_RankTransformed.rds")) #true protein levels
protein_list <- episcore_df$protein #extract list of 133 proteins for investigation
diseases <- fread(here("data/disease_codes_combined.csv")) #Disease phenotypes
deaths1 <- readRDS(here(data_path, "proteomics_mortality_phenodata.rds")) #CVD deaths
gsappt <- fread(here("data/GS_appt.txt")) #appointment dates

#subset protein_df to 133 unique proteins and ids for test group
protein_df <- protein_df %>% 
  dplyr::select(id, protein_list) %>%
  dplyr::filter(id %in% test_target$id)

#scale protein data
scaled_pheno <- protein_df
scaled_pheno[,2:134] <- scale(scaled_pheno[,2:134])  

#add in Sample_Sentrix_ID for ease of further wrangling
ids <- test_target %>% dplyr::select(id, Sample_Sentrix_ID, age, sex)
scaled_pheno <- left_join(scaled_pheno, ids, by = "id")

#________________________________________________________________________________________
# Generate outcome df for individuals in test set: either cvd diagnosis or cvd death prior to censor date 202308
# Individuals who die of other causes prior to 202308 should also be censored
# Individuals who have a diagnosis of cvd prior to study recruitment should be excluded
# CVD event = composite cvd outcome (for now) = cvd-related diagnosis (stroke, CHD-NOS and myocardial infarction) OR cvd-related death 
#____________Details: - Secondary care data only
#___________________: - Same conditions used as Ola


#subset disease phenotypes to ids in test target
diseases <- diseases %>% dplyr::filter(.,id %in% test_target$id)
diseases <- diseases %>% left_join(ids, by = "id") #add Sample_IDs

#set-up composite cvd outcome
disease_list <- unique(diseases$Disease)
cvd_list <- c("CHD_NOS", "Isch_stroke", "myocardial_infarction")
#identical(test_target$id, diseases$id)

#create subset dataframe of only cvd
cvd_only <- diseases %>% 
  dplyr::filter(Disease %in% cvd_list & Source == "Secondary_Care" & incident == "1")
dim(cvd_only) #263 cases incident diagnosis of cvd

#Look at prevalent disease
prevalent_cvd <- diseases %>%
  dplyr::filter(Disease %in% cvd_list & Source == "Secondary_Care" & incident == "0" ) 

#check for any missing values
sapply(cvd_only, function(x) sum(is.na(x))) #0

#check diagnosis dates are after appointment date
cvd_only %>% dplyr::filter(dt1_ym < gs_appt) #0 entries
cvd_only %>% dplyr::filter(dt1_ym > gs_appt) %>% dim() 
length(unique(cvd_only$id)) 

#retain earliest diagnosis date where there are multiple entries
unique(cvd_only$Disease) #check all are in list of cvd diseases
cvd_only_filt <- cvd_only %>%
 dplyr::group_by(id) %>% #for each id, filter to earliest diagnosis date
 filter(dt1_ym == min(dt1_ym))  
length(unique(cvd_only_filt$id)) 
cvd_only_filt %>% group_by(id) %>% dplyr::filter(n()>1) %>% arrange(id) 
cvd_only_filt <- cvd_only_filt %>%
  dplyr::select(-Disease) %>%
  distinct()  

cvd_only_filt <- cvd_only_filt %>% mutate(cvd = c("1"))  #create new column to indicate presence of cvd
#check for any missing values
sapply(cvd_only_filt, function(x) sum(is.na(x))) #0

#generate cvd deaths info & check if any overlap with cvd diagnoses
heart_deaths <- c("I00", "I01", "I02", "I03", "I04", "I05", "I06", "I07", "I08", "I09", "I11",
         "I13", "I20", "I21", "I22", "I23", "I24", "I25", "I26", "I27", "I28", "I29",
         "I30", "I31", "I32", "I33", "I34", "I35", "I36", "I37", "I38", "I39", "I40",
         "I41", "I42", "I43", "I44", "I45", "I46", "I47", "I48", "I49", "I50", "I51")
hypertension_deaths <- c("I10", "I12", "I15")
cerebrovascular_deaths <- c("I60", "I61", "I62", "I63", "I64", "I65", "I66", "I67", "I68", "I69")
#add in remaining codes to encompass all I00-I99
complete_list <- c("I14", "I16", "I17", "I18", "I19", "I52", "I53", "I54", "I55", "I56", "I57", "I58", "I59",
                    "I70", "I71", "I72", "I73", "I74", "I75", "I76", "I77", "I78", "I79",
                    "I80", "I81", "I82", "I83", "I84", "I85", "I86", "I87", "I88", "I89",
                    "I90", "I91", "I92", "I93", "I94", "I95", "I96", "I97", "I98", "I99")
all <- c(heart_deaths, hypertension_deaths, cerebrovascular_deaths, complete_list)  #length = 100                

#create cvd and alternative info
deaths1$id <- as.character(deaths1$id)
deaths_filt <- deaths1 %>% 
  filter(id %in% test_target$id) #subset to test IDs
  
deaths_filt <- deaths_filt %>% #create new variable for cardiac_death - checks for ICD code in cause of death columns
  mutate(cardiac_death = ifelse(
          if_any(everything(), ~. %in% c(heart_deaths, hypertension_deaths, cerebrovascular_deaths, complete_list)),
                    "1", "0")) %>% 
    dplyr::select(id,dod_ym,dead, cardiac_death)
  
#check numbers
deaths_filt %>% dplyr::filter(dead == "1") %>% dim() 
deaths_filt %>% dplyr::filter(dead == "1" & cardiac_death == "1") %>% dim() 
deaths_filt %>% dplyr::filter(dead == "1" & cardiac_death == "0") %>% dim() 

#compare to original deaths file
deaths1 %>% dplyr::filter(id %in% test_target$id & if_any(everything(), ~. %in% c(heart_deaths, hypertension_deaths, cerebrovascular_deaths))) %>% dim() 

#add in appointment date and/or diagnosis date for all test set individuals

#first add gsappt to target df
sapply(ids, function(x) sum(is.na(x))) #0 missing values
sapply(gsappt, function(x) sum(is.na(x))) #0 missing values
ids <- ids %>% left_join(gsappt, by = "id")
sapply(ids, function(x) sum(is.na(x))) #

ids <- ids %>% 
  dplyr::select(id, Sample_Sentrix_ID, age, sex, appt) %>%
  drop_na() #remove those with missing appt data

cvd_merge <- cvd_only_filt %>% dplyr::select(id,dt1_ym,cvd) #prep cvd df for merging
sapply(cvd_merge, function(x) sum(is.na(x)))  #0 missing 


#add cleaned cvd data to id df
df_cox <- left_join(ids, cvd_merge, by = "id")
sapply(df_cox, function(x) sum(is.na(x)))   

#add deaths data
deaths_filt$id <- as.character(deaths_filt$id)
df_cox$id <- as.character(df_cox$id)
df_cox <- df_cox %>% left_join(deaths_filt, by = "id")
sapply(df_cox, function(x) sum(is.na(x)))
 
#check for overlap between cvd cases and deaths & any data entry anomalies
df_cox %>% filter(cvd == "1" & dead == "1") %>% dim() 
df_cox %>% filter(cvd == "1" & dead == "1" & cardiac_death == "1") %>% dim() 
df_cox %>% filter(cvd == "1" & dead == "1" & dt1_ym > dod_ym) 
df_cox %>% filter(cvd == "1" & dead == "1" & cardiac_death == "1" & dt1_ym > dod_ym) 

#remove any individuals with prevalent cvd
head(prevalent_cvd)
df_cox <- df_cox %>% filter(!id %in% prevalent_cvd$id) 
df_cox %>% dplyr::filter(cvd == "1") %>% dim() 
table(df_cox$cvd) 
table(df_cox$cardiac_death) 
table(df_cox$dead) 

#set up variables for cox analysis
#create composite outcome, where 1 signifies either cardiac disease or cardiac death and 0 signifies neither
df_cox <- df_cox %>% mutate(composite_cardiac_outcome = 
ifelse(cvd == "1" | cardiac_death == "1", "1", "0"))


#add in censor date: 202308 as this is when event data is from
df_cox <- df_cox %>% mutate(censor_date = c(ym(202308)))

#convert dates to correct format
df_cox2 <- df_cox
df_cox2$appt <- ym(str_sub(df_cox$appt, 1,7))
df_cox2$dt1_ym <- ym(df_cox$dt1_ym)
df_cox2$dod_ym <- ym(df_cox2$dod_ym)
sapply(df_cox2, function(x) sum(is.na(x)))

#check for any data anomalies
# - any CVD events after censor date?
df_cox2 %>% dplyr::filter(dt1_ym > censor_date) #0
# - any deaths prior to cvd event
df_cox2 %>% dplyr::filter(dod_ym < dt1_ym) #0
df_cox2 %>% dplyr::filter(composite_cardiac_outcome == "1" & dod_ym > censor_date) 
df_cox3 <- df_cox2 %>%
  mutate(composite_cardiac_outcome = ifelse(composite_cardiac_outcome == "1" & dod_ym > censor_date, "0", composite_cardiac_outcome),
          cardiac_death = ifelse(cardiac_death == "1" & dod_ym > censor_date, "0", cardiac_death))


df_cox2 <- df_cox3
df_cox2 <- df_cox2 %>%
  mutate(event_date = pmin(censor_date, dt1_ym, dod_ym, na.rm = TRUE)) #event date is the earliest out of diagnosis/death/censor
sapply(df_cox2, function(x) sum(is.na(x)))

#do some checks
noncvd <- df_cox2 %>% dplyr::filter(dead == "1" & cardiac_death == "0") #
cvd_anddeath <- df_cox2 %>% dplyr::filter(cvd == "1" & cardiac_death == "1") #

#Calculate time to event
df_cox2 <- df_cox2 %>%
  mutate(tte = difftime(event_date, appt, units = "days" ))
#check no negative results
df_cox2 %>% dplyr::filter(tte < 0) #0


#create time to event in years
df_cox2 <- df_cox2 %>%
  mutate(tte_years = as.numeric(tte)/365.25)

#set status as censored or event
df_cox2$composite_cardiac_outcome <- replace(df_cox2$composite_cardiac_outcome, is.na(df_cox2$composite_cardiac_outcome), "0")
df_cox2 <- df_cox2 %>%
  mutate(status = ifelse(composite_cardiac_outcome == "1" & event_date < censor_date, "1", "0"))
outcomes <- df_cox2 %>% dplyr::filter(composite_cardiac_outcome == "1") 

#check for missing values
sum(is.na(df_cox2$tte_years)) 

#save out file
fwrite(df_cox2, file = paste0(here(data_path), "/cox_df.csv"))
df_cox2 <- fread(paste0(data_path,"/cox_df.csv"))

df_cox2 %>% dplyr::filter(dead == "1" & cardiac_death == "0" & dod_ym < censor_date & event_date == dod_ym) %>% dim()


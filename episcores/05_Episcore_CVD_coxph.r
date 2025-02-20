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
# CVD event = composite cvd outcome (for now) = cvd-related diagnosis (stroke, CHD-NOS and myocardial infarction) OR cvd-related death (also cvd-related procedure?)
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
  dplyr::filter(Disease %in% cvd_list & Source == "Secondary_Care" & incident == "0" ) #160 cases, 116 prevalent cvd, 2 no baseline appt

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
deaths1 %>% dplyr::filter(id %in% test_target$id & if_any(everything(), ~. %in% c(heart_deaths, hypertension_deaths, cerebrovascular_deaths))) %>% dim() #42

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
sapply(df_cox, function(x) sum(is.na(x)))  #only missing data is for those without cvd 

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
df_cox2 %>% dplyr::filter(composite_cardiac_outcome == "1" & dod_ym > censor_date) #1 individual
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

#________Run Cox Models_________________________________________________________
###EpiScores

#add in episcores
episcore_df <- fread("/data_dir/ewas/bayesr/episcore/episcore_protein_predictions_0407.csv")
episcore_df[1:5, 1:5]

episcore_df <- transpose(episcore_df, keep.names = "protein", make.names = "protein")
episcore_df$protein <- gsub("X", "", episcore_df$protein)
episcore_df <- episcore_df %>% dplyr::rename(Sample_Sentrix_ID = protein)

#rank transform episcores (using same method as for rank transformation of protein data)
rint <- function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
episcore_df <- as.data.frame(episcore_df)
episcore_df[1:5, 1:5]
episcore_df <- episcore_df %>% modify_if(., is.numeric, rint)

#scale the episcores
episcore_df[,2:134] <- scale(episcore_df[,2:134])
df_cox2$status <- as.numeric(df_cox2$status)

#Combine with cox data frame
episcore_coxdf <- df_cox2 %>% dplyr::left_join(episcore_df, by = "Sample_Sentrix_ID")

#Look at CoxPH models for EpiScores with R > 0.1 and p < 0.05
sig_episcores <- fread(paste0(data_path, "significant_episcores.csv"))
protein_list <- sig_episcores$protein


#Run cox models and extract model output and assumption tests
assump_test_episcore <- list()
cox_results <- list()
cox_raw <- list() #for non-tidy output
aic_results <- list()
for (i in protein_list) {
formula_str <- paste("Surv(tte_years, status) ~", i, "+ age + as.factor(sex)")
formula <- as.formula(formula_str)
model <- coxph(formula, data = episcore_coxdf)
cox_raw[[i]] <- model
cox_results[[i]] <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
assump_test_episcore[[i]] <- cox.zph(model)
aic_results[[i]] <- extractAIC(model)
}
# head(cox_results)
# head(assump_test_episcore)
# head(aic_results)

#combine model results into dataframe
cox_df <- bind_rows(cox_results, .id = "protein")
head(cox_df)
fwrite(cox_df, file = here("ewas/bayesr/episcore/epi_coxresults.csv"))


#_______________________________________________________________________________________________________
#Summarise results

#Create summary table of results
p_thresholds <- c(0.05, 0.05/112) #112 EpiScores/Proteins included

# Function to summarize results at different p-value thresholds
summarize_cox_results <- function(data, thresholds) {
  summary_list <- lapply(thresholds, function(thresh) {
    subset_data <- data[data$p.value < thresh, ]
    subset_data <- subset_data %>% filter(protein == term)
    summary <- data.frame(
      threshold = thresh,
      count = nrow(subset_data),
      mean_HR = mean(subset_data$estimate, na.rm = TRUE)
    )
    return(summary)
  })
  summary_df <- do.call(rbind, summary_list)
  return(summary_df)
}

# Summarize results
summary_results <- summarize_cox_results(cox_df, p_thresholds)
print(summary_results) 
fwrite(summary_results, here("ewas/bayesr/episcore/summary_epiresults.csv"))

#________Run Cox Models_________________________________________________________
##Proteins - run for proteins in episcore list

#Compare to true protein levels
scaled_pheno[1:5, 1:5]
df_cox2$id <- as.character(df_cox2$id)
scaled_pheno$id <- as.character(scaled_pheno$id)
scaled_pheno <- scaled_pheno %>% dplyr::select(-c(age, sex))
true_test <- df_cox2 %>% left_join(scaled_pheno, by = "id")
true_test <- true_test %>% dplyr::filter(protein %in% episcore_list)

assump_test <- list()
cox_true_results <- list()
aic_protein <- list()

for (protein in protein_list) {
formula_str <- paste("Surv(tte_years, status) ~", protein, "+ age + as.factor(sex)")
formula <- as.formula(formula_str)
model <- coxph(formula, data = true_test)
cox_true_results[[protein]] <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
assump_test[[protein]] <- cox.zph(model)
aic_protein[[protein]] <- extractAIC(model)
}
head(cox_true_results)
head(assump_test)
head(aic_protein)

true_cox_df <- bind_rows(cox_true_results, .id = "protein")

#Create summary table of results
p_thresholds <- c(0.05, 0.05/112)

# Function to summarize results at different p-value thresholds
summarize_cox_results <- function(data, thresholds) {
  summary_list <- lapply(thresholds, function(thresh) {
    subset_data <- data[data$p.value < thresh, ]
    subset_data <- subset_data %>% filter(protein == term)
    summary <- data.frame(
      threshold = thresh,
      count = nrow(subset_data),
      mean_HR = mean(subset_data$estimate, na.rm = TRUE)
    )
    return(summary)
  })
  summary_df <- do.call(rbind, summary_list)
  return(summary_df)
}

# Summarize results
summary_results <- summarize_cox_results(true_cox_df, p_thresholds)
print(summary_results)
fwrite(summary_results, here("ewas/bayesr/episcore/summary_protresults.csv"))

#combine two result dataframes together
head(cox_df)
episcore_cox <- cox_df %>% mutate(predictor = c("Episcore"))
head(episcore_cox)
true_cox_df <- bind_rows(cox_true_results, .id = "protein")
true_cox <- true_cox_df %>% mutate(predictor = c("True"))

combined <- rbind(episcore_cox, true_cox)
fwrite(combined, file = paste0(here(data_path), "/cox_results_combined.csv"))



#__________proteins____________________________________________________________
#Checking proportional hazards
extract_matrix <- function(x) {
  return(as.data.frame(x$table)) #function to extract dataframes from cox.zph lists 
}

#for true protein levels
assump_df <- lapply(assump_test, extract_matrix) #convert .zph outputs to dataframes

prophaz <- bind_rows(assump_df, .id = "protein") # combine to one dataframe
prophaz <- data.table(prophaz, keep.rownames = TRUE)
prophaz2 <- prophaz %>% group_by(protein) %>% 
  dplyr::filter(rn == protein | grepl("GLOBAL", rn)) #retain details for protein & global 
head(prophaz2)

#identify models failing prop. hazards test
prophaz_protein_sig <- prophaz2 %>% dplyr::filter(p < 0.05) %>% dplyr::select(protein)
prophaz_protein <- unique(prophaz_protein_sig$protein)
prophaz2 %>% dplyr::filter(p < 0.05/112) %>% dim() #0
fwrite(prophaz2, file = paste0(here(data_path), "/cox_prophaz_true.csv"))

prophaz_failed <- prophaz2 %>% 
  dplyr::filter(p < 0.05) %>%
  dplyr::select(rn, protein, p) %>%
  dplyr::rename(term = "rn")
length(unique(prophaz_failed$protein)) 

fwrite(prophaz_failed, file = paste0(here(data_path), "/cox_prophaz_failedtrue.csv"))
prophaz_protein <- unique(prophaz_failed$protein)

#__________EpiScores____________________________________________________________
assump_df_epi <- lapply(assump_test_episcore, extract_matrix) #convert .zph outputs to dataframes

prophaz_epi <- bind_rows(assump_df_epi, .id = "protein")
prophaz_epi <- data.table(prophaz_epi, keep.rownames = TRUE)
prophaz_epi <- prophaz_epi %>% group_by(protein) %>% dplyr::filter(rn == protein | grepl("GLOBAL", rn))
head(prophaz_epi)

#identify models failing prop. hazards test
prophaz_epi_sig <- prophaz_epi %>% dplyr::filter(p < 0.05) %>% dplyr::select(protein)
prophaz_epilist <- unique(prophaz_epi_sig$protein)

prophaz_failed <- prophaz_epi %>% 
  dplyr::filter(p < 0.05) %>%
  dplyr::select(rn, protein, p) %>%
  dplyr::rename(term = "rn")
length(unique(prophaz_failed$protein)) 

fwrite(prophaz_failed, file = paste0(here(data_path), "/cox_prophaz_failedepi.csv"))
prophaz_protein <- unique(prophaz_failed$protein)
prophaz_epi %>% dplyr::filter(p < 0.05) %>% dim() 
prophaz_epi %>% dplyr::filter(p < 0.05/112) %>% dim() 
fwrite(prophaz_epi, file = paste0(here(data_path), "/cox_prophaz_episcore.csv"))


#_______________________________________________________________________________________
#Explore results
combined <- fread(paste0((data_path), "/cox_results_combined.csv"))

#count those significant at nominal significance threshold
combined %>%
  dplyr::filter(p.value < 0.05) %>%
  aggregate(. ~predictor, summary) #52 Episcore, 21 protein

combined %>%
  dplyr::filter(p.value <0.05/112) %>%
  aggregate(. ~predictor, summary) #17 Episcore, 1 protein
fwrite(combined, here("ewas/bayesr/episcore/select_coxresults.csv"))

#evaluate strength of associations for positive results
combined2 <- combined %>%
  dplyr::select(protein,estimate,predictor) %>%
  pivot_wider(names_from = predictor,
              values_from = estimate) %>%
  mutate(relationship = ifelse(
                          Episcore > 1 & Episcore > True, "stronger_pos",
                          ifelse(
                              Episcore < 1 & Episcore < True, "stronger_neg",
                                "weaker"
    )
  ))
table(combined2$relationship)
#_________________________________________________________________________________


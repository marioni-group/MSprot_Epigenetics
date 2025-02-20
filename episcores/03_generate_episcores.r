#Episcore Analysis

#Summary of data-prep for input to this:
#Data split into training and test subsets:

#______ Test = wave3, n = 3463
#______ Training = Remainder - any related to wave3, n = 6816
#______ For training set: Methylation M-values pre-regressed for age, sex, batch, residuals taken, then scaled
#______ For test set: Methylation M-values scaled only
#______ Protein data (matched to training or test subsets): pre-regressed for age, sex, logbmi, logpackyears, residuals taken, then scaled
#Then training subsets run through bayesR to generate CpG weights for each protein

#Load required libraries
library(tidyverse)
library(data.table)
library(ggpmisc)
library(broom)
library(plyr)
library(rstatix)
library(forcats)
library(here)
library(readr)

#Look at step2 output from gmrm-omi - contains effect sizes for each CpG in --out-dir /home_dir/gmrm_scores/output/${A}/
# #Load contents of mlma file
here()
data_path <- "/output/"

#generate list of files from step 2 gmrm-omi
file_list <- list.files(data_path, recursive = TRUE, pattern = "*.mlma", full.names = TRUE)

#Run through files to create list of dfs

results_list <- list() #initiate emtpy list

for (file in file_list) {
  filename <- tools::file_path_sans_ext(basename(file))
  data <- fread(file)
  results_list[[filename]] <- data
}
head(results_list)


results_list <- purrr::map(results_list, ~dplyr::rename(., marker_no = V1,
                effect_size = V2,
                pip = V3))

#add in CpGs
datadir <- "data_dir/"
CpG_order <- fread(here(datadir, "methylation_data", "cpg_order_training.csv"))
CpG_order <- CpG_order %>%
  mutate(marker_no = c(0:752721)) #marker should start from 0 
results_list <- purrr::map(results_list, ~ dplyr::left_join(., CpG_order, by = "marker_no"))

results_df <- bind_rows(results_list, .id = "protein")
head(results_df)

#save out results as csv and df
saveRDS(results_list, file = "data_dir/episcore/training_effect_sizes.rds")
fwrite(results_df, file = "/data_dir/episcore/training_effect_sizes.csv")
results_df <- fread(here(datadir, "training_effect_sizes.csv"))

#_____________________________________________________________________________________________________________________________
### Make EpiScores###

#Set up for creating predicted protein levels using the CpG weights created during training step
# 1. Load CpG weights and test methylation data
# 2. Check CpG weight orders and CpGs in methylation data match
# 3. For each individual, create predicted protein levels using the weights x CpG methylation values
# 4. Compare with actual protein data: pre-regressed/scaled, true, true + co-variates

#1. Load Data
#test set methylation data is in directory below (as prepared in episcore_data_prep script): 
datadir 
methdir <- "methylation_data/test_data"

#Load test methylation data
test_meth <- read_csv(here(datadir, methdir, "test_scaled_meth_data.csv")) 

#Load training effect sizes
results_df <- fread(here(datadir, "training_effect_sizes.csv"))

# 2. Check orders
#Check orders of CpGs between weights and test data, 
#transpose methylation data, so CpGs are rows
test_meth1 <- transpose(test_meth, keep.names = "rn", make.names = "rn")
test_meth1[1:5, 1:5]

#Create episcore data-frame
weights_df <- results_df %>% dplyr::select(CpGs, protein, effect_size)
weights_wide <- pivot_wider(weights_df,
                            names_from = protein,
                            values_from = effect_size)
identical(test_meth1$rn, weights_wide$cpgs) #check cpg orders of weights and methylation data frame match
#TRUE

#Check SSID orders match & extract protein list
protein_list <- unique(weights_df$protein)
SSID_list <- colnames(test_meth1)
SSID_list <- SSID_list[-1] #remove 'rn'


#extra analysis - how many CpGs have weights >0 or <0 for each protein?
included_cpgs <- weights_df %>%
  dplyr::filter(effect_size != 0)

episcore_count_cpgs<- included_cpgs %>% dplyr::count(protein, sort = TRUE)
fwrite(episcore_count_cpgs, file = "/data_dir/bayesr/cpgs_byepiscore.csv")

#Calculate EpiScores_______________________________________________
#Function to calculate episcore (multiplies CpGs by CpG weights)
calculate_episcore <- function (x,y) { 
  sum(x*y)
}

#Loop through all proteins
#Methylation data should have rows as cpgs, SSIDs as columns
input_meth <- test_meth1[,-1] #remove cpg column for input to calculation

episcore_results <- list()  #initialise empty list

for (protein in protein_list) {
 y <- weights_wide[[protein]] #extracts CpG weights for protein 
 print(protein)
episcore_results[[protein]] <- purrr::map(input_meth, ~(calculate_episcore(.x,y)),
                     .id = "SSID") #
}

episcore_df <- episcore_results %>% lmap(bind_rows, .id = "protein") %>% bind_rows()  
fwrite(episcore_df, file = here(datadir, "episcore_protein_predictions_0407.csv"))
#____________________________________________________________________________________________________________________




#Full ewas results

#Load required packages
library(tidyverse)
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
library(limma)
library(data.table)
library(qqman)
library(QCEWAS)
library(purrr)
library(readr)

#-------------------LOAD DATA - FUll RUN - BASIC MODEL -------------------------------#
#This should list all files in batches folder and subfolders

file_list <- list.files(path = "/data_dir/ewas/output_files/basic/fast_linear/", pattern = "*.linear", full.names = TRUE)
results_list <- list()

for (file in file_list) {
  filename <- tools::file_path_sans_ext(basename(file))
  data <- read_delim(file, delim ="\t")
  results_list[[filename]] <- data
}
head(results_list)

#_______________________________________________________________________________________________________
#Function to summarise counts
p_bonf <- 3.6e-8/133 #set Bonferroni-correct p value

results_list <- results_list[!grepl("\\.", names(results_list))] #select uniquely mapped proteins
results_filtered <- purrr::map(results_list, ~ dplyr::filter(., p < p_bonf)) #list of significant results only

#function to summarise osca results (to be run on the results from each model)
summarise_osca_results <- function(results_list)  {
 results_filtered <- purrr::map(results_list, ~ dplyr::filter(., p < p_bonf)) #filters to bonferroni significance threshold
 results_filtered_df <- bind_rows(results_filtered, .id = "protein") #creates df of significant results
 summary <- data.frame(
  model = model,
  dataset = "unique_proteins",
  unique_proteins = length(unique(results_filtered_df$protein)), #number proteins with assoc.
  unique_CpGs = length(unique(results_filtered_df$Probe)), #number proteins
  Total_significant = results_filtered_df %>% nrow()) #total significant results
}

#_______________________________________________________________________________________________________

model <- "basic" #set model type
basic_summary <- summarise_osca_results(results_list)

#write out these results
saveRDS(results_filtered_df, file = "/data_dir/ewas/output_files/basic/summary/fl_pbonf_significant_results.RDS")
fwrite(basic_summary, file = "/data_dir/ewas/output_files/basic/summary/basic_summary_133pbonf.csv")

#_______________________________________________________________________________________________________
#Function to summarise results by protein 

summary_byprot <- function(results_list) {

#summarise significant p-values for each protein: No. significant probes, p-value range,
no_sig_cpgs <- purrr::map(results_filtered, 
  ~nrow(.)) #generates count of number of significant cpgs per protein

#make dataframe of count of significant cpgs per protein
sig_cpgs <- bind_rows(no_sig_cpgs, .id = "protein") %>% #create df from list
pivot_longer(everything(),  #pivot so protein name becomes new column
  names_to = "protein", 
  values_to = "count_sig_CpGs")

#Calculate lambda for each protein
lambda <- purrr::map(results_list, 
                        ~dplyr::select(., p) %>%
                         P_lambda()) 

lambda_df <- bind_rows(lambda, .id = "protein") %>% #create df from list
pivot_longer(everything(),  #pivot so protein name becomes new column
  names_to = "protein", 
  values_to = "Lambda")

#Combine Lambda and sig cpgs
sig_cpgs <- merge(sig_cpgs, lambda_df, by = "protein")

#get summary numbers
sig_cpgs <- dplyr::arrange(sig_cpgs, desc(count_sig_CpGs))
return(sig_cpgs)
}
#_______________________________________________________________________________________________________
#Run on basic model

byprot_df <- summary_byprot(results_list)
byprot_df <- byprot_df %>%
  mutate(model = "basic")

#rename to retain for future work
basic_byprot_df <- byprot_df
fwrite(byprot_df, "/data_dir/ewas/output_files/basic/summary/basic_byprot_df_133pbonf.csv")
#-----------------------------------------------------------#
#MODEL2

file_list <- list.files(path = "/data_dir/ewas/output_files/model_2/fast_linear/", pattern = "*.linear", full.names = TRUE)
results_list <- list()

for (file in file_list) {
  filename <- tools::file_path_sans_ext(basename(file))
  data <- read_delim(file, delim ="\t")
  results_list[[filename]] <- data
}
head(results_list)
results_list <- results_list[!grepl("\\.", names(results_list))]
results_filtered <- purrr::map(results_list, ~ dplyr::filter(., p < p_bonf))
results_filtered_df <- bind_rows(results_filtered, .id = "protein")
fwrite(results_filtered_df, file = "/data_dir/ewas/output_files/model_2/summary/pbonf_results.csv")
model <- "model_2"

model2_summary <- summarise_osca_results(results_list)
model2_byprot_df <- summary_byprot(results_list)

fwrite(model2_summary, "/data_dir/ewas/output_files/model_2/summary/model2_summary_133pbonf.csv")
fwrite(model2_byprot_df, "/data_dir/ewas/output_files/model_2/summary/model2_byprot_summary_133pbonf.csv")

#-----------------------------------------------------------#
#MODEL3

file_list <- list.files(path = "/data_dir/ewas/output_files/model_3/fast_linear/", pattern = "*.linear", full.names = TRUE)
results_list <- list()

for (file in file_list) {
  filename <- tools::file_path_sans_ext(basename(file))
  data <- read_delim(file, delim ="\t")
  results_list[[filename]] <- data
}
head(results_list)
results_list <- results_list[!grepl("\\.", names(results_list))]

results_filtered <- purrr::map(results_list, ~ dplyr::filter(., p < p_bonf))
results_filtered_df <- bind_rows(results_filtered, .id = "protein")
fwrite(results_filtered_df, file = "/data_dir/ewas/output_files/model_3/summary/pbonf_results.csv")

model <- "model_3"

model3_summary <- summarise_osca_results(results_list)
model3_byprot_df <- summary_byprot(results_list)
head(model3_summary)
head(model3_byprot_df)

fwrite(model3_summary, file = "/data_dir/ewas/output_files/model_3/summary/model3summary_133pbonf.csv")
fwrite(model3_byprot_df, file = "/data_dir/ewas/output_files/model_3/summary/model3summary_byprot_133pbonf.csv")

#-----------------------------------------------------------#
#MODEL4

file_list <- list.files(path = "/data_dir/ewas/output_files/model_4/fast_linear/", pattern = "*.linear", full.names = TRUE)
results_list <- list()

for (file in file_list) {
  filename <- tools::file_path_sans_ext(basename(file))
  data <- readr::read_delim(file, delim ="\t")
  results_list[[filename]] <- data
}
head(results_list)
results_list <- results_list[!grepl("\\.", names(results_list))] #used for filtering to 133 proteins
results_filtered <- purrr::map(results_list, ~ dplyr::filter(., p < p_bonf))  #for 133
fwrite(results_filtered_df, file = "/data_dir/ewas/output_files/model_4/summary/133pbonf_results.csv")


model <- "model_4"

model4_summary <- summarise_osca_results(results_list)
model4_byprot_df <- summary_byprot(results_list)
head(model4_summary)
head(model4_byprot_df)

model4_byprot_df$protein <- gsub("_PC", "", model4_byprot_df$protein) #clean up names


fwrite(model4_summary, file = "/data_dir/ewas/output_files/model_4/model4summary_133pbonf.csv")
fwrite(model4_byprot_df, file = "/data_dir/ewas/output_files/model_4/summary/model4summary_byprot_133pbonf.csv")

#_______________________________________________________________________________________________
#Combine all summary results and write out
df <- rbind(basic_summary, model2_summary, model3_summary, model4_summary)
head(df)
fwrite(df, "/data_dir/ewas/output_files/model_4/summary/full_summary_133pbonf.csv")



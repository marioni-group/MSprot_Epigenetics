#BayesR output processing parallel

library(tidyverse)
library(data.table)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(biomaRt)
library(UpSetR)
library(RColorBrewer)

#Step1_______________________________________________

#Read in step 1 output files & remove test or repeated files
file_list <- list.files(path = "/home_dir/gmrm/output", recursive = TRUE, pattern = "*.csv", full.names = TRUE)
results_list <- list()

#extract data from files in file-list & create a list of dataframes
for (file in file_list) {
  filename <- tools::file_path_sans_ext(basename(file))
  data <- fread(file)
  results_list[[filename]] <- data
}
head(results_list)

#Create dataframe of all results
results_df <- bind_rows(results_list, .id = "protein")
head(results_df)

#rename output columns
results_df <- results_df %>% 
  dplyr::rename(Iteration = V1,
         Trait_indicator = V2,
         Variance_explained = V3, # SigmaG
         Residual_variance = V4, # SigmaE
         Heritability = V5,
         No_incl_markers = V6,
         Mixture_group_indicator = V7,
         No_prior_mixture_comp = V8,
         prob_prior_group_1 = V9,
         prob_prior_group_2 = V10,
         prob_prior_group_3 = V11,
         prob_prior_group_4 = V12)
head(results_df)

#Calculate mean heritability on 1250 iterations
heritability <- results_df %>% 
  dplyr::filter(., Iteration >750) %>%
  dplyr::select(Heritability, protein) %>%
  dplyr::group_by(protein) %>%
  summarise(mean_heritability = mean(Heritability),
          sd_heritability = sd(Heritability))  

#write out heritability results          
fwrite(heritability, file = "/data_dir/ewas/bayesr/output/new_step2/protein_herit_burnin.csv")

#write out full results
fwrite(results_df, file = "/data_dir/ewas/bayesr/output/step1_full_results.csv")


#Step2_______________________________________________________________________________________________

#.yest file = protein predictions
file_list <- list.files(path = "/home_dir/gmrm/output", recursive = TRUE, pattern = "*.yest", full.names = TRUE)
results_list <- list()

for (file in file_list) {
  filename <- tools::file_path_sans_ext(basename(file))
  data <- fread(file)
  results_list[[filename]] <- data
}
head(results_list)
results_df <- bind_rows(results_list, .id = "protein")
head(results_df)

predict_output <- results_df
fwrite(predict_output, file = "/data_dir/ewas/bayesr/output/step2_full_results.csv")


#mlma file =  effect sizes, pip etc
file_list <- list.files(path = "/home_dir/gmrm/output/", recursive = TRUE, pattern = "*.mlma", full.names = TRUE)
results_list <- list()

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
CpG_order <- fread("/data_dir/ewas/bayesr/order_files/CpG_order_bayesR.csv")
CpG_order <- CpG_order %>%
  mutate(marker_no = c(0:752721)) #marker should start from 0 
results_list <- purrr::map(results_list, ~ dplyr::left_join(., CpG_order, by = "marker_no"))
results_filtered <- purrr::map(results_list, ~dplyr::filter(., pip >0.95))#Select lead CpGs with PIP > 0.95

results_df <- bind_rows(results_list, .id = "protein")
head(results_df)


#Create df for lead CpGs
highpip <- bind_rows(results_filtered, .id = "protein")
Cpg_counts <- highpip %>% dplyr::count(protein, sort = TRUE)
protein_counts <- highpip %>% dplyr::count(CpGs, sort = TRUE)

#write these out:
fwrite(Cpg_counts, file = "/data_dir/ewas/bayesr/Cpg_counts_0310.csv")
fwrite(protein_counts, file = "/data_dir/ewas/bayesr/protein_counts_0310.csv")
fwrite(highpip, file ="/data_dir/ewas/bayesr/highpipcpgs_0310.csv")
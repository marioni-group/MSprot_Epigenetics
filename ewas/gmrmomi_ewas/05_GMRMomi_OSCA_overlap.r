#GMRMomi, OSCA results comparison

library(tidyverse)
library(data.table)
library(rstatix)

#Load GMRMomi results
all_anno <- fread("/data_dir/ewas/bayesr/output/new_step2/cis_trans_anno_trial_1604.csv")
head(all_anno)
dim(all_anno)

#Load osca results
osca_results <- fread("/data_dir/ewas/output_files/model_4/summary/fast_linear_summary/133pbonf_results.csv")
dim(osca_results)
head(osca_results)

osca_results <- osca_results %>%
  mutate(protein = gsub("_PC", "", protein)) #clean up data

#set up dfs for comparison
osca <- osca_results %>% dplyr::select(protein, Probe)
bayesr <- all_anno %>% dplyr::select(protein, CpGs) %>%
                          dplyr::rename(Probe = CpGs)

#identify common associations
common_associations <- dplyr::intersect(osca, bayesr) 
head(common_associations)
dim(common_associations) #505

#extract full osca details for overlap
osca_overlap_results <- left_join(common_associations, osca_results, by = c("protein", "Probe"))
dim(osca_overlap_results)

#extract full gmrmomi details for overlap
all_anno <- all_anno %>% dplyr::rename(Probe = CpGs)
bayesr_overlap_results <- left_join(common_associations, all_anno, by = c("protein", "Probe"))
colnames(bayesr_overlap_results) <- paste(colnames(bayesr_overlap_results), "bayes", sep = "_")
head(bayesr_overlap_results)
bayesr_overlap_results <- bayesr_overlap_results %>%
  dplyr::rename(protein = protein_bayes,
                Probe = Probe_bayes)

#set up to combine
colnames(osca_overlap_results) <- paste(colnames(osca_overlap_results), "osca", sep = "_")
osca_overlap_results <- osca_overlap_results %>%
  dplyr::rename(protein = protein_osca,
                Probe = Probe_osca)
head(osca_overlap_results)

#create df for both sets of results
common_results <- left_join(osca_overlap_results, bayesr_overlap_results, by = c("protein", "Probe"))
head(common_results)

#check effect directions
common_results2 <- common_results %>% 
  mutate(effect_direction = 
          ifelse(b_osca > 0 & effect_size_bayes > 0, "same", 
            ifelse(b_osca < 0 & effect_size_bayes < 0, "same", "different")))
head(common_results2)
common_results2 %>% dplyr::filter(effect_direction == "different") %>% nrow() #0 - #all effect directions are the same

#check correlation of effect sizes
cor(common_results2$b_osca, common_results2$effect_size_bayes) #0.9423316
cor.test(common_results2$b_osca, common_results2$effect_size_bayes, method = "pearson")

#save out
fwrite(common_results2, file = "/data_dir/ewas/output_files/model_4/summary/osca_overlap_results_pbonf133.csv" )
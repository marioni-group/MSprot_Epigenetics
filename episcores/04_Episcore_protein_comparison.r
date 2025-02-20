### Compare EpiScores to scaled protein data ###

#Load required libraries
library(tidyverse)
library(data.table)

#1. Compare with raw/scaled protein data________________________________________________________________________________________________
datadir <- "data_dir/"
test_target <- fread(here(datadir, "test_target.csv"))
protein_df <- readRDS(here("data", "GS_ProteinGroups_RankTransformed_23Aug2023.rds"))
episcore_df <- fread(here(datadir,"episcore_protein_predictions_0407.csv")) 

#subset protein_df to 133 unique proteins and ids for test group
protein_df <- protein_df %>% dplyr::select(id, protein_list)
protein_df <- protein_df %>% dplyr::filter(id %in% test_target$id)
ids <- test_target %>% dplyr::select(id, Sample_Sentrix_ID)
protein_df <- left_join(protein_df, ids, by = "id")
protein_df <- protein_df %>% dplyr::select(Sample_Sentrix_ID, everything())
protein_df <- protein_df %>% dplyr::select(-id)
# identical(protein_df$Sample_Sentrix_ID, test_target$Sample_Sentrix_ID)

#scale protein data
scaled_pheno <- protein_df
scaled_pheno[,2:134] <- scale(scaled_pheno[,2:134])

#some clearning of episcore_df
episcore_df1 <- data.table::transpose(episcore_df, keep.names = "protein", make.names = "protein")
episcore_df1$protein <- gsub("X", "", episcore_df1$protein)
episcore_df1 <- episcore_df1 %>% dplyr::rename(Sample_Sentrix_ID = protein)

#match orders
scaled_pheno <- scaled_pheno[match(episcore_df1$Sample_Sentrix_ID, scaled_pheno$Sample_Sentrix_ID),]
identical(scaled_pheno$Sample_Sentrix_ID, episcore_df1$Sample_Sentrix_ID)

#Prepare df for plotting/regression
true_df <- scaled_pheno
true_df <- true_df %>%
  mutate(type = c("True"))
true_df <- true_df %>% dplyr::select(Sample_Sentrix_ID, type, everything())
true_df1 <- true_df %>% pivot_longer(.,
    !c(Sample_Sentrix_ID, type),
    names_to = "protein",
    values_to = "value")  

predict_df <- episcore_df1
predict_df <- predict_df %>%
  mutate(type = c("Predicted"))
predict_df <- predict_df %>% dplyr::select(Sample_Sentrix_ID, type, everything())ß
predict_df1 <- predict_df %>% pivot_longer(.,
  !c(Sample_Sentrix_ID, type),
  names_to = 'protein',
  values_to = "value")

#Join the two 
plot_df <- rbind(true_df1, predict_df1)
plot_df1 <- plot_df %>% pivot_wider(.,
                                    names_from = "type",
                                    values_from = value)


#also add in covariates
cov <- test_target %>% dplyr::select(Sample_Sentrix_ID, age, sex, Batch)
plot_df2 <- left_join(plot_df1, cov, by = "Sample_Sentrix_ID")
plot_df2$sex <- as.factor(plot_df2$sex)
plot_df2$Batch <- as.factor(plot_df2$Batch)

#save out file
fwrite(plot_df2, file = here(datadir,)"comp_episcore_true_rawprotein_cov.csv")


#________________________________________________________________________________
#Compare EpiScores to measured, scaled proteins

plot_df2 <- fread(here(datadir, "comp_episcore_true_rawprotein_cov.csv"))
plot_df1 <- plot_df2 %>% dplyr::select(-c(age, sex, Batch))

#Function calculating r
cor_predictr <-   function(x){
  return(
          tidy(cor.test(x$True, x$Predicted, method = "pearson"))) 
}
correlation_results <- plyr::ddply(plot_df1, .(protein), cor_predictr)
head(correlation_results)
fwrite(correlation_results, file = paste0(datadir, "correlation_true_epi_2407.csv"))

#Filter to significantly correlated EpiScores
correlation_results <- correlation_results %>% arrange(desc(estimate))
corr_results_sig <- correlation_results %>% dplyr::filter(estimate > 0.1 & p.value < 0.05) #112

#save out list of significantly correlated episcores
fwrite(corr_results_sig, file = paste0(datadir, "significant_episcores.csv"))
sig_episcores <- corr_results_sig$protein

#plot correlations
level_order <- sig_episcores
p <- corr_results_sig %>%
  ggplot(aes(x = factor(protein, level = level_order), y = estimate)) +
  geom_point() +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw() +
  labs(x = "MS protein EpiScore", y = "Correlation with true protein level") +
  theme(axis.title.x = element_text(margin = margin(t = 30), size = 10),
        axis.title.y = element_text(margin = margin(r = 20), size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8)) +
  coord_fixed(ratio = 60)
ggsave("/data_dir/ewas/bayesr/episcore/episcore_correlation_plot_0310.png", height = 15, width = 30, units = "cm")
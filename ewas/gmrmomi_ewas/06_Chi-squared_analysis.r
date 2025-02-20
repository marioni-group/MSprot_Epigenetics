#Chi-squared analysis
library(tidyverse)
library(data.table)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#read in results 
anno <- fread("/data_dir/ewas/bayesr/output/new_step2/cis_trans_anno_trial_1604.csv")

#get annotations for full epic array
epic_anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
epic_anno <- data.frame(epic_anno)

#Load in CpGs from epic array used in analysis
CpG_order <- fread("/data_dir/ewas/bayesr/order_files/CpG_order_bayesR.csv")

#Chi-squared to compare positions relative to epic array positions
data <- table(anno$Relation_to_Island)
epic_anno <- setDT(epic_anno, keep.rownames = "probe")
epic_anno <- epic_anno %>% dplyr::filter(probe %in% CpG_order$CpGs)
epic <- table(epic_anno$Relation_to_Island)
chi_table <- rbind(data, epic)
chi_table <- chi_table %>% as.data.frame %>% cbind(rn = row.names(chi_table))
chi_table <- setDT(chi_table, keep.rownames = "source")
chi_table <- chi_table %>% dplyr::select(-rn)

chi_table
chi_table2 <- data.table::transpose(chi_table, keep.names = "source", make.names = "source")
chi_table2 <- chi_table2 %>% dplyr::rename(region = "source")
chi_table2 <- chi_table2 %>%
  mutate(exp_prop = epic/sum(epic),
         actual_prop = data/sum(data),
         expected_res = sum(data) * exp_prop)
chi_table2

#goodness of fit test with proportions
proportions <- chi_table2$exp_prop
observed <- chi_table2$data
expected_counts <- chi_table2$expected_res
chi_test <- chisq.test(observed, p = proportions)

# Calculate standardized residuals
std_residuals <- (observed - expected_counts) / sqrt(expected_counts)

# Calculate Z-scores for each category
z_scores <- (observed - expected_counts) / sqrt(expected_counts * (1 - proportions))

# Convert Z-scores to two-tailed p-values
p_values_z <- 2 * (1 - pnorm(abs(z_scores)))

# Adjust p-values for multiple comparisons (e.g., Bonferroni correction)
adjusted_p_values_z <- p.adjust(p_values_z, method = "bonferroni")

# Identify significant Z-scores after adjustment
significant_z <- which(adjusted_p_values_z < 0.05)

# Output significant Z-scores and their p-values
list(significant_z = significant_z,
     adjusted_p_values_z = adjusted_p_values_z[significant_z])

results_df <- data.frame(
  Category = names(data),
  No_bayesr_results = observed,
  Expected_Counts = round(expected_counts),
  Observed_Proportions = round(observed / sum(observed), 2),
  Expected_Proportions = round(proportions,2),
  Standardized_Residuals = round(std_residuals,2),
  z_scores = z_scores,
  P_Values = p_values_z,
  Adj_pvalues = adjusted_p_values_z,
  Significance = ifelse(adjusted_p_values_z < 0.05, "significant", "non-significant")
)
head(results_df)
fwrite(results_df, file = "/data_dir/ewas/bayesr/output/new_step2/chi_sq_results.csv")

#Set up data-frame for plot
head(chi_table)
chi_table <- chi_table %>%
  pivot_longer(., !source, names_to = "Relation_to_Island",
    values_to = "probes")

chi_table <- chi_table %>%
  dplyr::rename(origin = source) %>%
  mutate(origin = ifelse(origin == "data", "results", "epic_array"))

fwrite(chi_table, "/data_dir/ewas/bayesr/output/new_step2/df_for_chisq_plot.csv")
#Plot created in gmrmomi_summary_plots.r script
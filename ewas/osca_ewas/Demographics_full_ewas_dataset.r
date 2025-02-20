#Demographics for full EWAS dataset (n = 14,671)
library(tidyverse)
library(data.table)
library(skimr)
library(arsenal)

#data paths
ext_data <- "/data_dir/"

#EWAS: Full dataset
pheno <- readRDS("data/GS_phenos_internal.rds")
order <- read.csv("ewas/input/methylation_data/order_resid_scaled_methylation_data.csv")

#Filter to correct IDs & covariates
pheno <- pheno %>% dplyr::filter(id %in% target$Sample_Name) #Filter to IDs in dataset
pheno2 <- pheno %>% dplyr::select(id, age, sex, bmi, pack_years) #select covariates of interest
pheno2$sex <- as.factor(pheno2$sex)

#Generate summary data
pheno2 %>%select(-id) %>% summary() %>% broom::tidy()
table <- tableby(sex ~ age + bmi + pack_years, data = pheno2, test = FALSE, numeric.stats = c("median", "q1q3"))
tab1 <- as.data.frame(summary(table))
fwrite(tab1, file = "ewas/ewas_demographics.csv")

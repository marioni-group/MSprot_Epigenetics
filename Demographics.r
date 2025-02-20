#Demographics for write-up
library(tidyverse)
library(data.table)
library(finalfit)

#data paths
ext_data <- ""

#EWAS: Full dataset
target <- readRDS(paste0(ext_data, "/GS20k_Targets_18869.rds"))
pheno <- readRDS("data/GS_phenos.rds")
order <- read.csv("ewas/input/methylation_data/order_resid_scaled_methylation_data.csv")
target <- target %>% dplyr::filter(Sample_Sentrix_ID %in% order$order)
pheno <- pheno %>% dplyr::filter(id %in% target$Sample_Name)
pheno2 <- pheno %>% dplyr::select(id, age, sex, bmi, pack_years)
pheno2$sex <- as.factor(pheno2$sex)

#create summary table
table <- tableby(sex ~ age + bmi + pack_years, data = pheno2)
tab1 <- as.data.frame(summary(table))
fwrite(tab1, file = "ewas/ewas_demographics.csv")

#EpiScores: Training dataset
target <-  fread(here("episcore/training_target.csv"))
pheno_training <- pheno %>% filter(id %in% target$id)
table <- tableby(sex ~ age + bmi + pack_years, data = pheno_training)
tab2 <- as.data.frame(summary(table))
fwrite(tab2, file = "ewas/episcore_training_demographics.csv")

#EpiScores: Test dataset
test_target <-  fread(here("/episcore/test_target.csv"))
pheno_test <- pheno %>% filter(id %in% target$id)
table <- tableby(sex ~ age + bmi + pack_years, data = pheno_test)
tab3 <- as.data.frame(summary(table))
fwrite(tab3, file = "ewas/episcore_test_demographics.csv")


#CoxPH models - events etc
data_path <- here::here("ewas", "bayesr", "episcore")
df_cox2 <- df_cox2 <- fread(file = paste0(here(data_path), "/cox_df.csv"))

#Set-up data for table
forsummary <- df_cox2 %>%
  dplyr::select(age,sex,cvd,dead,cardiac_death,composite_cardiac_outcome,status, tte_years)

# Convert 0/1 variables to factors
forsummary$cvd <- replace(forsummary$cvd, is.na(forsummary$cvd), "0")
forsummary <- forsummary %>%
  mutate(across(where(~ all(. %in% c(0, 1))), ~ factor(., levels = c("0", "1"), labels = c("No", "Yes"))))

explanatory = c("age", "sex", "tte_years", "cvd", "cardiac_death")
dependent = "status" #cardiac event counted as event for coxph

#table construction
table2 <- forsummary %>%
    mutate(
        age = ff_label(age, "Age (years)"),
        sex = ff_label(sex, "Sex"),
        tte_years = ff_label(tte_years, "Time to event (years)"),
        cvd = ff_label(cvd, "Incident CVD"),
        cardiac_death = ff_label(cardiac_death, "CVD-related death"),
        composite_cardiac_outcome = ff_label(status, "CVD diagnosis or death")
    ) %>%
    summary_factorlist(dependent, explanatory, cont = "median",
                                         total_col = TRUE,
                                         column = FALSE,
                                         add_row_total = TRUE,
                                         include_row_missing_col = FALSE,
                                         add_dependent_label = TRUE, 
                                         dependent_label_prefix = "") %>%
                                         knitr::kable()

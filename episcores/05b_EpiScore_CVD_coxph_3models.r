#Run cox models with CVD covariates


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

#Load coxdf prepared in 05_Episcore_CVD_coxph.r
data_path <- ""
df_cox2 <- fread(paste0(data_path,"/cox_df.csv"))

#________Run Cox Models_________________________________________________________
###EpiScores

#add in episcores
episcore_df <- fread("episcore_protein_predictions_0407.csv")
episcore_df[1:5, 1:5]

episcore_df <- data.table::transpose(episcore_df, keep.names = "protein", make.names = "protein")
episcore_df$protein <- gsub("X", "", episcore_df$protein)
episcore_df <- episcore_df %>% dplyr::rename(Sample_Sentrix_ID = protein)

#filter to coxdf ids
dim(episcore_df)
dim(df_cox2)
episcore_df <- episcore_df %>% dplyr::filter(Sample_Sentrix_ID %in% df_cox2$Sample_Sentrix_ID)

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

#add covariates: smoking, bmi, alcohol units, SIMD rank, diabetes, HTN, BP, cholesterol
pheno <- readRDS("GS_phenos.rds") 
diseases <- fread(here("disease_codes_combined.csv")) #Disease phenotypes

#create disease data for diabetes and hypertension
diseases <- diseases %>% filter(Disease == "diabetes" | Disease == "hypertension")

#keep only prevalent diseases
diseases <- diseases %>% filter(incident == 0)

#use secondary -care level data
diseases <- diseases %>% filter(Source == "Secondary_Care")

#turn into data suitable for merging with coxdf
diabdf <- diseases %>%
  filter(Disease == "diabetes") %>%
  mutate(diabetes = ifelse(Disease == "diabetes", 1, 0)) %>%
  select(id, diabetes)

htndf <- diseases %>%
  filter(Disease == "hypertension") %>%
  mutate(htn = ifelse(Disease == "hypertension", 1, 0)) %>%
  select(id, htn)


covs <- pheno %>%
  select(id, bmi, units, pack_years, HDL_cholesterol, Total_cholesterol, avg_sys, avg_dia)

dim(covs)
covs <- covs %>% filter(id %in% episcore_coxdf$id)
sapply(covs, function(x) sum(is.na(x)))

#                id               bmi             units        pack_years 
#                 0                16               254                 4 
#   HDL_cholesterol Total_cholesterol           avg_sys           avg_dia 
#                25                21                 3                 3 

#add disease outcomes
covs <- left_join(covs, diabdf, by = "id")
covs <- left_join(covs, htndf, by = "id")
hist(covs$bmi) 
hist(covs$pack_years)
hist(covs$units)
hist(covs$HDL_cholesterol)
hist(covs$Total_cholesterol)
hist(covs$avg_sys)
hist(covs$avg_dia)

#Transform skewed variables
covs <- covs %>%
  mutate(logbmi = log(bmi),
         logunits = log(units + 1),
         logsmok = log(pack_years + 1),
         diabetes = coalesce(diabetes, 0),
         htn = coalesce(htn, 0))

#save out before data transformations or imputation
fordemo <- left_join(df_cox2, covs, by = "id")
fwrite(fordemo, file = "updated_cox_cov_fortable.csv" )

#below have mostly normal distributions
hist(covs$HDL_cholesterol)
hist(covs$Total_cholesterol)
hist(covs$avg_sys)
hist(covs$avg_dia)

#Prep data for scaling
str(covs)
covs$id <- as.character(covs$id)
covs$htn <- as.character(covs$htn)
covs$diabetes <- as.character(covs$diabetes)
str(covs)

covs2 <- covs %>% mutate(across(where(is.numeric), scale))

#impute missing
library(VIM)
#remove extra cov columns
covs2 <- covs2 %>% select(-c(bmi, units, pack_years))

# > sapply(covs2, function(x) sum(is.na(x)))
#                id   HDL_cholesterol Total_cholesterol           avg_sys 
#                 0                25                21                 3 
#           avg_dia          diabetes               htn            logbmi 
#                 3                 0                 0                16 
#          logunits           logsmok 
#               254                 4 

covs3 <- kNN(covs2, k = 5)

#add to episcore_coxdf
episcore_coxdf$id <- as.character(episcore_coxdf$id)
coxdf <- left_join(episcore_coxdf, covs3, by = "id")

#remove extra impute columns
coxdf2 <- coxdf %>% select(!contains("_imp"))
# coxdf2 %>% filter(heart_disease_Y == "1") %>% select(status)

#save out coxdf2
fwrite(coxdf2, file = "episcore_coxdf_rev_cvdmodels.csv")

#create equivalent protein df
protdf <- readRDS("ProteinGroups_RankTransformed_23Aug2023.rds")

#filter to proteins with significant episcores
protdf <- protdf %>% select(c(id, all_of(episcore_list)))
protdf <- protdf %>% filter(id %in% coxdf2$id)
str(protdf)
protdf$id <- as.character(protdf$id)
protdf <- protdf %>% mutate(across(where(is.numeric), scale))

#create coxdf
protcoxdf <- coxdf2 %>% select(-all_of(episcore_list))
protcoxdf <- left_join(protcoxdf, protdf, by = "id")
fwrite(protcoxdf, file = "prot_coxdf_rev_cvdmodels.csv")

#____RUN COX MODELS______#
#Load coxdfs

protcoxdf <- fread("prot_coxdf_rev_cvdmodels.csv")
epicoxdf <- fread("episcore_coxdf_rev_cvdmodels.csv")


#set up list to iterate through
sig_episcores <- fread(paste0(data_path, "significant_episcores.csv"))
episcore_list <- sig_episcores$protein



#Loop to run 3 models for proteins and then episcores
cox_results <- list()
assump_test <- list()
aic_results <- list()

#Run cox models for proteins

for (protein in episcore_list) {
    print(protein)

  #model list  
formula_list <- list(
    paste("Surv(tte_years, status) ~", protein, "+ age + as.factor(sex)"),
    paste("Surv(tte_years, status) ~", protein, "+ age + as.factor(sex) + logbmi + logsmok + logunits"),
    paste("Surv(tte_years, status) ~", protein, "+ age + as.factor(sex) + logbmi + logsmok + logunits + diabetes + htn + HDL_cholesterol + Total_cholesterol + avg_sys + avg_dia"))

  # Convert formulas to formula objects
    formula_list <- lapply(formula_list, as.formula)
   
   # Model indexing
    model_index <- 1

    #initialise lists for models
    assump_results <- list()
    model_results <- list()
    aic_list <- list()

for (formula in formula_list) {
model <- coxph(formula, data = protcoxdf)
model_name <- paste("model", model_index, sep = "_")
model_results[[model_name]] <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = model_name)
assump_results[[model_name]] <- cox.zph(model)$table %>% as.data.frame() %>% mutate(model = model_name)
aic_list[[protein]] <- extractAIC(model) 

  # Increment model index
      model_index <- model_index + 1  
}
assump_test[[protein]] <- bind_rows(assump_results) #***** May be better to use coxph output for this instead
cox_results[[protein]] <- bind_rows(model_results) 
aic_results[[protein]] <- bind_rows(aic_list)
}

head(cox_results)
head(assump_test)
head(aic_protein)

#combine model results into dataframe
protresults <- bind_rows(cox_results, .id = "protein")
protresults %>% filter(protein == term & p.value < 0.05/112)  #2 protein models significant

assumpdf <- bind_rows(assump_test, .id = "protein")

assumpdf %>% filter(p < 0.05)
assumpdf$term <- rownames(assumpdf)

aic <- bind_rows(aic_protein, .id = "protein") #model 3 only?

#save
fwrite(protresults, file = "prot_cvd_coxmodels_0108.csv")
fwrite(assumpdf, file = "prot_cvd_coxassump_0108.csv")
fwrite(aic, file = "prot_cvd_coxaic_0108.csv")


#
#Run the same for episcores
#

#Loop to run 3 models for proteins and then episcores
cox_results <- list()
assump_test <- list()
aic_results <- list()

#Run cox models for proteins

for (protein in episcore_list) {
    print(protein)

  #model list  
formula_list <- list(
    paste("Surv(tte_years, status) ~", protein, "+ age + as.factor(sex)"),
    paste("Surv(tte_years, status) ~", protein, "+ age + as.factor(sex) + logbmi + logsmok + logunits"),
    paste("Surv(tte_years, status) ~", protein, "+ age + as.factor(sex) + logbmi + logsmok + logunits + diabetes + htn + HDL_cholesterol + Total_cholesterol + avg_sys + avg_dia"))

  # Convert formulas to formula objects
    formula_list <- lapply(formula_list, as.formula)
   
   # Model indexing
    model_index <- 1

    #initialise lists for models
    assump_results <- list()
    model_results <- list()
    aic_list <- list()

for (formula in formula_list) {
model <- coxph(formula, data = epicoxdf)
model_name <- paste("model", model_index, sep = "_")
model_results[[model_name]] <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = model_name)
assump_results[[model_name]] <- cox.zph(model)$table %>% as.data.frame() %>% mutate(model = model_name)
aic_list[[protein]] <- extractAIC(model) 

  # Increment model index
      model_index <- model_index + 1  
}
assump_test[[protein]] <- bind_rows(assump_results) 
cox_results[[protein]] <- bind_rows(model_results) 
aic_results[[protein]] <- bind_rows(aic_list)
}

head(cox_results)
head(assump_test)
head(aic_protein)

#combine model results into dataframe
epiresults <- bind_rows(cox_results, .id = "protein")
epiresults %>% filter(protein == term & p.value < 0.05/112)  #

assumpdf <- bind_rows(assump_test, .id = "protein")
assumpdf$term <- rownames(assumpdf)
assumpdf %>% filter(p < 0.05)

aic <- bind_rows(aic_protein, .id = "protein") 

#save
fwrite(epiresults, file = "epi_cvd_coxmodels_0108.csv")
fwrite(assumpdf, file = "epi_cvd_coxassump_0108.csv")
fwrite(aic, file = "epi_cvd_coxaic_0108.csv")

#Read in 
protresults <- fread("prot_cvd_coxmodels_0108.csv")
epiresults <- fread("epi_cvd_coxmodels_0108.csv")

#create combined df
protresults <- protresults %>% mutate(predictor = "measured")
epiresults <- epiresults %>% mutate(predictor = "episcore")

comb <- rbind(protresults, epiresults)

#2. Create plot df
df <- comb %>% dplyr::filter(term == protein) #remove age and sex terms

#save out for supplementary
# fwrite(df, file = "comb_coxresults_forsupp.csv")
# df <- fread("comb_coxresults_forsupp.csv")


#Plots
library(patchwork)
library(egg)


#add in protein annotations
anno <- fread("protein_annos_forpub_nodup.csv")
anno <- anno %>% filter(protein %in% df$protein)
df <- left_join(df, anno, by = "protein")

df <- df %>%
  mutate(plot_label = ifelse(gene_names_simpl == "", protein, gene_names_simpl)) #TF is duplicated (serotransferrin and transferrin)

df <- df %>%
  mutate(plot_label = ifelse(plot_label == "TF" & protein == "C9JB55", protein, plot_label))

#filter to either predictor with Bonferroni signif. 
df <- df %>% dplyr::group_by(protein) %>% 
  dplyr::filter(any(p.value[model == "model_1"] < 0.05 / 112)) %>% ungroup()

df <- df %>% group_by(predictor) %>%
  mutate(protein = forcats::fct_reorder(protein, estimate)) #reorder so plot is easier to read
df <- df %>% 
  mutate(significant = ifelse(p.value < 0.05/112, "yes", "no"))  

df <- df %>%
  mutate(plot_lg = ifelse(predictor == "episcore" & significant == "yes", "Significant_Episcore",
                          ifelse(predictor == "episcore" & significant == "no", "Non-significant_Episcore",
                                 ifelse(predictor == "measured" & significant == "yes", "Significant_protein", "Non-significant_protein"))))


paired <- RColorBrewer::brewer.pal(12, "Paired")
GnBu <- RColorBrewer::brewer.pal(9, "GnBu")
df$plot_lg <- factor(df$plot_lg, levels = c("Significant_Episcore", "Non-significant_Episcore", "Significant_protein", "Non-significant_protein"))
df <- df %>% mutate(plot_fill = plot_lg)


df <- df %>% group_by(predictor) %>%
  mutate(plot_label = forcats::fct_reorder(plot_label, estimate)) #reorder so plot is easier to read

combp <- ggplot(data=df, aes(y=plot_label, x=estimate, xmin=conf.low, xmax=conf.high, colour=plot_lg, shape=plot_lg)) +
  geom_point(position=position_dodge(width=0.75), size=2) + 
  geom_errorbarh(position=position_dodge(width=0.75), height=.1) +
  scale_color_manual(
    values = c(
      "Significant_Episcore" = paired[6],
      "Non-significant_Episcore" = paired[5],
      "Significant_protein" = GnBu[9],
      "Non-significant_protein" = GnBu[6]
    ),
    labels = c(
      "Significant EpiScore",
      "Non-significant EpiScore",
      "Significant protein",
      "Non-significant protein"
    ),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "Significant_Episcore" = 19,
      "Non-significant_Episcore" = 1,
      "Significant_protein" = 19,
      "Non-significant_protein" = 1
    )
  ) +
guides(
    colour = guide_legend("Model & Bonferroni Significance", 
                          override.aes = list(shape = c(19, 1))),  # Override shapes to match significance
    shape = "none") +
  labs(title='', x='Hazard Ratio [95% Confidence Interval]', y='Protein') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  coord_cartesian(xlim=c(0.5, 2)) +
  theme_minimal() + 
  theme(
    legend.position = "top",
    legend.background = element_rect(fill="white", size=.4),
    legend.title = element_text(size=14),
    legend.text = element_text(size=10),
    strip.text = element_text(size=14),
    panel.spacing = unit(2, "lines"),
    axis.text = element_text(size=12),
    axis.title = element_text(size=14),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.title.y = element_text(margin = margin(r = 20))
  ) +
  facet_wrap(~model, labeller = as_labeller(c(model_1 = "Model 1", model_2 = "Model 2", model_3 = "Model 3")))

ggsave(combp, file = "ewas/bayesr/episcore/revision_cvd_models/cvd_combmodel_hrplot_0108.pdf", width = 30, height = 30, units = "cm")



#Model checks
#EpiScores
epiassump <- fread("epi_cvd_coxassump_0108.csv")
epifailed <- epiassump %>% filter(p < 0.05) %>% select(protein) %>% unique()
episig <- epiresults %>% filter(protein == term & p.value < 0.05/112) %>% select(protein) %>% unique()
episig %>% filter(protein %in% epifailed$protein) #Q99459 - only significant in 1st model
epiresults %>% filter(protein == "Q99459" & protein == term)

#Proteins
protassump <- fread("prot_cvd_coxassump_0108.csv")
protfailed <- protassump %>% filter(p < 0.05) %>% select(protein) %>% unique()
protsig <- protresults %>% filter(protein == term & p.value < 0.05/112) %>% select(protein) %>% unique()
protsig %>% filter(protein %in% protfailed$protein) #P04004 - significant in first model only

#combine
head(protassump)
protassump <- protassump %>% mutate(predictor = "protein")
head(epiassump)
epiassump <- epiassump %>% mutate(predictor = "episcore")
cb <- rbind(epiassump, protassump)
cb <- cb %>% mutate(model_component = ifelse(df == 1, "protein", "global"))
fwrite(cb, file = "prot_cvd_combassump_1108.csv")
st <- cb %>%
  group_by(model, predictor, model_component) %>%
  filter(p < 0.05) %>%
  summarise(
    n_failed = n(),
    failed_uniprots = paste(unique(protein), collapse = ","),
    .groups = "drop"
  )

fwrite(st, file = "prop_hazard_summary.csv")
# st <- fread("prop_hazard_summary.csv")


#Create summary table of results
head(df)
df <- df %>% filter(protein %in% episcore_list)
df <- df %>% mutate(p_threshold = ifelse(p.value < 0.05/112, "0.05/112", 
                                    ifelse(p.value < 0.05 & p.value > 0.05/112, "0.05", "non-significant")))

# Define your p-value thresholds
thresholds <- c(0.05/112, 0.05)

# Create a summary table for multiple thresholds
st <- df %>%
  group_by(predictor, model) %>%
  summarize(
    "P < 0.05/112" = sum(p.value < thresholds[1], na.rm = TRUE),
    "P < 0.05" = sum(p.value < thresholds[2], na.rm = TRUE)
  ) %>%
  arrange(model, predictor)

# View the summary table
print(st)
fwrite(st, file = "summary_sig_results_table.csv")
st <- fread("summary_sig_results_table.csv")

#add in proportional hazards
epiassump <- fread("epi_cvd_coxassump_0108.csv")
epiassump <- epiassump %>% mutate(predictor = "episcore")
protassump <- fread("prot_cvd_coxassump_0108.csv")
protassump <- protassump %>% mutate(predictor = "protein")

cb <- rbind(epiassump, protassump)
cb$term <- gsub("\\...\\d+", "", cb$term)
cb2 <- cb %>% group_by(protein) %>% 
  dplyr::filter(term == protein | grepl("GLOBAL", term))
head(cb2)
tail(cb2)

fwrite(cb2, "prophaz_all_1108.csv")

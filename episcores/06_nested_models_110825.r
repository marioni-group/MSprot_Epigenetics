#Nested models - no age filtering, any model with nominal significance included#
library(data.table)
library(survival)
library(survminer)
library(here)
#for plot ordering
library(forcats)

#Look at nested models: ~ protein + age + sex vs. ~ episcore + protein + age + sex
# 1. Create list of coxph results for first model
# 2. Create list of coxph results for second model
# 3. Use anova.coxph function to compare results



#1. Create list of coxph results for first model
datadir <- ""
data_path <- here::here("")

#Create list of proteins/episcores to include
coxmodelresults <- fread("comb_coxresults_forsupp.csv") #complete results for unfiltered dataset
results <- coxmodelresults %>% 
  dplyr::filter(model == "model_1") #244


df2 <- results %>% dplyr::group_by(protein) %>%
  dplyr::filter(any(p.value < 0.05)) #n = 56

df <- results %>% dplyr::group_by(protein) %>% 
  dplyr::filter(all(p.value < 0.05)) #filter to both predictor with nominal signif
protein_list <- unique(df$protein) %>% unlist() #17 proteins
  

#Load prepared cox dataframe 
df_cox2 <- fread(here(datadir, "cox_df_2907.csv")) 
protein_df <- readRDS(here("GS_ProteinGroups_RankTransformed_23Aug2023.rds")) 


#subset protein_df to ids for test group
protein_df <- protein_df %>% 
  dplyr::select(id, protein_list) %>%
  dplyr::filter(id %in% df_cox2$id)


#scale protein data
scaled_pheno <- protein_df
scaled_pheno[,2:19] <- scale(scaled_pheno[,2:19])
scaled_pheno$id <- as.character(scaled_pheno$id)
df_cox2$id <- as.character(df_cox2$id)
true_test <- df_cox2 %>% left_join(scaled_pheno, by = "id")

#Run cox models
assump_test <- list()
cox_true_results <- list()
cox_models_1  <- list()

for (protein in protein_list) {
formula_str <- paste("Surv(tte_years, status) ~", protein, "+ age + as.factor(sex)")
formula <- as.formula(formula_str)
model <- coxph(formula, data = true_test)
cox_models_1[[protein]] <- coxph(formula, data = true_test)
cox_true_results[[protein]] <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
assump_test[[protein]] <- cox.zph(model)
}
head(cox_true_results)
head(assump_test)

#summarise model1 output
model1_df <- bind_rows(cox_true_results, .id = "protein")
model1_sig <- model1_df %>% 
  dplyr::filter(term == "estimate" & p.value <0.05/17) #0

#Run cox model with episcore as additional variable
#load and subset episcore_df
episcore_df <- fread(here(datadir, "episcore_protein_predictions_0407.csv"))
episcore_df <- data.table::transpose(episcore_df, keep.names = "protein", make.names = "protein")
episcore_df$protein <- gsub("X", "", episcore_df$protein)
episcore_df <- episcore_df %>% dplyr::rename(Sample_Sentrix_ID = protein)

#rank transform episcores (using same method as for rank transformation of protein data)
episcore_df <- episcore_df %>% dplyr::filter(Sample_Sentrix_ID %in% true_test$Sample_Sentrix_ID)
rint <- function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
episcore_df <- as.data.frame(episcore_df)
episcore_df[1:5, 1:5]
episcore_df <- episcore_df %>% modify_if(., is.numeric, rint)
episcore_df[,2:134] <- scale(episcore_df[,2:134])

identical(episcore_df$Sample_Sentrix_ID, true_test$Sample_Sentrix_ID) #TRUE
episcore_coxdf <- df_cox2 %>% dplyr::left_join(episcore_df, by = "Sample_Sentrix_ID")

assump_test_model2 <- list()
cox_model2_results <- list()
cox_models_2 <- list()

for (protein in protein_list) {
data_with_estimate <- true_test
data_with_estimate$estimate <-  episcore_coxdf[[protein]] #***need to check orders****
formula_str <- paste("Surv(tte_years, status) ~ estimate +", protein, "+ age + as.factor(sex)")
formula <- as.formula(formula_str)
model <- coxph(formula, data = data_with_estimate)
cox_models_2[[protein]] <- coxph(formula, data = data_with_estimate)
cox_model2_results[[protein]] <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
assump_test_model2[[protein]] <- cox.zph(model)
}
head(cox_model2_results)
head(assump_test_model2)

#summarise model2 output
model2_df <- bind_rows(cox_model2_results, .id = "protein")
nested_sig <- model2_df %>% 
  dplyr::filter(term == "estimate" & p.value  < (0.05/17)) #10

episcore_prot_model <- model2_df %>% dplyr::filter(term == "estimate") 
episcore_prot_model <- episcore_prot_model %>%
  mutate(predictor = "Full model") 
prot_model <- bind_rows(cox_true_results, .id = "protein")
prot_model <- prot_model %>% dplyr::filter(protein == term)
prot_model <- prot_model %>%
  mutate(predictor = "Protein only model")


#add them all together
nested_models <- rbind(episcore_prot_model, prot_model)


#compare models
anova_results <- list()
for (protein in protein_list) {
  anova_results[[protein]] <- stats::anova(
    cox_models_1[[protein]], cox_models_2[[protein]])
}
head(anova_results)

#tidy results
tidy_anova <- function(anova_result, i) {
  data.frame(
    Predictor = i,
    loglik1 = anova_result$loglik[1],
    loglik2 = anova_result$loglik[2],
    Chisq = anova_result$Chisq[2],
    Df = anova_result$Df[2],
    PValue = anova_result$`Pr(>Chi)`[2]
  )
}

tidy_results <- list()
for (protein in protein_list) {
tidy_results[[protein]] <- data.frame(
    Predictor = protein,
    loglik1 = anova_results[[protein]]$loglik[1],
    loglik2 = anova_results[[protein]]$loglik[2],
    Chisq = anova_results[[protein]]$Chisq[2],
    Df = anova_results[[protein]]$Df[2],
    PValue = anova_results[[protein]]$`Pr(>|Chi|)`[2])
}


head(tidy_results)
tidy_df <- bind_rows(tidy_results, .id = "protein")
# anova_sig <- tidy_df %>% dplyr::filter(PValue < 0.05/18) #18
fwrite(tidy_df, file = "anova_nested_models_110825.csv")
head(tidy_df)
tidy_df2 <- tidy_df %>%
  mutate(improved_fit = ifelse(PValue < 0.05/17, "yes", "no"))
head(tidy_df2)
table(tidy_df2$improved_fit)
models_with_impr_fit <- tidy_df2 %>% dplyr::filter(improved_fit == "yes") %>% dplyr::select(protein) %>% as.list()
models_with_impr_fit <- models_with_impr_fit$protein

fwrite(tidy_df2, file = "anova_nested_models_110825.csv")

#add to nested_models df
nested_models2 <- nested_models %>%
  mutate(improved_fit = ifelse(protein %in% models_with_impr_fit & predictor == "Full model", "yes", 
                               ifelse(predictor == "Protein only model", "protein",
                                     ifelse(predictor == "Episcore", "episcore", "no"))))

fwrite(nested_models2, file = "nested_models_110825.csv")

nested_models2 <- fread("nested_models_110825.csv")
#Using final fit to explore nested models
#First explore using example protein & Episcore

#1. Set up new data frame with selected proteins
length(protein_list) #17
epi_list <- tidy_df %>% dplyr::filter(PValue < 0.05/17) %>% dplyr::select(protein)
epi_list <- epi_list$protein

head(df_cox2)
identical(scaled_pheno$id, df_cox2$id) #FALSE
scaled_pheno <- scaled_pheno[match(df_cox2$id, scaled_pheno$id),]
identical(scaled_pheno$id, df_cox2$id) #TRUE

#nb prepped episcore df
identical(episcore_df$Sample_Sentrix_ID, df_cox2$Sample_Sentrix_ID) #TRUE


#2. set up dependent and explanatory variables for finalfit
library(finalfit)
explanatory <- c("episcore", "true","age", "sex")
explanatory_multi <- c("true", "age", "sex")
dependent <- "Surv(tte_years, status)"
# fitdf %>%
#   finalfit(dependent, explanatory, explanatory_multi, keep_models = TRUE)

finalfitlist <- list()
fitdf <- df_cox2
for (protein in protein_list) {
fitdf$true <- scaled_pheno[[protein]]
fitdf$episcore <- episcore_df[[protein]]
finalfitlist[[protein]] <- fitdf %>%
  finalfit(dependent, explanatory, explanatory_multi, keep_models = TRUE)
}  

head(finalfitlist)
finalfitlist_df <- bind_rows(finalfitlist, .id = "protein")

view(finalfitlist_df)
fwrite(finalfitlist_df, file = here("finalfit_nestedmodels_hrs_1108.csv"))



#Nested model plot:

nested_models2 <- fread("nested_models_110825.csv")
nested_models2 <- nested_models2 %>% mutate(
  plot_label = ifelse(predictor == "Full model" & improved_fit == "yes", "Full model - improved fit", 
                  ifelse(predictor == "Full model" & improved_fit == "no", "Full model - no fit improvement", "Protein only model")))

colours <- c("#420F75FF", "#E7A8FBFF", "#429AFCFF")

#plot HRs by predictor 


#Add gene-name annotations
anno <- fread("protein_annos_forpub_nodup.csv")
anno <- anno %>% filter(protein %in% nested_models2$protein)
anno <- anno %>% select(protein, gene_names_simpl)
nested_models2 <- left_join(nested_models2, anno, by = "protein")

nested_models2 <- nested_models2 %>%
 mutate(gene_names_simpl = forcats::fct_reorder(gene_names_simpl, estimate))



p1 <- ggplot(data=nested_models2, aes(y=gene_names_simpl, x=estimate, xmin=conf.low, xmax=conf.high, colour = plot_label)) +
  geom_point(position = position_dodge(width = 0.75)) + 
  geom_errorbarh(position = position_dodge(width = 0.75), height=.1) +
    scale_color_manual(values = c("Full model - improved fit" = "#420F75FF",
                                "Full model - no fit improvement"= "#420F7560",
                                "Protein only model"= "#429AFC60")) +
                                guides(
    colour = guide_legend("Model")) + 
  labs(title='', x='Hazard Ratio [95% Confidence Interval]', y = 'Protein') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  coord_cartesian(xlim = c(0.5, 2)) +
  theme_minimal() +
     theme(legend.position = c(.8, .2),
        legend.background = element_rect(fill="white", size=.15),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) 


ggsave(filename = "nested_models_120825.tiff", 
        plot = egg::set_panel_size(p=p1, width=unit(8, "cm"), height=unit(16, "cm")),
        compression = "lzw",
        dpi = 300)

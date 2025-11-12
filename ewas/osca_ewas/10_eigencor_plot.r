#PCA eigencor plot

#first use micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate R
R


library(PCAtools)
library(tidyverse)
library(data.table)
library(reshape2)
library(Hmisc)

#Load PCs
methPC <- fread("gs20k_pcs.csv")


#Select first 20 PCs
methPC[1:5, 1:5]
methPC <- methPC[,1:21]

#Load full OSCA covariates
cov <- fread("quant_cov.cov")
order <- read.csv("order_resid_scaled_methylation_data.csv")

#add age and sex
pheno <- readRDS("GS20k_Targets_18869.rds")

#Combine files
names(order) <- "ID"
cov <- cov %>% select(-FID) %>% rename(ID = IID)
df <- left_join(order, methPC, by = "ID")
df <- left_join(df, cov, by = "ID")
pheno <- pheno %>% select(Sample_Sentrix_ID, age, sex, Batch) %>% rename(ID = Sample_Sentrix_ID)

df2 <- left_join(df, pheno, by = "ID")
df3 <- df2 %>% select(-ID)

#create eigencor plot
df3 <- df3 %>%
  mutate(sex = as.factor(sex),
         Batch = as.factor(Batch))


#for now remove categorical data
df4 <- df3 %>% select(-c(sex, Batch))

#create separate dfs for pcs and covs
pcdf <- df4[1:20]
covsdf <- df4[21:29] 

#cor matrix
cor_matrix <- cor(pcdf, covsdf)
cor_matrix <- as.data.frame(cor_matrix)
setDT(cor_matrix, keep.rownames = TRUE)
fwrite(cor_matrix, file = "cor_matrix_pcs_covs.csv")


# Melt the correlation matrix for plotting
cor_melted <- melt(cor_matrix)
cor_melted$label <- round(cor_melted$value, 2)

# Plot using ggplot2
p1 <- ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = label), color = "black", size = 2) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1)) +
  theme_bw() +
  labs(title = "Eigencor Plot", x = "PCs", y = "Covariates", fill = "Correlation") +
  theme(axis.text = element_text(size = 7))
ggsave(p1, file = "eigencor.pdf", width = 20, height = 10, units = "cm")


#Now create an eigencorr plot with the proteins
protein_df <- readRDS("GS_ProteinGroups_RankTransformed_23Aug2023.rds")
pheno <- readRDS("GS20k_Targets_18869.rds")
cvs <- pheno %>% 
    select(Sample_Name, Sample_Sentrix_ID) %>%
    rename(id = Sample_Name,
           ID = Sample_Sentrix_ID) %>%
    filter(ID %in% df2$ID)

protein_df <- protein_df %>% filter(id %in% cvs$id)
protein_select <- names(protein_df)[!grepl("\\.", names(protein_df))]
protein_df <- protein_df[, protein_select]
protein_df$id <- as.character(protein_df$id)

df <- left_join(cvs, protein_df, by = "id")

#add PCs
df2 <- left_join(df, methPC, by = "ID")


#Make plot
#remove id columns
df3 <- df2 %>% select(-c(ID, id))

#create separate dfs for pcs and proteins
pcdf <- df3[134:153]
protdf <- df3[1:133] 

#cor matrix
cor_matrix <- cor(pcdf, protdf)

#save out for supplementary
cor_matrix[1:5, 1:5]
cor_matrix <- as.data.frame(cor_matrix)
setDT(cor_matrix, keep.rownames = TRUE)
fwrite(cor_matrix, file = "cor_matrix_pcs_proteins.csv")



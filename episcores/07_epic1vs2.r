#compare EPIC v1 CpGs (in training), EPIC v2 CpGs, Episcore CpGs

library(tidyverse)
library(data.table)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
epic_anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
epic2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
epic2 <- setDT(as.data.frame(epic2))
epic2 <- epic2 %>% mutate(on_v1 = ifelse(EPICv1_Loci == "", "no", EPICv1_Loci))
rownames(epic2) <- NULL
overlap <- intersect()

#CpGs from EPICv1 in EPICv2 = 727,162 = 78%

# CpGs from EPICv1 in training
traincpgs <- fread("cpgs_tokeep.txt", header=F)

overlap2 <- traincpgs %>% filter(V1 %in% epic2$EPICv1_Loci)
dim(overlap2) #686130 of the training CpGs are in the EPICv1_Loci on EPICv2
686130/752722
# [1] 0.9115317

686130/930075
# [1] 0.7377147    #74% of CpGs used in training are on the EPICv2 array

#Load EpiScore Weights
#Load training effect sizes
results_df <- fread(here("training_effect_sizes.csv"))


#Create episcore data-frame
weights_df <- results_df %>% dplyr::select(CpGs, protein, effect_size)
weights_wide <- pivot_wider(weights_df,
                            names_from = protein,
                            values_from = effect_size)


#Create df of CpGs for each EpiScore
#extra analysis - how many CpGs have weights >0 or <0 for each protein?
included_cpgs <- weights_df %>%
  dplyr::filter(effect_size != 0)

episcore_count_cpgs<- included_cpgs %>% dplyr::count(protein, sort = TRUE)
fwrite(episcore_count_cpgs, file = "cpgs_byepiscore.csv")
# episcore_count_cpgs <- fread("cpgs_byepiscore.csv")


#Look at proportion of these CpGs that are EPICv2 CpGs
included_cpgs <- included_cpgs %>%
 mutate(EPICv2_locus = ifelse(CpGs %in% epic2$EPICv1_Loci, "yes", "no"))

episcore_cpgs <- included_cpgs %>% 
    group_by(protein) %>%
    summarise(
        cpgs = n(),
        epic_v2loci = count(EPICv2_locus == "yes"))

episcore_cpgs <- episcore_cpgs %>%
  mutate(prop_epic2 = epic_v2loci/cpgs)
episcore_cpgs <- episcore_cpgs %>% arrange(prop_epic2)
fwrite(episcore_cpgs, file = "epicv1_epiv2_episcores.csv")


#__________________________________________________________________________________________-
#Create EpiScores with  EpiC1&2 Loci only

loci <- fread("epic1_epic2_cpgs.csv")
head(loci)
dim(loci)
head(weights_df)


weights_df2 <- weights_df %>%
  mutate(EPICv2_locus = ifelse(CpGs %in% epic2$EPICv1_Loci, "yes", "no"))
head(weights_df2)

weights_df2 <- weights_df2 %>%
  mutate(
    effect_size = ifelse(EPICv2_locus == "no", 0, effect_size)
  )
head(weights_df2)

weights_wide <- pivot_wider(weights_df2,
                            names_from = protein,
                            values_from = effect_size)
head(weights_wide)  

#Make EpiScores for loci in both Epic1 and Epic2 arrays:
#1. Load Data
#test set methylation data is in directory below (as prepared in episcore_data_prep script): 
datadir <- ""
methdir <- ""

#Load test methylation data
test_meth <- read_csv(here(datadir, methdir, "test_scaled_meth_data.csv"))
test_meth[1:5, 1:5]


# 2. Check orders
#Check orders of CpGs between weights and test data, 
#transpose methylation data, so CpGs are rows
test_meth1 <- transpose(test_meth, keep.names = "rn", make.names = "rn")
test_meth1[1:5, 1:5]

#Create episcore data-frame
identical(test_meth1$rn, weights_wide$CpGs) #check cpg orders of weights and methylation data frame match
#TRUE
weights_wide <- weights_wide[match(test_meth1$rn, weights_wide$CpGs),]
identical(test_meth1$rn, weights_wide$CpGs)
#TRUE

#Calculate EpiScores_______________________________________________
#Function to calculate episcore (multiplies CpGs by CpG weights)
calculate_episcore <- function (x,y) { 
  sum(x*y)
}

#Loop through all proteins
#Methylation data should have rows as cpgs, SSIDs as columns
input_meth <- test_meth1[,-1] #remove cpg column for input to calculation

episcore_results <- list()  #initialise empty list
protein_list <- colnames(weights_wide)
protein_list <- protein_list[-1]
protein_list <- protein_list[-1]

for (protein in protein_list) {
 y <- weights_wide[[protein]] #extracts CpG weights for protein 
 print(protein)
episcore_results[[protein]] <- purrr::map(input_meth, ~(calculate_episcore(.x,y)),
                     .id = "SSID") #
}

episcore_v1_2 <- episcore_results %>% lmap(bind_rows, .id = "protein") %>% bind_rows()  
fwrite(episcore_v1_2, file = "episcores_loci_v12_only.csv")
head(episcore_v1_2)


##Also run correlation with the protein data
protein_df <- readRDS(here("GS_ProteinGroups_RankTransformed_23Aug2023.rds"))
protein_df[1:5, 1:5]
datadir <- ""
test_target <- fread(here(datadir, "test_target.csv"))

# subset protein_df to 133 unique proteins and ids for test group
protein_df <- protein_df %>% dplyr::select(id, all_of(protein_list))
protein_df <- protein_df %>% dplyr::filter(id %in% test_target$id)
ids <- test_target %>% dplyr::select(id, Sample_Sentrix_ID)
protein_df <- left_join(protein_df, ids, by = "id")
protein_df <- protein_df %>% dplyr::select(Sample_Sentrix_ID, everything())
protein_df <- protein_df %>% dplyr::select(-id)
dim(protein_df)

#scale protein data
scaled_pheno <- protein_df
scaled_pheno[,2:134] <- scale(scaled_pheno[,2:134])
scaled_pheno[1:5, 1:5]
head(episcore_v12)
dim(episcore_v12)
episcore_v12[,2:134] <- scale(episcore_v12[,2:134])
dim(scaled_pheno)
identical(scaled_pheno$Sample_Sentrix_ID, episcore_v12$Sample_Sentrix_ID) #FALSE
scaled_pheno <- scaled_pheno[match(episcore_v12$Sample_Sentrix_ID, scaled_pheno$Sample_Sentrix_ID),]
identical(scaled_pheno$Sample_Sentrix_ID, episcore_v12$Sample_Sentrix_ID) #TRUE

#Run correlation
protein_corrs <- list()

for (protein in protein_list) {
  print(protein)
  x <- episcore_v12 %>% select(all_of(protein)) %>% unlist()
  y <- scaled_pheno %>% select(all_of(protein)) %>% unlist()
  res <- cor.test(x,y)
  protein_corrs[[protein]] <- broom::tidy(res, conf.int = TRUE)
}

protein_corrdf <- bind_rows(protein_corrs, .id = "protein")
fwrite(protein_corrdf, file = "episcore_v12_proteincorr.csv")
head(protein_corrdf)
summary(protein_corrdf$estimate)

#identify significant correlations
sigcor <- protein_corrdf %>% filter(estimate > 0.1 & p.value < 0.05)
dim(sigcor) #109
sigcor %>% filter(protein %in% sig_epis) %>% dim()
fwrite(sigcor, file = "episcore_v12_proteincorr_sigonly.csv")
head(sigcor)

#Create df with both next to each other
head(protein_corrdf)
protein_corrdf <- protein_corrdf %>% mutate(type = "EPICv1&2")

origcorr <- fread(paste0(datadir, "correlation_true_epi_2407.csv"))
origcorr <- origcorr %>% mutate(type = "EPICv1")
head(origcorr)
combdf <- rbind(protein_corrdf, origcorr)
combdf <- combdf %>% filter(protein %in% sig_epis)

anno <- fread("protein_annos_forpub_nodup.csv")
anno <- anno %>% filter(protein %in% sig_epis)

#set up data for plotting
combdf <- left_join(combdf, anno, by = "protein")

combdf <- combdf %>%
  mutate(plot_label = ifelse(gene_names_simpl == "", protein, gene_names_simpl)) #TF is duplicated (serotransferrin and transferrin)
combdf <- combdf %>%
  mutate(plot_label = ifelse(plot_label == "TF" & protein == "C9JB55", protein, plot_label))
combdf <- combdf %>% mutate(
  plot_label = fct_reorder(plot_label, estimate, .desc = TRUE))
str(combdf$plot_label)
combdf <- combdf %>%
  mutate(plot_label = fct_reorder(plot_label, estimate * (type == "EPICv1"), .fun = mean, .desc = TRUE))


library(paletteer)
p <- combdf %>%
  ggplot(aes(x = plot_label, y = estimate, colour = type)) +
  geom_point(position = position_dodge(width = 1)) +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.8)) +
  scale_colour_manual(values = c("EPICv1" = "#003464FF",
                                "EPICv1&2"= "#7CB0C1FF")) +
                                guides(
    colour = guide_legend("CpG Loci")) + 
  theme_bw() +
  labs(x = "MS protein EpiScore", y = "Correlation with true protein level") +
  theme(axis.title.x = element_text(margin = margin(t = 30), size = 10),
        axis.title.y = element_text(margin = margin(r = 20), size = 10),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 8)) +
  coord_fixed(ratio = 60)

ggsave(p, file = "episcore_comp_proteincorr.pdf", width = 30, height = 15, units = "cm", dpi = 300)



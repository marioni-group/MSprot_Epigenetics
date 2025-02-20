#Bayesian results summary plots
#____________________________________________
library(ggridges)
library(paletteer)
library(tidyverse)
library(data.table)
library(gridExtra)
library(scales)
library(data.table)
library(patchwork)
library(cowplot)
library(egg)

#load annotated association file
highpip <- fread("/data_dir/ewas/bayesr/summary_data/cis_trans_anno_trial_1604.csv")

#summary plots: Generate counts of associations by protein & by CpG
Cpg_counts <- highpip %>% dplyr::count(protein, sort = TRUE)
protein_counts <- highpip %>% dplyr::count(CpGs, sort = TRUE)

#First plot: Number of CpG associations & number of proteins
t1 <- Cpg_counts %>%
  ggplot(aes(x = n)) +
  geom_histogram(binwidth = 1, fill = "#85A6F8", colour = "black") +
  scale_x_continuous(breaks = c(0:22)) +
  xlab("CpG association count") +
  ylab("Protein number") +
  theme_bw() +
  theme(axis.title.y = element_text(margin = margin(r = 5)),
        axis.title.x = element_text(margin = margin(t = 5)))

#Second plot: Number of protein associations & number of CpGs
t2 <- protein_counts %>%
  ggplot(aes(x = n)) +
  geom_histogram(binwidth = 1, fill = "#85A6F8", colour = "black") +
  scale_x_continuous(breaks = c(0:22)) +
  xlab("Protein association count") +
  ylab("CpG count") +
  theme_bw() +
  theme(axis.title.y = element_text(margin = margin(r = 5)),
        axis.title.x = element_text(margin = margin(t = 5)))

#Third plot: Heritability of proteins & number of associated CpGs
t3 <- heritability %>%
  drop_na() %>%
  ggplot(aes(x = mean_heritability, y = cpg_count)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = "Mean protein heritability", y = "No. CpGs with PIP > 0.95") +
  theme_bw() +
  theme(axis.title.y = element_text(margin = margin(r = 5)),
        axis.title.x = element_text(margin = margin(t = 5))) 
cor.test(heritability$mean_heritability, heritability$cpg_count)

#Fourth plot: Distribution of results compared with epic array
#load chisq df
chisq_df <- fread("/data_dir/ewas/bayesr/output/new_step2/df_for_chisq_plot.csv")

t4 <- chisq_df %>%
    ggplot(aes(x = origin, y = probes, fill = Relation_to_Island)) +
  geom_col(position = "fill", colour = "black") + 
  ylab("Proportion") + 
  xlab("Source") +
  guides(fill = guide_legend(title = "CpG Position")) +
  theme_bw() + 
  scale_fill_brewer(palette = "Blues") + 
   theme(legend.position = "right",
        axis.text.x = element_text(),
        axis.title.y = element_text(margin = margin(r = 5)),
        axis.title.x = element_text(margin = margin(t = 5))) +
  scale_x_discrete(labels=c('Epic Array', 'EWAS results'))

tspace <- ggplot() + theme_void()

#summary plots together
sumplot2 <- ggarrange(t1, tspace, t2, tspace, tspace, tspace, t3, tspace, t4, 
                      ncol = 3, nrow = 3, widths = c(2.5,0.3,2.5), 
                      heights = c(2, 0.3, 2.5),
                      labels = c("A", "", "B",
                                 "", "", "",
                                 "C", "", "D"))
ggsave("/data_dir/ewas/bayesr/sumplot2.png", sumplot2, height = 20, width = 25, units = "cm" )


#_______________________________________________________________________________________________________
#Effect-size plot

table_cis_trans <- all_anno %>% dplyr::count(protein, CIS_TRANS, sort = FALSE)
table_cis_trans$number <- table_cis_trans$n
table_cis_trans$CIS_TRANS[table_cis_trans$CIS_TRANS==""] <- NA

labels <- c("Cis", "Trans", "Unassigned")
p1 <- all_anno %>%
  select(effect_size, CIS_TRANS) %>%
  ggplot(aes(x = effect_size, y = CIS_TRANS)) +
  geom_density_ridges(aes(fill = CIS_TRANS), 
                      jittered_points = TRUE,
                      position = "raincloud",
                      alpha = 0.7,
                      scale = 0.55) +
  ylab("Association type") +
  xlab("Mean effect size") +
  scale_fill_brewer(palette = "Blues") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.title.x = element_text(margin = margin(t = 15))) +
  scale_y_discrete(label = labels)
#_______________________________________________________________________________________________________
#Cis_Trans plot

#load annotated association file
all_anno <- fread("/data_dir/ewas/bayesr/summary_data/cis_trans_anno_trial_1604.csv")

#Generate position of CpG Probe on Chromosomes
all_anno$probe_chr <- as.numeric(all_anno$probe_chr)
all_anno <- all_anno %>%  #gives position of probe on the chromosome
  mutate(., probe_pos_plot = as.character(pos/CHR_Length))
all_anno$probe_pos_plot_2 <- gsub("0[.]",".", all_anno$probe_pos_plot) #this removes the 0 from the start of each value in pos_plot 
all_anno$probe_pos_plot_2 = as.numeric(paste0(all_anno$probe_chr, all_anno$probe_pos_plot_2))

#Generate position of protein gene on chromosomes
all_anno$start <- as.numeric(all_anno$start)
all_anno$protein_pos_plot = as.character(as.numeric(all_anno$start/all_anno$protein_chr_length))
all_anno$protein_pos_plot <- gsub("0[.]",".", all_anno$protein_pos_plot)
all_anno$protein_pos_plot2 = as.numeric(paste0(all_anno$protein_chr, all_anno$protein_pos_plot))


p2 <- ggplot(all_anno, aes(probe_pos_plot_2, protein_pos_plot2)) +
ggplot2::geom_jitter(aes(colour = as.factor(CIS_TRANS)), size = 1) + 
xlab("CpG Position") +
ylab("Protein Position") +
scale_x_continuous(limits = c(1,23), 
                     breaks = c(1:23), 
                     labels = c(1:22,"X")) + 
scale_y_continuous(limits = c(1,23),
                    breaks = c(1:23),
                    labels = c(1:22,"X")) + 
theme_bw() +
theme(legend.title = element_blank(),
      legend.position = "top",
      axis.title.y = element_text(margin = margin(r = 15)),
      axis.title.x = element_text(margin = margin(t = 15))) + 
geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.1, color = "darkslategray4")


ewas_summary <- ggarrange(p1, ggplot() + theme_void(), p2, ncol = 3, widths = c(2.25, 1, 5), labels = c("E", "", "F")) 
ggsave("/data_dir/ewas/bayesr/cis_trans_effect_size.png", ewas_summary, height = 14, width = 25, units = "cm")   

#All plots together
bayes_summary2 <-  grid.arrange(sumplot2, tspace, ewas_summary,
                      ncol = 1, nrow = 3, 
                      heights = c(5, 0.3, 4)
                     )
ggsave("/data_dir/ewas/bayesr/bayes_results_summary2.png", bayes_summary2, height = 35, width = 25, units = "cm")
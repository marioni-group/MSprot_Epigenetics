#Episcore - CVD HR plot

library(tidyverse)
library(gridExtra)
library(scales)
library(data.table)
library(patchwork)
library(cowplot)
library(egg)
#____________________________________________
#HR forest plot for cox models (112 episcores)
#____________________________________________

#Load data: this df is already filtered to the 112 well performing episcores
combined <- fread("/data_dir/select_coxresults.csv")
df <- combined

#2. Create plot df
df <- combined %>% dplyr::filter(term == protein) #remove age and sex terms
df <- df %>% dplyr::group_by(protein) %>% 
  dplyr::filter(any(p.value < 0.05)) #filter to either predictor with nominal signif. 
df <- df %>% group_by(predictor) %>%
  mutate(protein = forcats::fct_reorder(protein, estimate)) %>%
  ungroup() #reorder so plot is easier to read
df <- df %>% 
  mutate(significant = ifelse(p.value < 0.05/112, "yes", "no"))  
head(df)  
df <- df %>%
  mutate(plot_lg = ifelse(predictor == "Episcore" & significant == "yes", "Significant_Episcore",
                          ifelse(predictor == "Episcore" & significant == "no", "Non-significant_Episcore",
                                 ifelse(predictor == "True" & significant == "yes", "Significant_protein", "Non-significant_protein"))))



paired <- RColorBrewer::brewer.pal(12, "Paired")
GnBu <- RColorBrewer::brewer.pal(9, "GnBu")


df <- df %>% dplyr::select(-std.error)
df$plot_lg <- factor(df$plot_lg, levels = c("Significant_Episcore", "Non-significant_Episcore", "Significant_protein", "Non-significant_protein"))
df <- df %>% mutate(plot_fill = plot_lg)

#Create plot HRs by predictor 
p1 <- ggplot(data=df, aes(y=protein, x=estimate, xmin=conf.low, xmax=conf.high, colour = plot_lg, shape = plot_lg)) +
  geom_point(position = position_dodge(width = 0.75), size = 2) + 
  geom_errorbarh(position = position_dodge(width = 0.75), height=.1) +
  scale_color_manual(values = c("Significant_Episcore" = paired[6],
                                "Non-significant_Episcore"= paired[5],
                                "Significant_protein"= GnBu[9],
                                "Non-significant_protein" = GnBu[6]),
                                labels = c("Significant EpiScore", "Non-significant EpiScore", "Significant protein", "Non-significant protein"))+
  scale_shape_manual(values = c("Significant_Episcore" = 19, 
                                "Non-significant_Episcore" = 1,
                                "Significant_protein" = 19, 
                                "Non-significant_protein" = 1)) +
  guides(
    colour = guide_legend("Model & Bonferroni Significance", 
                          override.aes = list(shape = c(19, 1))),  # Override shapes to match significance
    shape = "none") +
  labs(title='', x='Hazard Ratio [95% Confidence Interval]', y = 'Protein') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  coord_cartesian(xlim = c(0.5, 2)) +
  theme_minimal() + 
  theme(legend.position = c(.8, .2),
        legend.background = element_rect(fill="white", size=.4),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) 


ggsave("/datadir/HR_epi.pdf", p1, height = 40, width = 30, units = "cm")



#____________________________________________
#HR forest plot for nested cox models (significant on log likelihood test)
#____________________________________________

nested_models <- fread("/datadir/nested_model_comp_forplot.csv")  
head(nested_models)

#plot HRs by predictor 
p1 <- ggplot(data=nested_models, aes(y=protein, x=estimate, xmin=conf.low, xmax=conf.high, colour = predictor)) +
  geom_point(position = position_dodge(width = 0.75)) + 
  geom_errorbarh(position = position_dodge(width = 0.75), height=.1) +
  labs(title='HR for composite cvd outcome by Protein', x='HR', y = 'Protein') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  coord_cartesian(xlim = c(0.5, 2)) +
  theme_minimal()

ggsave(file = here("/nested_model_plot.svg"))

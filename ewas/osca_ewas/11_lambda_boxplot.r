#Boxplots of Lambda for OSCA models

library(tidyverse)
library(data.table)
library(stringr)

#Make boxplots of Lambda across 4 different models included in the osca EWAS

#Load summary results & clean
model1 <- fread("basic_byprot_df_133pbonf.csv")
model1 <- model1 %>% mutate(model = "model_1")

model2 <- fread("model2_byprot_summary_133pbonf.csv")
model2 <- model2 %>% 
    mutate(model = "model_2") %>%
    mutate(protein = str_remove_all(protein, "_wcc"))

model3 <- fread("model3summary_byprot_133pbonf.csv")
model3 <- model3 %>% 
    mutate(model = "model_3") %>%
    mutate(protein = str_remove_all(protein, "_smok"))

model4 <- fread("model4summary_byprot_133pbonf.csv")
model4 <- model4 %>% 
    mutate(model = "model_4")

#Combine for plotting
df <- rbind(model1, model2, model3, model4)

#plot
df %>%
ggplot(aes(x = model, y = Lambda)) +
 geom_boxplot() + 
 geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
 scale_y_continuous(breaks = c(0:16)) +
 theme_bw()
ggsave("lambda_boxplot.pdf")

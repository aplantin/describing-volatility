### This document compares distances/dissimilarities between consecutive    ###
### samples for a subject by time lag, study, and (for moving pictures) the ###
### specific distance or dissimilarity metric used.                          ### 

library(tidyverse)
library(gridExtra)
library(dplyr)
library(MBVolDescrip)

distvol_rarefy <- readRDS("VolSumms/distvol_rarefy.rds")
distvol_rarefy_mp <- readRDS("VolSumms/distvol_rarefy_mp.rds")


# Average dissimilarity by study, for different time lags 
distvol_rarefy$Study[distvol_rarefy$Study == "Student Microbiome Project (Gut)"] <- "SMP (Gut)"
p1 <- distvol_rarefy %>% 
  mutate(Study = factor(Study, 
                        levels = c("Moving Pictures (Gut)", "SMP (Gut)", 
                                   "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  mutate(Rarefaction = case_when(Rarefy == "None" ~ "Original", 
                                 Rarefy == "P100" ~ "100%", 
                                 Rarefy == "P80" ~ "80%", 
                                 Rarefy == "P60" ~ "60%")) %>% 
  mutate(Rarefaction = factor(Rarefaction, levels = c("Original", "100%", "80%", "60%"))) %>% 
  ggplot() + 
  geom_violin(aes(x=Rarefaction, y=BrayCurtis), width=1) + 
  geom_boxplot(aes(x=Rarefaction, y=BrayCurtis), outlier.size=0.5, width=0.1) + 
  ggtitle("Intra-Subject Bray-Curtis Dissimilarity By Rarefaction")  + 
  theme_bw() + 
  facet_wrap(vars(Study), nrow=1) + 
  theme(text=element_text(size=32)) + 
  xlab("Rarefaction") + ylab("Bray-Curtis Dissimilarity") 

png(filename = paste0(pre, "Figures/rarefy_braycurtis.png"),
    width = 1600, height = 500)
p1
dev.off() 


# Different dissimilarities & time lags 
p2 <- distvol_rarefy_mp %>% 
  pivot_longer(unweightedUF:BrayCurtis, names_to="Dissimilarity", values_to="Value") %>% 
  mutate(Distance = case_when(Dissimilarity == "unweightedUF" ~ "Unweighted UniFrac", 
                              Dissimilarity == "weightedUF" ~ "Weighted UniFrac", 
                              Dissimilarity == "genUF0.5" ~ "Gen. UniFrac (alpha=0.5)", 
                              Dissimilarity == "BrayCurtis" ~ "Bray-Curtis")) %>% 
  mutate(Distance = factor(Distance, levels = c("Unweighted UniFrac", "Gen. UniFrac (alpha=0.5)", 
                                                "Weighted UniFrac", "Bray-Curtis"))) %>% 
  mutate(Rarefaction = case_when(Rarefy == "None" ~ "Original", 
                                 Rarefy == "P100" ~ "100%", 
                                 Rarefy == "P80" ~ "80%", 
                                 Rarefy == "P60" ~ "60%")) %>% 
  mutate(Rarefaction = factor(Rarefaction, levels = c("Original", "100%", "80%", "60%"))) %>% 
  ggplot() + 
  geom_violin(aes(x=Rarefaction, y=Value), width=1) + 
  geom_boxplot(aes(x=Rarefaction, y=Value), outlier.size=0.5, width=0.1, color="grey", alpha=0.2) + 
  ggtitle("Moving Pictures Intra-Subject Dissimilarity By Distance Metric and Rarefaction")  + 
  theme_bw() + 
  facet_wrap(vars(Distance), nrow=1) + 
  theme(text=element_text(size=32)) + 
  xlab("Rarefaction") + ylab("Dissimilarity") 

png(filename = paste0(pre, "Figures/rarefy_distmetric_mp.png"),
    width = 1600, height = 500)
p2
dev.off() 


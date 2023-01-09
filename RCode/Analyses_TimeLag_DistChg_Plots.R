### This document compares distances/dissimilarities between consecutive    ###
### samples for a subject by time lag, study, and (for moving pictures) the ###
### specific distance or dissimilarity metric used.                          ### 

library(tidyverse)
library(gridExtra)
library(dplyr)
library(MBVolDescrip)

pre <- "~/AnnaDocuments/Research/Methods_indiv/2022_DescribingVolatility/"

mp_distvol <- readRDS(paste0(pre, "VolSumms/mp_distvol_rarefyNone.rds"))
smp_distvol <- readRDS(paste0(pre, "VolSumms/smp_distvol_rarefyNone.rds"))
gaj_distvol <- readRDS(paste0(pre, "VolSumms/gaj_distvol_rarefyNone.rds"))
rav_distvol <- readRDS(paste0(pre, "VolSumms/rav_distvol_rarefyNone.rds"))

# Merged dataset 
distvol_BC_all <- rbind(cbind(Study = "Moving Pictures (Gut)", 
                              mp_distvol[, c("subjID", "changeID", "GapSize", "BrayCurtis")]), 
                        cbind(Study = "Student Microbiome Project (Gut)", smp_distvol), 
                        cbind(Study = "Gajer (Vaginal)", gaj_distvol), 
                        cbind(Study = "Ravel (Vaginal)", rav_distvol)) %>% 
  mutate(TimeLag = paste0(GapSize, " day")) %>%
  mutate(TimeLag = factor(TimeLag, levels = c("1 day", "3 day", "7 day", "28 day")))  %>% 
  filter(!is.na(TimeLag))


# Average dissimilarity by study, for different time lags 
distvol_BC_all$Study[distvol_BC_all$Study == "Student Microbiome Project (Gut)"] <- "SMP (Gut)"
p1 <- distvol_BC_all %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", "SMP (Gut)", 
                                          "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  ggplot() + 
  geom_violin(aes(x=TimeLag, y=BrayCurtis), width=1) + 
  geom_boxplot(aes(x=TimeLag, y=BrayCurtis), outlier.size=0.5, width=0.1) + 
  ggtitle("Intra-Subject Bray-Curtis Dissimilarity By Time Lag")  + 
  theme_bw() + 
  facet_wrap(vars(Study), nrow=1) + 
  theme(text=element_text(size=32)) + 
  xlab("Time Lag") + ylab("Bray-Curtis Dissimilarity") 

png(filename = paste0(pre, "Figures/timelag_braycurtis.png"),
    width = 1600, height = 500)
p1
dev.off() 


# Different dissimilarities & time lags 
p2 <- mp_distvol %>% 
  pivot_longer(unweightedUF:BrayCurtis, names_to="Dissimilarity", values_to="Value") %>% 
  mutate(Distance = case_when(Dissimilarity == "unweightedUF" ~ "Unweighted UniFrac", 
                              Dissimilarity == "weightedUF" ~ "Weighted UniFrac", 
                              Dissimilarity == "genUF0.5" ~ "Gen. UniFrac (alpha=0.5)", 
                              Dissimilarity == "BrayCurtis" ~ "Bray-Curtis")) %>% 
  mutate(Distance = factor(Distance, levels = c("Unweighted UniFrac", "Gen. UniFrac (alpha=0.5)", 
                                                "Weighted UniFrac", "Bray-Curtis"))) %>% 
  mutate(TimeLag = paste0(GapSize, " day")) %>%
  mutate(TimeLag = factor(TimeLag, levels = c("1 day", "3 day", "7 day", "28 day")))  %>% 
  ggplot() + 
  geom_violin(aes(x=TimeLag, y=Value), width=1) + 
  geom_boxplot(aes(x=TimeLag, y=Value), outlier.size=0.5, width=0.1, color="grey", alpha=0.2) + 
  ggtitle("Moving Pictures Intra-Subject Dissimilarity By Distance Metric and Time Lag")  + 
  theme_bw() + 
  facet_wrap(vars(Distance), nrow=1) + 
  theme(text=element_text(size=32)) + 
  xlab("Time Lag") + ylab("Dissimilarity") 

png(filename = paste0(pre, "Figures/timelag_distmetric_mp.png"),
    width = 1600, height = 500)
p2
dev.off() 


### This document compares distances/dissimilarities between consecutive    ###
### samples for a subject by time lag, study, and (for moving pictures) the ###
### specific distance or dissimilarity metric used.                          ### 

library(tidyverse)
library(gridExtra)
library(dplyr)
library(MBVolDescrip)

pre <- "~/AnnaDocuments/Research/Methods_indiv/2022_DescribingVolatility/"


#source(paste0(pre, "Analyses_SamplingDensity_VolCalcs.R"))

# DistVol objects from VolCalcs 
head(mp_distvol)
head(smp_distvol)
head(gaj_distvol)
head(rav_distvol)

# Merged dataset 
distvol_BC_all <- rbind(cbind(Study = "Moving Pictures (Gut)", 
                              mp_distvol[, c("subjID", "changeID", "GapSize", "BrayCurtis")]), 
                        cbind(Study = "Student Microbiome Project (Gut)", smp_distvol), 
                        cbind(Study = "Gajer (Vaginal)", gaj_distvol), 
                        cbind(Study = "Ravel (Vaginal)", rav_distvol)) %>% 
  mutate(TimeLag = case_when(GapSize == 1 ~ "1 day", 
                             GapSize == 3 ~ "3 days", 
                             GapSize == 7 ~ "7 days", 
                             GapSize == 28 ~ "28 days")) %>%
  mutate(TimeLag = factor(TimeLag, levels = c("1 day", "3 days", "7 days", "28 days"))) 


# Average dissimilarity by study, for different time lags 
distvol_BC_all %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  ggplot() + 
  geom_violin(aes(x=TimeLag, y=BrayCurtis), width=1) + 
  geom_boxplot(aes(x=TimeLag, y=BrayCurtis), outlier.size=0.5, width=0.1) + 
  ggtitle("Intra-Subject Bray-Curtis Dissimilarity By Time Lag")  + 
  theme_bw() + 
  facet_wrap(vars(Study), nrow=1) + 
  theme(text=element_text(size=14)) + 
  xlab("Time Lag") + ylab("Bray-Curtis Dissimilarity") 


# Different dissimilarities & time lags 
mp_distvol %>% 
  pivot_longer(unweightedUF:BrayCurtis, names_to="Dissimilarity", values_to="Value") %>% 
  mutate(Distance = case_when(Dissimilarity == "unweightedUF" ~ "Unweighted UniFrac", 
                              Dissimilarity == "weightedUF" ~ "Weighted UniFrac", 
                              Dissimilarity == "genUF0.5" ~ "Generalized UniFrac (alpha=0.5)", 
                              Dissimilarity == "BrayCurtis" ~ "Bray-Curtis")) %>% 
  mutate(Distance = factor(Distance, levels = c("Unweighted UniFrac", "Generalized UniFrac (alpha=0.5)", 
                                                "Weighted UniFrac", "Bray-Curtis"))) %>% 
  mutate(TimeLag = case_when(GapSize == 1 ~ "1 day", 
                             GapSize == 3 ~ "3 days", 
                             GapSize == 7 ~ "7 days", 
                             GapSize == 28 ~ "28 days")) %>%
  mutate(TimeLag = factor(TimeLag, levels = c("1 day", "3 days", "7 days", "28 days")))  %>% 
  ggplot() + 
  geom_violin(aes(x=TimeLag, y=Value), width=1) + 
  geom_boxplot(aes(x=TimeLag, y=Value), outlier.size=0.5, width=0.1, color="grey", alpha=0.2) + 
  ggtitle("Moving Pictures Intra-Subject Dissimilarity By Distance Metric and Time Lag")  + 
  theme_bw() + 
  facet_wrap(vars(Distance), nrow=1) + 
  theme(text=element_text(size=14)) + 
  xlab("Time Lag") + ylab("Dissimilarity") 


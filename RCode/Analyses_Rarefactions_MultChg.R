### This document .....  ###

library(tidyverse)
library(gridExtra)
library(MBVolDescrip)

## Read in rarefications volatility data
vol_all_rarefy <- readRDS(file="VolSumms/vol_all_rarefy.rds") 

## Rarefaction box plot - avg additive change vs. rarefaction % 
vol_all_rarefy %>% 
  mutate(Study = factor(Study, 
                        levels = c("Moving Pictures (Gut)", "Student Microbiome Project (Gut)", 
                                   "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  mutate(Rarefaction = case_when(Rarefy == "None" ~ "Original", 
                                 Rarefy == "P100" ~ "100%", 
                                 Rarefy == "P80" ~ "80%", 
                                 Rarefy == "P60" ~ "60%")) %>% 
  mutate(Rarefaction = factor(Rarefaction, levels = c("Original", "100%", "80%", "60%"))) %>% 
  ggplot() + 
  geom_violin(aes(x=Rarefaction, y=(AvgNZLogFC_day7)), width=1) + 
  geom_boxplot(aes(x=Rarefaction, y=(AvgNZLogFC_day7)), 
               outlier.size=0.5, width=0.1) + 
  xlab("Rarefaction (Percent of Minimum Read Count)") + 
  ylab("Log Fold Change in Nonzero Abundance") + 
  facet_wrap(vars(Study), scales="free", nrow=1) + 
  theme_bw() + 
  theme(text=element_text(size=14), 
        legend.position = "bottom") + 
  ggtitle("Log Fold Change in Nonzero Abundance by Rarefaction (7-Day Lag)")


## Avg additive change vs. abundance, grouped by rarefaction % 
vol_all_rarefy %>% 
  mutate(Study = factor(Study, 
                        levels = c("Moving Pictures (Gut)", "Student Microbiome Project (Gut)", 
                                   "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  mutate(Rarefaction = case_when(Rarefy == "None" ~ "Original", 
                                 Rarefy == "P100" ~ "100%", 
                                 Rarefy == "P80" ~ "80%", 
                                 Rarefy == "P60" ~ "60%")) %>% 
  mutate(Rarefaction = factor(Rarefaction, levels = c("Original", "100%", "80%", "60%"))) %>% 
  ggplot() + 
  xlab("Log Average Relative Abundance") + 
  ylab("Log Fold Change in Nonzero Abundance") + 
  ylim(-5,5) + 
  geom_point(aes(x=log(AvgTaxAbund_day7), y=(AvgNZLogFC_day7),
                 group=Rarefy, color=Rarefy), size=0.1, alpha=0.15) + 
  geom_smooth(aes(x=log(AvgTaxAbund_day7), y=(AvgNZLogFC_day7),
                  group=Rarefy, color=Rarefy),  
              method = "loess", se=T) + 
  geom_abline(slope=0, intercept=0, color="black", lty=2) + 
  facet_wrap(vars(Study), scales="free", nrow=1) + 
  theme_bw() + 
  theme(text=element_text(size=14), 
        legend.position = "bottom") + 
  ggtitle("Log Fold Change in Nonzero Abundance by Rarefaction and Taxon Abundance")




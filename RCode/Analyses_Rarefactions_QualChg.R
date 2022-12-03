### This document extracts average qualitative change between time points ###
### and compares these averages by taxon abundance, time lag, and study   ###

library(tidyverse)
library(gridExtra)
library(dplyr)
library(MBVolDescrip)

## Read in rarefications volatility data
vol_all_rarefy <- readRDS(file="VolSumms/vol_all_rarefy.rds") 

# Average qualitative changes: all lags vs. taxon abundance 
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
  theme_bw() + 
  xlab("Log Average Relative Abundance") + 
  ylab("Proportion Qualitative Changes") + 
  geom_point(aes(x=log(AvgTaxAbund_day7), y=PropQualChange_day7,
                 group=Rarefaction, color=Rarefaction), size=0.1, alpha=0.15) + 
  geom_smooth(aes(x=log(AvgTaxAbund_day7), y=PropQualChange_day7,
                  group=Rarefaction, color=Rarefaction),  
              method = "loess", se=T) + 
  geom_abline(slope=0, intercept=0, color="black", lty=2) + 
  facet_wrap(vars(Study), nrow=1) + 
  theme(text=element_text(size=14), 
        legend.position = "bottom") + 
  ggtitle("Frequency of Presence/Absence Change by Rarefaction and Abundance")


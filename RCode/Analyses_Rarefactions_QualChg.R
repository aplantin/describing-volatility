### This document extracts average qualitative change between time points ###
### and compares these averages by taxon abundance, time lag, and study   ###

library(tidyverse)
library(gridExtra)
library(dplyr)

## Read in long-format rarefied volatility data
qual_longchg <- readRDS(paste0(pre, "VolSumms/vol_long_day7_allstudy_allrarefy.rds")) %>% 
  filter(AvgTaxAbund != 0) %>% 
  mutate(Study = factor(Study, 
                        levels = c("Moving Pictures (Gut)", "Student Microbiome Project (Gut)", 
                                   "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  mutate(Rarefaction = case_when(Rarefy == "None" ~ "Original", 
                                 Rarefy == "P100" ~ "100%", 
                                 Rarefy == "P80" ~ "80%", 
                                 Rarefy == "P60" ~ "60%")) %>% 
  mutate(Rarefaction = factor(Rarefaction, levels = c("Original", "100%", "80%", "60%"))) %>% 
  group_by(Study, taxID, AvgTaxAbund, Rarefaction) %>% 
  summarize(PropQualChange = mean(QualChg != 0)) 


# Make plot 
qual_longchg %>% 
  ggplot() + 
  theme_bw() + 
  ylim(0, 0.65) + 
  xlab("Log Average Relative Abundance") + 
  ylab("Proportion Qualitative Changes") + 
  geom_point(aes(x=log(AvgTaxAbund), y=PropQualChange,
                 group=Rarefaction, color=Rarefaction), size=0.1, alpha=0.15) + 
  geom_smooth(aes(x=log(AvgTaxAbund), y=PropQualChange,
                  group=Rarefaction, color=Rarefaction),  
              method = "loess", se=T) + 
  geom_abline(slope=0, intercept=0, color="black", lty=2) + 
  facet_wrap(vars(Study), nrow=1, scales="free_x") + 
  theme(text=element_text(size=14), 
        legend.position = "bottom") + 
  ggtitle("Proportion of Time Points with Presence/Absence Change")


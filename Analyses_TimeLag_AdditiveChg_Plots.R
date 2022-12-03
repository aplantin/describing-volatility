### This document extracts average additive change between time points  ###
### (not excluding pairs of observations with qualitative changes)      ###
### and compares these averages by taxon abundance, time lag, and study ###

library(tidyverse)
library(gridExtra)
library(MBVolDescrip)
pre <- "~/AnnaDocuments/Research/Methods_indiv/2022_DescribingVolatility/"

#source(paste0(pre, "Analyses_SamplingDensity_VolCalcs.R"))

# Moving Pictures 
mp_vol_all <- readRDS(paste0(pre, "VolSumms/mp_vol_taxlevel.rds")) 
mp_vol_long <- mp_vol_all %>% 
  select(taxID, 
         AvgTaxAbund = AvgTaxAbund_day1, 
         AvgAbsAdditive_day1, AvgAbsAdditive_day3, 
         AvgAbsAdditive_day7, AvgAbsAdditive_day28) %>% 
  pivot_longer(AvgAbsAdditive_day1:AvgAbsAdditive_day28, 
               names_to = "GapSize", 
               names_prefix = "AvgAbsAdditive_", 
               values_to = "AvgAbsAdditive")

# Student Microbiome Project
smp_vol_all <- readRDS(paste0(pre, "VolSumms/smp_vol_taxlevel.rds"))
smp_vol_long <- smp_vol_all %>% 
  select(taxID, 
         AvgTaxAbund = AvgTaxAbund_day7, 
         AvgAbsAdditive_day7, AvgAbsAdditive_day28) %>% 
  pivot_longer(AvgAbsAdditive_day7:AvgAbsAdditive_day28, 
               names_to = "GapSize", 
               names_prefix = "AvgAbsAdditive_", 
               values_to = "AvgAbsAdditive") %>% 
  filter(!is.na(AvgAbsAdditive))

# Gajer 
gaj_vol_all <- readRDS(paste0(pre, "VolSumms/gaj_vol_taxlevel.rds"))
gaj_vol_long <- gaj_vol_all %>% 
  select(taxID, 
         AvgTaxAbund = AvgTaxAbund_day3, 
         AvgAbsAdditive_day3, AvgAbsAdditive_day7, AvgAbsAdditive_day28) %>% 
  pivot_longer(AvgAbsAdditive_day3:AvgAbsAdditive_day28, 
               names_to = "GapSize", 
               names_prefix = "AvgAbsAdditive_", 
               values_to = "AvgAbsAdditive")

# Ravel 
rav_vol_all <- readRDS(paste0(pre, "VolSumms/rav_vol_taxlevel.rds"))
rav_vol_long <- rav_vol_all %>% 
  select(taxID, 
         AvgTaxAbund = AvgTaxAbund_day1, 
         AvgAbsAdditive_day1, AvgAbsAdditive_day3, AvgAbsAdditive_day7) %>% 
  pivot_longer(AvgAbsAdditive_day1:AvgAbsAdditive_day7, 
               names_to = "GapSize", 
               names_prefix = "AvgAbsAdditive_", 
               values_to = "AvgAbsAdditive")


# Merged dataset 
vol_all <- rbind(cbind(Study = "Moving Pictures (Gut)", mp_vol_long), 
                 cbind(Study = "Student Microbiome Project (Gut)", smp_vol_long), 
                 cbind(Study = "Gajer (Vaginal)", gaj_vol_long), 
                 cbind(Study = "Ravel (Vaginal)", rav_vol_long)) %>% 
  mutate(TimeLag = case_when(GapSize == "day1" ~ "1 day", 
                             GapSize == "day3" ~ "3 days", 
                             GapSize == "day7" ~ "7 days", 
                             GapSize == "day28" ~ "28 days")) %>%
  mutate(TimeLag = factor(TimeLag, levels = c("1 day", "3 days", "7 days", "28 days"))) %>% 
  filter(!is.na(AvgAbsAdditive))



## 
## Magnitude of Additive Changes at Different Time Scales (LOESS Curve)
## 

vol_all %>% 
  mutate(Study = factor(Study, 
                        levels = c("Moving Pictures (Gut)", "Student Microbiome Project (Gut)", 
                                   "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  ggplot() + 
  xlab("Log Average Relative Abundance") + 
  ylab("Log Average Absolute Additive Change") + 
  geom_point(aes(x=log(AvgTaxAbund), y=log(AvgAbsAdditive),
                 group=TimeLag, color=TimeLag), size=0.1, alpha=0.15) + 
  geom_smooth(aes(x=log(AvgTaxAbund), y=log(AvgAbsAdditive),
                  group=TimeLag, color=TimeLag),  
              method = "loess", se=T) + 
  geom_abline(slope=1, intercept=0, color="black", lty=2) + 
  facet_wrap(vars(Study), scales="free", nrow=1) + 
  theme_bw() + 
  theme(text=element_text(size=14), 
        legend.position = "bottom") + 
  ggtitle("Additive Changes in Abundance by Time Lag and Taxon Abundance")



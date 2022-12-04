### This document extracts average log fold change between time points   ###
### (excluding pairs of observations such that the taxon has zero        ### 
### abundance during at least one of the two, by necessity) and compares ###
### these averages by taxon abundance, time lag, and study               ###

library(tidyverse)
library(gridExtra)
library(dplyr)
pre <- "~/AnnaDocuments/Research/Methods_indiv/2022_DescribingVolatility/"

source(paste0(pre, "Analyses_TimeLag_VolCalcs.R"))

# Moving Pictures 
mp_vol_long <- mp_vol_all %>% 
  select(taxID, subjID, 
         AvgTaxAbund = AvgTaxAbund_day1, 
         AvgNZLogFC_day1, AvgNZLogFC_day3, 
         AvgNZLogFC_day7, AvgNZLogFC_day28) %>% 
  pivot_longer(AvgNZLogFC_day1:AvgNZLogFC_day28, 
               names_to = "GapSize", 
               names_prefix = "AvgNZLogFC_", 
               values_to = "AvgNZLogFC")

# Student Microbiome Project
smp_vol_long <- smp_vol_all %>% 
  select(taxID, subjID, 
         AvgTaxAbund = AvgTaxAbund_day7, 
         AvgNZLogFC_day7, AvgNZLogFC_day28) %>% 
  pivot_longer(AvgNZLogFC_day7:AvgNZLogFC_day28, 
               names_to = "GapSize", 
               names_prefix = "AvgNZLogFC_", 
               values_to = "AvgNZLogFC")

# Gajer 
gaj_vol_long <- gaj_vol_all %>% 
  select(taxID, subjID, 
         AvgTaxAbund = AvgTaxAbund_day3, 
         AvgNZLogFC_day3, AvgNZLogFC_day7, AvgNZLogFC_day28) %>% 
  pivot_longer(AvgNZLogFC_day3:AvgNZLogFC_day28, 
               names_to = "GapSize", 
               names_prefix = "AvgNZLogFC_", 
               values_to = "AvgNZLogFC")

# Ravel 
rav_vol_long <- rav_vol_all %>% 
  select(taxID, subjID, 
         AvgTaxAbund = AvgTaxAbund_day1, 
         AvgNZLogFC_day1, AvgNZLogFC_day3, AvgNZLogFC_day7) %>% 
  pivot_longer(AvgNZLogFC_day1:AvgNZLogFC_day7, 
               names_to = "GapSize", 
               names_prefix = "AvgNZLogFC_", 
               values_to = "AvgNZLogFC")


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
  filter(!is.na(AvgNZLogFC))


## Moving pictures, SMP, and Gajer: Day 28 vs. Day 7
## Moving Pictures, Gajer, and Ravel: Day 7 vs. Day 3

vol_wide <- vol_all %>% 
  pivot_wider(id_cols=Study:AvgTaxAbund, names_from=TimeLag, values_from = AvgNZLogFC) 

# Average log fold changes: 7 days vs. 28 days 
vol_wide %>% 
  filter(Study != "Ravel (Vaginal)") %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)"))) %>% 
  ggplot(aes(x=`7 days`, y = `28 days`)) + 
  geom_point(size=0.1) + 
  geom_abline(slope=1, intercept=0, lty=2) + 
  geom_smooth(aes(col="OLS"), method="lm", se=F) + 
  facet_wrap(vars(Study), scales="free") + 
  xlab("7 Day Lag") + ylab("28 Day Lag") + 
  ggtitle("Log Fold Change in Nonzero Abundance by Time Lag")  + 
  theme_bw() + 
  theme(legend.position="none", text=element_text(size=14)) 

# Average log fold changes: all lags vs. taxon abundance 
vol_all %>% 
  filter(Study != "Ravel (Vaginal)") %>% 
  mutate(Study = factor(Study, 
                        levels = c("Moving Pictures (Gut)", "Student Microbiome Project (Gut)", 
                                   "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  ggplot() + 
  ylim(-3,3) + 
  xlab("Log Average Relative Abundance") + 
  ylab("Log Fold Change in Nonzero Abundance") + 
  geom_point(aes(x=log(AvgTaxAbund), y=AvgNZLogFC,
                 group=TimeLag, color=TimeLag), size=0.1, alpha=0.15) + 
  geom_smooth(aes(x=log(AvgTaxAbund), y=AvgNZLogFC,
                  group=TimeLag, color=TimeLag),  
              method = "loess", se=T) + 
  geom_abline(slope=0, intercept=0, color="black", lty=2) + 
  facet_wrap(vars(Study), scales="free", nrow=1) + 
  theme_bw()  + 
  theme(text=element_text(size=14), 
        legend.position = "bottom") + 
  ggtitle("Log Fold Changes by Time Lag and Taxon Abundance")


# Table of SD of log fold changes 
## 201,026 nonzero pairs out of all 2,616,040 pairs of changes across the studies/time lags/etc. 
longfc_all <- rbind(cbind(Study = "Moving Pictures (Gut)", TimeLag = "1 day", mp_longchange_1),
                     cbind(Study = "Moving Pictures (Gut)", TimeLag = "3 days", mp_longchange_3), 
                     cbind(Study = "Moving Pictures (Gut)", TimeLag = "7 days", mp_longchange_7), 
                     cbind(Study = "Moving Pictures (Gut)", TimeLag = "28 days", mp_longchange_28), 
                     cbind(Study = "Student Microbiome Project (Gut)", TimeLag = "7 days", smp_longchange_7),  
                     cbind(Study = "Student Microbiome Project (Gut)", TimeLag = "28 days", smp_longchange_28),
                     cbind(Study = "Gajer (Vaginal)", TimeLag = "3 days", gaj_longchange_3),  
                     cbind(Study = "Gajer (Vaginal)", TimeLag = "7 days", gaj_longchange_7),  
                     cbind(Study = "Gajer (Vaginal)", TimeLag = "28 days", gaj_longchange_28),  
                     cbind(Study = "Ravel (Vaginal)", TimeLag = "1 day", rav_longchange_1),  
                     cbind(Study = "Ravel (Vaginal)", TimeLag = "3 days", rav_longchange_3),  
                     cbind(Study = "Ravel (Vaginal)", TimeLag = "7 days", rav_longchange_7)) %>% 
  mutate(TimeLag = factor(TimeLag, levels = c("1 day", "3 days", "7 days", "28 days"))) %>% 
  select(Study, TimeLag, taxID, changeID, subjID, NZLogFC, NZPairInd, AvgTaxAbund, PropTaxNZ) %>% 
  filter(NZPairInd == 1) 
longfc_all %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  group_by(Study, TimeLag) %>% 
  summarize(sdLogFC = sd(NZLogFC), 
            nPairs = n()) %>% 
  pivot_wider(names_from=TimeLag, values_from = sdLogFC) %>% 
  select(Study, `1 day`, `3 days`, `7 days`, `28 days`) 

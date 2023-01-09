library(tidyverse)
library(gridExtra)
library(dplyr)
pre <- "~/AnnaDocuments/Research/Methods_indiv/2022_DescribingVolatility/"


# Long-format changes --> overall SDs 
longchg_overall <- readRDS(paste0(pre, "VolSumms/vol_long_day7_allstudy_allrarefy.rds")) %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", "Ravel (Vaginal)")), 
         Rarefy = factor(Rarefy, levels = c("None", "P100", "P80", "P60"))) %>% 
  select(Study, Rarefy, taxID, changeID, subjID, AdditiveChg, NZLogFC, CLRChg, AvgTaxAbund) %>% 
  mutate(TaxAbundCat = cut(AvgTaxAbund, c(0, 1e-5, 1e-4, 1e-3, 1), include.lowest=T)) %>% 
  group_by(Study, Rarefy) %>% 
  summarize(sdAddChg = sd(AdditiveChg), 
            nAddChg = sum(!is.na(AdditiveChg)), 
            sdLogFC = sd(NZLogFC, na.rm=T), 
            nLogFC = sum(!is.na(NZLogFC)), 
            sdCLR = sd(CLRChg), 
            nCLR = sum(!is.na(CLRChg))) %>% 
  mutate(TaxAbundCat = "All", 
         AddChg = paste0(round(sdAddChg, 6), " (", nAddChg, ")"), 
         LogFC = paste0(round(sdLogFC, 2), " (", nLogFC, ")"), 
         CLRChg = paste0(round(sdCLR, 2), " (", nCLR, ")")) %>% 
  select(Study, Rarefy, TaxAbundCat, AddChg, LogFC, CLRChg) %>% 
  pivot_longer(AddChg:CLRChg, names_to="ChangeMeasure", values_to="ChangeSD") %>% 
  mutate(MeasureOfChange = case_when(ChangeMeasure == "AddChg" ~ "Additive", 
                                     ChangeMeasure == "LogFC" ~ "Log Fold Change", 
                                     ChangeMeasure == "CLRChg" ~ "CLR-Based")) %>% 
  pivot_wider(names_from=Rarefy, values_from=ChangeSD) 


# Long-format changes --> taxon abundance specific SDs 
longchg_abund <- readRDS(paste0(pre, "VolSumms/vol_long_day7_allstudy_allrarefy.rds")) %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", "Ravel (Vaginal)")), 
         Rarefy = factor(Rarefy, levels = c("None", "P100", "P80", "P60"))) %>% 
  select(Study, Rarefy, taxID, changeID, subjID, AdditiveChg, NZLogFC, CLRChg, AvgTaxAbund) %>% 
  mutate(TaxAbundCat = cut(AvgTaxAbund, c(0, 1e-5, 1e-4, 1e-3, 1), include.lowest=T)) %>% 
  group_by(Study, Rarefy, TaxAbundCat) %>% 
  summarize(sdAddChg = sd(AdditiveChg), 
            nAddChg = sum(!is.na(AdditiveChg)), 
            sdLogFC = sd(NZLogFC, na.rm=T), 
            nLogFC = sum(!is.na(NZLogFC)), 
            sdCLR = sd(CLRChg), 
            nCLR = sum(!is.na(CLRChg))) %>% 
  mutate(AddChg = paste0(round(sdAddChg, 6), " (", nAddChg, ")"), 
         LogFC = paste0(round(sdLogFC, 2), " (", nLogFC, ")"), 
         CLRChg = paste0(round(sdCLR, 2), " (", nCLR, ")")) %>% 
  select(Study, Rarefy, TaxAbundCat, AddChg, LogFC, CLRChg) %>% 
  pivot_longer(AddChg:CLRChg, names_to="ChangeMeasure", values_to="ChangeSD") %>% 
  mutate(MeasureOfChange = case_when(ChangeMeasure == "AddChg" ~ "Additive", 
                                     ChangeMeasure == "LogFC" ~ "Log Fold Change", 
                                     ChangeMeasure == "CLRChg" ~ "CLR-Based")) %>% 
  pivot_wider(names_from=Rarefy, values_from=ChangeSD) 

# Merge and sort 
longchg_all <- rbind(longchg_overall, longchg_abund) %>% 
  mutate(TaxAbundCat = factor(TaxAbundCat, levels = c("All", "[0,1e-05]", 
                                                      "(1e-05,0.0001]", "(0.0001,0.001]",
                                                      "(0.001,1]"))) %>% 
  select(-ChangeMeasure) %>% 
  arrange(MeasureOfChange, Study, TaxAbundCat)

# Write to file 
write.csv(longchg_all, file = paste0(pre, "Tables/rarefy_sds.csv"), row.names = F)

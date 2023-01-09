
library(tidyverse)
library(gridExtra)
library(dplyr)
pre <- "~/AnnaDocuments/Research/Methods_indiv/2022_DescribingVolatility/"

# Long-format additive changes 
mp_longchange <- readRDS(paste0(pre, "VolSumms/mp_longchg_rarefyP100.rds")) %>% 
  do.call(rbind, .)
smp_longchange <- readRDS(paste0(pre, "VolSumms/smp_longchg_rarefyP100.rds")) %>% 
  do.call(rbind, .) 
gaj_longchange <- readRDS(paste0(pre, "VolSumms/gaj_longchg_rarefyP100.rds")) %>% 
  do.call(rbind, .) 
rav_longchange <- readRDS(paste0(pre, "VolSumms/rav_longchg_rarefyP100.rds")) %>% 
  do.call(rbind, .)



longfc_all <- rbind(cbind(Study = "Moving Pictures (Gut)", mp_longchange),
                    cbind(Study = "Student Microbiome Project (Gut)", smp_longchange),  
                    cbind(Study = "Gajer (Vaginal)", gaj_longchange),  
                    cbind(Study = "Ravel (Vaginal)", rav_longchange))%>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", "Ravel (Vaginal)")), 
         TimeLag = factor(paste0(NominalLag, " day"), 
                          levels = c("1 day", "3 day", "7 day", "28 day")), 
         TaxAbundCat = cut(AvgTaxAbund, c(0, 1e-5, 1e-4, 1e-3, 1), include.lowest=T)) %>% 
  group_by(Study, TimeLag, TaxAbundCat) %>% 
  summarize(sdAddChg = sd(AdditiveChg), 
            nAddChg = sum(!is.na(AdditiveChg)), 
            sdLogFC = sd(NZLogFC, na.rm=T), 
            nLogFC = sum(!is.na(NZLogFC)), 
            sdCLR = sd(CLRChg), 
            nCLR = sum(!is.na(CLRChg))) %>% 
  mutate(AddChg = paste0(round(sdAddChg, 6), " (", nAddChg, ")"), 
         LogFC = paste0(round(sdLogFC, 2), " (", nLogFC, ")"), 
         CLRChg = paste0(round(sdCLR, 2), " (", nCLR, ")")) %>% 
  select(Study, TimeLag, TaxAbundCat, AddChg, LogFC, CLRChg) %>% 
  pivot_longer(AddChg:CLRChg, names_to="ChangeMeasure", values_to="ChangeSD") %>% 
  mutate(MeasureOfChange = case_when(ChangeMeasure == "AddChg" ~ "Additive", 
                                     ChangeMeasure == "LogFC" ~ "Log Fold Change", 
                                     ChangeMeasure == "CLRChg" ~ "CLR-Based")) %>% 
  pivot_wider(names_from=TimeLag, values_from=ChangeSD) 


# Long-format changes --> overall SDs 
longchg_overall <- rbind(cbind(Study = "Moving Pictures (Gut)", mp_longchange),
                         cbind(Study = "Student Microbiome Project (Gut)", smp_longchange),  
                         cbind(Study = "Gajer (Vaginal)", gaj_longchange),  
                         cbind(Study = "Ravel (Vaginal)", rav_longchange)) %>% 
  mutate(TimeLag = factor(paste0(NominalLag, " day"), 
                          levels = c("1 day", "3 day", "7 day", "28 day"))) %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  group_by(Study, TimeLag) %>% 
  summarize(sdAddChg = sd(AdditiveChg), 
            nAddChg = sum(!is.na(AdditiveChg)), 
            sdLogFC = sd(NZLogFC, na.rm=T), 
            nLogFC = sum(!is.na(NZLogFC)), 
            sdCLR = sd(CLRChg), 
            nCLR = sum(!is.na(CLRChg))) %>% 
  mutate(AddChg = paste0(round(sdAddChg, 6), " (", nAddChg, ")"), 
         LogFC = paste0(round(sdLogFC, 2), " (", nLogFC, ")"), 
         CLRChg = paste0(round(sdCLR, 2), " (", nCLR, ")")) %>% 
  mutate(TaxAbundCat = "All") %>% 
  select(Study, TimeLag, TaxAbundCat, AddChg, LogFC, CLRChg) %>% 
  pivot_longer(AddChg:CLRChg, names_to="ChangeMeasure", values_to="ChangeSD") %>% 
  mutate(MeasureOfChange = case_when(ChangeMeasure == "AddChg" ~ "Additive", 
                                     ChangeMeasure == "LogFC" ~ "Log Fold Change", 
                                     ChangeMeasure == "CLRChg" ~ "CLR-Based")) %>% 
  pivot_wider(names_from=TimeLag, values_from=ChangeSD) 
  
 
# Long-format changes --> taxon abundance specific SDs 
longchg_abund <- rbind(cbind(Study = "Moving Pictures (Gut)", mp_longchange),
                         cbind(Study = "Student Microbiome Project (Gut)", smp_longchange),  
                         cbind(Study = "Gajer (Vaginal)", gaj_longchange),  
                         cbind(Study = "Ravel (Vaginal)", rav_longchange))%>% 
  mutate(TimeLag = factor(paste0(NominalLag, " day"), 
                          levels = c("1 day", "3 day", "7 day", "28 day"))) %>% 
  select(Study, TimeLag, taxID, changeID, subjID, AdditiveChg, NZLogFC, CLRChg, AvgTaxAbund) %>% 
  mutate(TaxAbundCat = cut(AvgTaxAbund, c(0, 1e-5, 1e-4, 1e-3, 1), include.lowest=T)) %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  group_by(Study, TimeLag, TaxAbundCat) %>% 
  summarize(sdAddChg = sd(AdditiveChg), 
            nAddChg = sum(!is.na(AdditiveChg)), 
            sdLogFC = sd(NZLogFC, na.rm=T), 
            nLogFC = sum(!is.na(NZLogFC)), 
            sdCLR = sd(CLRChg), 
            nCLR = sum(!is.na(CLRChg))) %>% 
  mutate(AddChg = paste0(round(sdAddChg, 6), " (", nAddChg, ")"), 
         LogFC = paste0(round(sdLogFC, 2), " (", nLogFC, ")"), 
         CLRChg = paste0(round(sdCLR, 2), " (", nCLR, ")")) %>% 
  select(Study, TimeLag, TaxAbundCat, AddChg, LogFC, CLRChg) %>% 
  pivot_longer(AddChg:CLRChg, names_to="ChangeMeasure", values_to="ChangeSD") %>% 
  mutate(MeasureOfChange = case_when(ChangeMeasure == "AddChg" ~ "Additive", 
                                     ChangeMeasure == "LogFC" ~ "Log Fold Change", 
                                     ChangeMeasure == "CLRChg" ~ "CLR-Based")) %>% 
  pivot_wider(names_from=TimeLag, values_from=ChangeSD) 


# Merge and sort 
longchg_all <- rbind(longchg_overall, longchg_abund) %>% 
  mutate(TaxAbundCat = factor(TaxAbundCat, levels = c("All", "[0,1e-05]", 
                                                      "(1e-05,0.0001]", "(0.0001,0.001]",
                                                      "(0.001,1]"))) %>% 
  select(-ChangeMeasure) %>% 
  arrange(MeasureOfChange, Study, TaxAbundCat)

# Write to file 
write.csv(longchg_all, file = paste0(pre, "Tables/timelag_sds.csv"), row.names = F)

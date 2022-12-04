### This document .....  ###

library(tidyverse)
library(vegan)
library(GUniFrac)
library(MBVolDescrip)

## days 1, 3, 7, 28: 3rd entry is day 7 for MP  
mp_longchange_none <- readRDS(paste0("VolSumms/mp_longchg_rarefyNone.rds"))[[3]]
mp_longchange_p100 <- readRDS(paste0("VolSumms/mp_longchg_rarefyP100.rds"))[[3]]
mp_longchange_p80 <- readRDS(paste0("VolSumms/mp_longchg_rarefyP80.rds"))[[3]]
mp_longchange_p60 <- readRDS(paste0("VolSumms/mp_longchg_rarefyP60.rds"))[[3]]

## day 7, day 28 (so first entry is day 7)
smp_longchange_none <- readRDS(paste0("VolSumms/smp_longchg_rarefyNone.rds"))[[1]]
smp_longchange_p100 <- readRDS(paste0("VolSumms/smp_longchg_rarefyP100.rds"))[[1]]
smp_longchange_p80 <- readRDS(paste0("VolSumms/smp_longchg_rarefyP80.rds"))[[1]]
smp_longchange_p60 <- readRDS(paste0("VolSumms/smp_longchg_rarefyP60.rds"))[[1]]

## days 3, 7, 28 (2nd entry) 
gaj_longchange_none <- readRDS(paste0("VolSumms/gaj_longchg_rarefyNone.rds"))[[2]]
gaj_longchange_p100 <- readRDS(paste0("VolSumms/gaj_longchg_rarefyP100.rds"))[[2]]
gaj_longchange_p80 <- readRDS(paste0("VolSumms/gaj_longchg_rarefyP80.rds"))[[2]]
gaj_longchange_p60 <- readRDS(paste0("VolSumms/gaj_longchg_rarefyP60.rds"))[[2]]

## days 1, 3, 7 (3rd entry) 
rav_longchange_none <- readRDS(paste0("VolSumms/rav_longchg_rarefyNone.rds"))[[3]]
rav_longchange_p100 <- readRDS(paste0("VolSumms/rav_longchg_rarefyP100.rds"))[[3]]
rav_longchange_p80 <- readRDS(paste0("VolSumms/rav_longchg_rarefyP80.rds"))[[3]]
rav_longchange_p60 <- readRDS(paste0("VolSumms/rav_longchg_rarefyP60.rds"))[[3]]

## combine 
vol_long_rarefy <- rbind(cbind(Study = "Moving Pictures (Gut)", Rarefy = "None", mp_longchange_none), 
                         cbind(Study = "Moving Pictures (Gut)", Rarefy = "P100", mp_longchange_p100), 
                          cbind(Study = "Moving Pictures (Gut)", Rarefy = "P80", mp_longchange_p80), 
                          cbind(Study = "Moving Pictures (Gut)", Rarefy = "P60", mp_longchange_p60), 
                          cbind(Study = "Student Microbiome Project (Gut)", Rarefy = "None", smp_longchange_none),
                          cbind(Study = "Student Microbiome Project (Gut)", Rarefy = "P100", smp_longchange_p100), 
                          cbind(Study = "Student Microbiome Project (Gut)", Rarefy = "P80", smp_longchange_p80), 
                          cbind(Study = "Student Microbiome Project (Gut)", Rarefy = "P60", smp_longchange_p60), 
                          cbind(Study = "Gajer (Vaginal)", Rarefy = "None", gaj_longchange_none), 
                          cbind(Study = "Gajer (Vaginal)", Rarefy = "P100", gaj_longchange_p100), 
                          cbind(Study = "Gajer (Vaginal)", Rarefy = "P80", gaj_longchange_p80), 
                          cbind(Study = "Gajer (Vaginal)", Rarefy = "P60", gaj_longchange_p60), 
                          cbind(Study = "Ravel (Vaginal)", Rarefy = "None", rav_longchange_none), 
                          cbind(Study = "Ravel (Vaginal)", Rarefy = "P100", rav_longchange_p100), 
                          cbind(Study = "Ravel (Vaginal)", Rarefy = "P80", rav_longchange_p80), 
                          cbind(Study = "Ravel (Vaginal)", Rarefy = "P60", rav_longchange_p60))  %>% 
  filter(!is.na(NZLogFC))

vol_long_rarefy %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  mutate(Rarefy = factor(Rarefy, levels = c("None", "P100", "P80", "P60"))) %>% 
  group_by(Study, Rarefy) %>%
  summarize(sdLogFC = sd(NZLogFC), 
            NinGrp = n()) %>% 
  View() 
  
  
  

library(tidyverse)
library(gridExtra)
library(dplyr)
library(ggallin) # for pseudolog10 scale 
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
  mutate(TimeLag = factor(paste0(NominalLag, " day"), 
                          levels = c("1 day", "3 day", "7 day", "28 day"))) %>% 
  select(Study, TimeLag, taxID, changeID, subjID, AdditiveChg, NZLogFC, CLRChg, AvgTaxAbund) %>% 
  mutate(TaxAbundCat = cut(AvgTaxAbund, c(0, 1e-5, 1e-4, 1e-3, 1), include.lowest=T)) %>% 
  group_by(Study, TimeLag, TaxAbundCat) %>% 
  pivot_longer(c(AdditiveChg, NZLogFC, CLRChg), names_to="ChangeMeasure", values_to="ChangeVal") %>% 
  mutate(MeasureOfChange = case_when(ChangeMeasure == "AdditiveChg" ~ "Additive", 
                                     ChangeMeasure == "NZLogFC" ~ "Log Fold Change", 
                                     ChangeMeasure == "CLRChg" ~ "CLR-Based"))
  
  
longfc_all %>% 
    mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                            "Student Microbiome Project (Gut)", 
                                            "Gajer (Vaginal)", 
                                            "Ravel (Vaginal)"))) %>% 
  ggplot(aes(x=TaxAbundCat, y=ChangeVal, fill=TimeLag)) + 
  geom_boxplot() + 
  facet_grid(cols=vars(Study), rows=vars(MeasureOfChange), scales="free_y") + 
  scale_x_discrete(guide = guide_axis(n.dodge=2)) + 
  scale_y_continuous(trans=pseudolog10_trans) + 
  scale_fill_grey() + 
  theme_bw() +
  theme(text=element_text(size=12)) + 
  xlab("Taxon Abundance Category") + 
  ylab("Change Measure")

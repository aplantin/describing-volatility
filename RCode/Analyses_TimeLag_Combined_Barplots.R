
library(tidyverse)
library(gridExtra)
library(dplyr)
pre <- "~/AnnaDocuments/Research/Methods_indiv/2022_DescribingVolatility/"

# Long-format changes 
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
  summarize(sdAddChg = sd(AdditiveChg), 
            sdLogFC = sd(NZLogFC, na.rm=T), 
            sdCLR = sd(CLRChg)) %>% 
  pivot_longer(sdAddChg:sdCLR, names_to="ChangeMeasure", values_to="ChangeSD") %>% 
  mutate(MeasureOfChange = case_when(ChangeMeasure == "sdAddChg" ~ "Additive", 
                                     ChangeMeasure == "sdLogFC" ~ "Log Fold Change", 
                                     ChangeMeasure == "sdCLR" ~ "CLR-Based"))
  
  
p <- longfc_all %>% 
    mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                            "Student Microbiome Project (Gut)", 
                                            "Gajer (Vaginal)", 
                                            "Ravel (Vaginal)"))) %>% 
  mutate(MeasureOfChange = factor(MeasureOfChange, levels = c("Additive", "Log Fold Change", "CLR-Based"))) %>% 
  ggplot(aes(x=TaxAbundCat, y=ChangeSD, fill=TimeLag)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_grid(cols=vars(Study), rows=vars(MeasureOfChange), scales="free_y") + 
  scale_x_discrete(guide = guide_axis(n.dodge=2)) + 
  scale_fill_grey() + 
  theme_bw() + 
  theme(text=element_text(size=24)) + 
  #geom_text(aes(label=nPairs, group=TimeLag), vjust=0.5, hjust=0, size=3, 
  #          position=position_dodge(width=0.9), angle=90) + 
  xlab("Taxon Abundance Category") + 
  ylab("Standard Deviation of Change Measure")


png(filename = paste0(pre, "Figures/combined_barplot_timelag.png"), width = 1600, height = 1200)
p 
dev.off() 


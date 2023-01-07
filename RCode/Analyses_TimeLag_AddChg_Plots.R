### This document extracts average additive change between time points   ###
### and compares these averages by taxon abundance, time lag, and study  ###

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
  mutate(TimeLag = factor(paste0(NominalLag, " day"), 
                          levels = c("1 day", "3 day", "7 day", "28 day"))) %>% 
  select(Study, TimeLag, taxID, changeID, subjID, AdditiveChg, AvgTaxAbund) 



# Bar graph of SD of additive changes, separated by taxon abundance 
longfc_all %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", 
                                          "Ravel (Vaginal)"))) %>% 
  mutate(TaxAbundCat = cut(AvgTaxAbund, c(0, 1e-5, 1e-4, 1e-3, 1), include.lowest=T)) %>% 
  group_by(Study, TimeLag, TaxAbundCat) %>% 
  summarize(sdAddChg = sd(AdditiveChg), 
            nPairs = n()) %>% 
  ggplot(aes(x=TaxAbundCat, y=sdAddChg, fill=TimeLag)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(vars(Study), nrow=1) + 
  scale_x_discrete(guide = guide_axis(n.dodge=2)) + 
  scale_fill_grey() + 
  theme_bw() + 
  #geom_text(aes(label=nPairs, group=TimeLag), vjust=0.5, hjust=0, size=3, 
  #          position=position_dodge(width=0.9), angle=90) + 
  xlab("Taxon Abundance Category") + 
  ylab("Standard Deviation of Additive Change")


# Scatterplot with loess showing changes are centered around 0 
longfc_avgs <- longfc_all %>% 
  select(Study, TimeLag, taxID, AdditiveChg, AvgTaxAbund) %>% 
  group_by(Study, TimeLag, taxID, AvgTaxAbund) %>% 
  summarize(meanChg = mean(AdditiveChg))
longfc_avgs %>% 
  mutate(Study = factor(Study, 
                        levels = c("Moving Pictures (Gut)", "Student Microbiome Project (Gut)", 
                                   "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  ggplot() + 
  xlab("Log Average Relative Abundance") + 
  ylab("Change in Abundance") + 
  geom_point(aes(x=log(AvgTaxAbund), y=meanChg,
                 group=TimeLag, color=TimeLag), size=0.1, alpha=0.15) + 
  geom_smooth(aes(x=log(AvgTaxAbund), y=meanChg,
                  group=TimeLag, color=TimeLag),  
              method = "loess", se=T) + 
  geom_abline(slope=0, intercept=0, color="black", lty=2) + 
  facet_wrap(vars(Study), nrow=2) + 
  theme_bw()  + 
  theme(text=element_text(size=14), 
        legend.position = "bottom") + 
  ggtitle("Additive Changes by Time Lag and Taxon Abundance")



# Table of SD of additive changes, separated by taxon abundance & overall 
sd.all <- longfc_all %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", 
                                          "Ravel (Vaginal)"))) %>% 
  group_by(Study, TimeLag) %>% 
  summarize(sdAddChg = sd(AdditiveChg), 
            nPairs = n()) %>% 
  mutate(TaxAbundCat = "All") %>% 
  select(Study, TimeLag, TaxAbundCat, sdAddChg, nPairs)
sd.abund <- longfc_all %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", 
                                          "Ravel (Vaginal)"))) %>% 
  mutate(TaxAbundCat = cut(AvgTaxAbund, c(0, 1e-5,  1e-3, 1), include.lowest=T)) %>% 
  group_by(Study, TimeLag, TaxAbundCat) %>% 
  summarize(sdAddChg = sd(AdditiveChg), 
            nPairs = n()) 

sd.add <- rbind(sd.all, sd.abund) %>% 
  mutate(TaxAbundCat = factor(TaxAbundCat, levels = c("All", "[0,1e-05]", "(1e-05,0.0001]", 
                                                      "(0.0001,0.001]", "(0.001,1]")), 
         SD.n = paste0(round(sdAddChg, 8), " (", nPairs, ")")) %>% 
  select(Study, TimeLag, TaxAbundCat, SD.n) %>% 
  pivot_wider(names_from=TimeLag, values_from=SD.n) %>% 
  arrange(Study, TaxAbundCat)
write.csv(sd.add, file=paste0(pre, "Tables/sd_additive.csv"), row.names=F)

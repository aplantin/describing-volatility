### This document extracts average qualitative change between time points ###
### and compares these averages by taxon abundance, time lag, and study   ###
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
                          levels = c("1 day", "3 day", "7 day", "28 day")))

# Qualitative volatility 
qualvol <- longfc_all %>% 
  filter(AvgTaxAbund != 0) %>% 
  select(Study, taxID, AvgTaxAbund, TimeLag, QualChg) %>% 
  group_by(Study, taxID, AvgTaxAbund, TimeLag) %>% 
  summarize(PropQualChg = mean(QualChg != 0)) 


# Average qualitative changes: all lags vs. taxon abundance 
p <- qualvol %>% 
  mutate(Study = factor(Study, 
                        levels = c("Moving Pictures (Gut)", "Student Microbiome Project (Gut)", 
                                   "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  ggplot(aes(x=log(AvgTaxAbund), y=PropQualChg, group=TimeLag, color=TimeLag)) + 
  theme_bw() + 
  xlab("Log Average Relative Abundance") + 
  ylab("Proportion Qualitative Changes") + 
  ylim(0, 0.75) + 
  geom_point(size=0.1, alpha=0.15) + 
  geom_smooth(method = "loess", se=T) + 
  geom_abline(slope=0, intercept=0, color="black", lty=2) + 
  facet_wrap(vars(Study), nrow=1) + 
  theme(text=element_text(size=24), 
        legend.position = "bottom") + 
  ggtitle("Proportion of Time Points with Presence/Absence Change")



png(filename = paste0(pre, "Figures/timelag_qualchange_vs_abundance.png"),
    width = 1600, height = 500)
p 
dev.off() 


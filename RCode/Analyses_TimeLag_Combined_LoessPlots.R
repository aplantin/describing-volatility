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
  select(Study, TimeLag, taxID, changeID, subjID, AdditiveChg, NZLogFC, CLRChg, AvgTaxAbund) 


# Scatterplot with loess showing changes are centered around 0 
longfc_avgs <- longfc_all %>% 
  select(Study, TimeLag, taxID, AdditiveChg, NZLogFC, CLRChg, AvgTaxAbund) %>% 
  group_by(Study, TimeLag, taxID, AvgTaxAbund) %>% 
  summarize(meanAddChg = mean(AdditiveChg), 
            meanLogFC = mean(NZLogFC, na.rm=T), 
            meanCLRChg = mean(CLRChg)) %>% 
  pivot_longer(meanAddChg:meanCLRChg, names_to="ChangeMeasure", values_to="ChangeMean") %>%
  mutate(MeasureOfChange = case_when(ChangeMeasure == "meanAddChg" ~ "Additive", 
                                     ChangeMeasure == "meanLogFC" ~ "Log Fold Change", 
                                     ChangeMeasure == "meanCLRChg" ~ "CLR-Based"))

longfc_avgs %>% 
  mutate(Study = factor(Study, 
                        levels = c("Moving Pictures (Gut)", "Student Microbiome Project (Gut)", 
                                   "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  ggplot() + 
  xlab("Log Average Relative Abundance") + 
  ylab("Change in Abundance") + 
  geom_point(aes(x=log(AvgTaxAbund), y=ChangeMean,
                 group=TimeLag, color=TimeLag), size=0.1, alpha=0.15) + 
  geom_smooth(aes(x=log(AvgTaxAbund), y=ChangeMean,
                  group=TimeLag, color=TimeLag),  
              method = "loess", se=T) + 
  geom_abline(slope=0, intercept=0, color="black", lty=2) + 
  facet_grid(rows=vars(MeasureOfChange), cols=vars(Study), scales="free") + 
  theme_bw()  + 
  theme(text=element_text(size=14), 
        legend.position = "bottom")


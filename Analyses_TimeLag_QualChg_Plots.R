### This document extracts average qualitative change between time points ###
### and compares these averages by taxon abundance, time lag, and study   ###

library(tidyverse)
library(gridExtra)
library(dplyr)
library(MBVolDescrip)

pre <- "~/AnnaDocuments/Research/Methods_indiv/2022_DescribingVolatility/"


#source(paste0(pre, "Analyses_SamplingDensity_VolCalcs.R"))

# Moving Pictures 
mp_vol_long <- mp_vol_all %>% 
  select(taxID, subjID, 
         AvgTaxAbund = AvgTaxAbund_day1, 
         PropQualChange_day1, PropQualChange_day3, 
         PropQualChange_day7, PropQualChange_day28) %>% 
  pivot_longer(PropQualChange_day1:PropQualChange_day28, 
               names_to = "GapSize", 
               names_prefix = "PropQualChange_", 
               values_to = "PropQualChange")

# Student Microbiome Project
smp_vol_long <- smp_vol_all %>% 
  select(taxID, subjID, 
         AvgTaxAbund = AvgTaxAbund_day7, 
         PropQualChange_day7, PropQualChange_day28) %>% 
  pivot_longer(PropQualChange_day7:PropQualChange_day28, 
               names_to = "GapSize", 
               names_prefix = "PropQualChange_", 
               values_to = "PropQualChange")


# Gajer 
gaj_vol_long <- gaj_vol_all %>% 
  select(taxID, subjID, 
         AvgTaxAbund = AvgTaxAbund_day3, 
         PropQualChange_day3, PropQualChange_day7, PropQualChange_day28) %>% 
  pivot_longer(PropQualChange_day3:PropQualChange_day28, 
               names_to = "GapSize", 
               names_prefix = "PropQualChange_", 
               values_to = "PropQualChange")


# Ravel 
rav_vol_long <- rav_vol_all %>% 
  select(taxID, subjID, 
         AvgTaxAbund = AvgTaxAbund_day1, 
         PropQualChange_day1, PropQualChange_day3, PropQualChange_day7) %>% 
  pivot_longer(PropQualChange_day1:PropQualChange_day7, 
               names_to = "GapSize", 
               names_prefix = "PropQualChange_", 
               values_to = "PropQualChange")


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
  filter(!is.na(PropQualChange))


## Moving pictures, SMP, and Gajer: Day 28 vs. Day 7
## Moving Pictures, Gajer, and Ravel: Day 7 vs. Day 3

vol_wide <- vol_all %>% 
  pivot_wider(id_cols=Study:AvgTaxAbund, names_from=TimeLag, values_from = PropQualChange) 

# Average log fold changes: 7 days vs. 3 days 
vol_wide %>% 
  filter(!Study %in% c("Student Microbiome Project (Gut)", "Ravel (Vaginal)")) %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", "Gajer (Vaginal)"))) %>% 
  pivot_longer(`7 days`:`28 days`, names_to = "SecondTimeLag", values_to = "SecondPropQual") %>% 
  ggplot(aes(x=`3 days`, y = SecondPropQual)) + 
  geom_point(size=0.1, alpha=0.25) + 
  geom_abline(slope=1, intercept=0, lty=2, col="blue") + 
  geom_smooth(aes(col="OLS"), method="lm", se=F) + 
  facet_grid(cols=vars(Study), rows=vars(SecondTimeLag), scales="free") + 
  xlab("3 Day Lag") + ylab("Longer Lag") + 
  ggtitle("Proportion of Time Points with Presence/Absence Change")  + 
  theme_bw() + 
  theme(legend.position="none", text=element_text(size=14)) 


# Average log fold changes: all lags vs. taxon abundance 
vol_all %>% 
  mutate(Study = factor(Study, 
                        levels = c("Moving Pictures (Gut)", "Student Microbiome Project (Gut)", 
                                   "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  ggplot() + 
  theme_bw() + 
  xlab("Log Average Relative Abundance") + 
  ylab("Proportion Qualitative Changes") + 
  geom_point(aes(x=log(AvgTaxAbund), y=PropQualChange,
                 group=TimeLag, color=TimeLag), size=0.1, alpha=0.15) + 
  geom_smooth(aes(x=log(AvgTaxAbund), y=PropQualChange,
                  group=TimeLag, color=TimeLag),  
              method = "loess", se=T) + 
  geom_abline(slope=0, intercept=0, color="black", lty=2) + 
  facet_wrap(vars(Study), nrow=1) + 
  theme(text=element_text(size=14), 
        legend.position = "bottom") + 
  ggtitle("Proportion of Time Points with Presence/Absence Change")


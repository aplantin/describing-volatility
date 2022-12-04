
### Checking whether volatility measures are associated with original read  ###
### count even after rarefying to the minimum read count in that study      ### 

library(tidyverse)
library(vegan)
library(GUniFrac)
library(MBVolDescrip)

source("Analyses_Rarefactions_AssocVolWithReadCt_DataPrep.R")

## Graph showing additive changes by read count and taxon abundance 
all_change_readct %>%
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  mutate(AbsAdditive = abs(AdditiveChg)) %>% 
  ggplot() + 
  xlab("Average Read Count for Time Point Pair") + 
  ylab("Log Absolute Additive Change") + 
  geom_point(aes(x=avg_readct, y=log(AbsAdditive),
                 group=TaxAbundCat, color=TaxAbundCat), size=0.1, alpha=0.15) + 
  geom_smooth(aes(x=avg_readct, y=log(AbsAdditive),
                  group=TaxAbundCat, color=TaxAbundCat),  
              method = "loess", se=T) + 
  theme_bw()  + 
  theme(text=element_text(size=14), 
        legend.position = "bottom") + 
  ggtitle("Log Absolute Additive Changes by Read Count and Taxon Abundance") + 
  facet_wrap(vars(Study), nrow=1, scales="free_x")


## Table of SD of log FC 
logfc_readct_table <- all_change_readct %>%
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  group_by(Study) %>% 
  mutate(ReadCtCat = cut(avg_readct, breaks=4, include.lowest=T, labels=FALSE)) %>% 
  mutate(`Read Count Quartile` = case_when(ReadCtCat == 1 ~ "1_Lowest_Quartile", 
                                           ReadCtCat == 2 ~ "2_Second_Quartile", 
                                           ReadCtCat == 3 ~ "3_Third_Quartile", 
                                           ReadCtCat == 4 ~ "4_Highest_Quartile")) %>% 
  group_by(Study, `Read Count Quartile`, TaxAbundCat) %>% 
  mutate(TaxAbundCat = factor(TaxAbundCat, levels = c("[0, 0.5]", "(0.5, 0.75]", "(0.75, 1]"))) %>% 
  filter(!is.na(NZLogFC)) %>% 
  summarize(sdLogFC = sd(NZLogFC), 
            NinGrp = n()) %>% 
  filter(Study != "Ravel (Vaginal)")
View(logfc_readct_table)


## bar plot showing proportion of qualitative changes by read count quartile and taxon abundance quantile 
all_props_readct = all_change_readct %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", "Ravel (Vaginal)"))) %>% 
  group_by(Study) %>% 
  mutate(ReadCtCat = cut(avg_readct, breaks=4, include.lowest=T, labels=FALSE)) %>% 
  mutate(`Read Count Quartile` = case_when(ReadCtCat == 1 ~ "1_Lowest", 
                               ReadCtCat == 2 ~ "2", 
                               ReadCtCat == 3 ~ "3", 
                               ReadCtCat == 4 ~ "4_Highest")) %>% 
  group_by(Study, `Read Count Quartile`, TaxAbundCat) %>% 
  mutate(TaxAbundCat = factor(TaxAbundCat, levels=c("[0, 0.5]", "(0.5, 0.75]", "(0.75, 1]"))) %>% 
  summarize(PropQualChg = mean(QualChg != 0), 
            NinGrp = n())
all_props_readct %>% 
  ggplot() + 
  geom_bar(aes(x=TaxAbundCat, y=PropQualChg, fill=`Read Count Quartile`), 
           stat="identity", position="dodge") + 
  facet_wrap(vars(Study), nrow=2) + 
  theme_bw() + 
  scale_fill_grey() + 
  xlab("Taxon Abundance Category") + 
  ylab("Proportion of Time Pairs with Qualitative Change") + 
  theme(text=element_text(size=14), legend.position="bottom") + 
  ggtitle("Rate of Qualitative Change by Taxon Abundance and Read Count")



  
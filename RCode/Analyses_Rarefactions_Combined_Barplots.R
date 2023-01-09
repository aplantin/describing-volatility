
library(tidyverse)
library(gridExtra)
library(dplyr)
pre <- "~/AnnaDocuments/Research/Methods_indiv/2022_DescribingVolatility/"

# Long-format changes
longchg_all <- readRDS(paste0(pre, "VolSumms/vol_long_day7_allstudy_allrarefy.rds")) %>% 
  select(Study, Rarefy, taxID, changeID, subjID, AdditiveChg, NZLogFC, CLRChg, AvgTaxAbund) %>% 
  mutate(TaxAbundCat = cut(AvgTaxAbund, c(0, 1e-5, 1e-4, 1e-3, 1), include.lowest=T)) %>% 
  group_by(Study, Rarefy, TaxAbundCat) %>% 
  summarize(sdAddChg = sd(AdditiveChg), 
            sdLogFC = sd(NZLogFC, na.rm=T), 
            sdCLR = sd(CLRChg)) %>% 
  pivot_longer(sdAddChg:sdCLR, names_to="ChangeMeasure", values_to="ChangeSD") %>% 
  mutate(MeasureOfChange = case_when(ChangeMeasure == "sdAddChg" ~ "Additive", 
                                     ChangeMeasure == "sdLogFC" ~ "Log Fold Change", 
                                     ChangeMeasure == "sdCLR" ~ "CLR-Based"))
  
  
p <- longchg_all %>% 
    mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                            "Student Microbiome Project (Gut)", 
                                            "Gajer (Vaginal)", 
                                            "Ravel (Vaginal)")), 
           MeasureOfChange = factor(MeasureOfChange, levels = c("Additive", 
                                                              "Log Fold Change", 
                                                              "CLR-Based")), 
           Rarefy = factor(Rarefy, levels = c("None", "P100", "P80", "P60"))) %>% 
  ggplot(aes(x=TaxAbundCat, y=ChangeSD, fill=Rarefy)) + 
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




png(filename = paste0(pre, "Figures/combined_barplot_rarefy.png"), width = 1600, height = 1200)
p 
dev.off() 


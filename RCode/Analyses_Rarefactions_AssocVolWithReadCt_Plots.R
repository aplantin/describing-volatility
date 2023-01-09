### Checking whether volatility measures are associated with original read  ###
### count even after rarefying to the minimum read count in that study      ### 
library(tidyverse)
library(vegan)
library(GUniFrac)
library(MBVolDescrip)
pre <- "~/AnnaDocuments/Research/Methods_indiv/2022_DescribingVolatility/"



## Get read counts for each sample 
data(mp_otu)
mp_readct <- data.frame(Study = "Moving Pictures (Gut)", 
                        sampID = colnames(mp_otu), 
                        readct = apply(mp_otu, 2, sum)) 


data(smp_otu)
smp_readct <- data.frame(Study = "Student Microbiome Project (Gut)", 
                         sampID = colnames(smp_otu), 
                         readct = apply(smp_otu, 2, sum))

data(gaj_otu)
gaj_readct <- data.frame(Study = "Gajer (Vaginal)", 
                         sampID = rownames(gaj_otu), 
                         readct = apply(gaj_otu, 1, sum)) 

data(rav_otu)
rav_readct <- data.frame(Study = "Ravel (Vaginal)", 
                         sampID = rownames(rav_otu), 
                         readct = apply(rav_otu, 1, sum)) 

all_readct <- do.call(rbind, 
                      list(mp_readct, smp_readct, gaj_readct, rav_readct))



## Add to volatility data (long format) 
longchange_day7 <- readRDS(paste0(pre, "VolSumms/vol_long_day7_allstudy_allrarefy.rds")) %>% 
  merge(., all_readct, by.x = c("Study", "sampID1"), by.y = c("Study", "sampID")) %>% 
  rename(samp1_readct = readct) %>% 
  merge(., all_readct, by.x = c("Study", "sampID2"), by.y = c("Study", "sampID")) %>% 
  rename(samp2_readct = readct)  %>% 
  mutate(avg_readct = (samp1_readct + samp2_readct) / 2) %>% 
  group_by(Study) %>% 
  mutate(avg_readct_cat = cut(avg_readct, 4, label=F)) %>% 
  mutate(`Read Count Quartile` = case_when(avg_readct_cat == 1 ~ "1 (Lowest)", 
                                           avg_readct_cat == 2 ~ "2", 
                                           avg_readct_cat == 3 ~ "3", 
                                           avg_readct_cat == 4 ~ "4 (Highest)"))
  


## SD bar chart (by taxon abundance and original read count quartile)
longchange_day7 %>% 
  select(Study, taxID, changeID, subjID, AdditiveChg, NZLogFC, CLRChg, 
         AvgTaxAbund, avg_readct, `Read Count Quartile`) %>% 
  mutate(TaxAbundCat = cut(AvgTaxAbund, c(0, 1e-5, 1e-4, 1e-3, 1), 
                           include.lowest=T)) %>% 
  group_by(Study, `Read Count Quartile`, TaxAbundCat) %>% 
  summarize(sdAddChg = sd(AdditiveChg), 
            sdLogFC = sd(NZLogFC, na.rm=T), 
            sdCLR = sd(CLRChg)) %>% 
  pivot_longer(sdAddChg:sdCLR, names_to="ChangeMeasure", values_to="ChangeSD") %>% 
  mutate(MeasureOfChange = case_when(ChangeMeasure == "sdAddChg" ~ "Additive", 
                                     ChangeMeasure == "sdLogFC" ~ "Log Fold Change", 
                                     ChangeMeasure == "sdCLR" ~ "CLR-Based")) %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", 
                                          "Ravel (Vaginal)"))) %>% 
  mutate(MeasureOfChange = factor(MeasureOfChange, 
                                  levels = c("Additive", "Log Fold Change", "CLR-Based"))) %>% 
  ggplot(aes(x=TaxAbundCat, y=ChangeSD, fill=`Read Count Quartile`)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_grid(cols=vars(Study), rows=vars(MeasureOfChange), scales="free_y") + 
  scale_x_discrete(guide = guide_axis(n.dodge=2)) + 
  scale_fill_grey() + 
  theme_bw() + 
  theme(text=element_text(size=12)) + 
  xlab("Taxon Abundance Category") + 
  ylab("Standard Deviation of Change Measure")


## line graph showing proportion of qualitative changes by read count quartile and taxon abundance quantile 
longchange_day7 %>% 
  mutate(Study = factor(Study, levels = c("Moving Pictures (Gut)", 
                                          "Student Microbiome Project (Gut)", 
                                          "Gajer (Vaginal)", 
                                          "Ravel (Vaginal)"))) %>%
  group_by(Study, `Read Count Quartile`, taxID, AvgTaxAbund) %>% 
  summarize(PropQualChg = mean(QualChg != 0), NinGrp = n())  %>% 
  ggplot() + 
  geom_point(aes(x=log(AvgTaxAbund), y=PropQualChg, group=`Read Count Quartile`, 
                 col=`Read Count Quartile`), size=0.1, alpha=0.15) + 
  geom_smooth(aes(x=log(AvgTaxAbund), y=PropQualChg, group=`Read Count Quartile`, 
                  col=`Read Count Quartile`),  
              method = "loess", se=T) + 
  geom_abline(slope=0, intercept=0, color="black", lty=2) + 
  facet_wrap(vars(Study), nrow=1, scales = "free_x") + 
  theme_bw() + 
  xlab("Taxon Abundance") + 
  ylab("Proportion Qualitative Changes") + 
  ggtitle("Proportion of Time Points with Presence/Absence Change") + 
  ylim(0,1) + 
  theme(text=element_text(size=14), legend.position="bottom") 



  
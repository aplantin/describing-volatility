
### Checking whether volatility measures are associated with original read  ###
### count even after rarefying to the minimum read count in that study      ### 

library(tidyverse)
library(vegan)
library(GUniFrac)
library(MBVolDescrip)


## Moving pictures 
data(mp_otu); data(mp_meta)
mp_meta$readct <- apply(mp_otu, 2, sum)
mp_meta$subjID = as.character(mp_meta$subjID)

# Rarefy to constant sampling depth 
set.seed(1) 
mp_otu_rarefy <- rrarefy(t(mp_otu), sample=min(mp_meta$readct))

# Temporal change metadata 
mp_changemeta <- temporalSubsampleMeta(mp_otu_rarefy, mp_meta, desired_spacing = 7, 
                        window_width = 1, taxaAreRows = FALSE) 

# Calculate microbiome changes at different time gaps 
mp_changemats <- calcMicrobiomeChanges(mp_otu_rarefy, mp_changemeta, Ds = NULL, taxaAreRows = FALSE)

# Long format - all changes  
mp_longchange <- longMicrobiomeChanges(mbchanges = mp_changemats)
mp_change_readct <- merge(mp_longchange, mp_meta[, c("subjID", "sampID", "readct")], 
                          by.x=c("subjID", "sampID1"), by.y=c("subjID", "sampID")) %>% 
  rename(samp1_readct = readct) %>% 
  merge(., mp_meta[, c("subjID", "sampID", "readct")], 
        by.x=c("subjID", "sampID2"), by.y=c("subjID", "sampID")) %>% 
  rename(samp2_readct = readct) %>% 
  mutate(avg_readct = (samp1_readct + samp2_readct)/2)

# Calculate volatility measures 
bks <- c(0, quantile(mp_change_readct$AvgTaxAbund, c(0.5, 0.75, 1)))
mp_change_readct <- mp_change_readct %>% 
  mutate(TaxAbundCuts = cut(AvgTaxAbund, breaks=bks)) %>% 
  mutate(TaxAbundCat = case_when(TaxAbundCuts == "(0,5.3e-06]" ~ "[0, 0.5]", 
                                 TaxAbundCuts == "(5.3e-06,3.41e-05]" ~ "(0.5, 0.75]", 
                                 TaxAbundCuts == "(3.41e-05,0.16]" ~ "(0.75, 1]"))
head(mp_change_readct)



## Student Microbiome Project
data(smp_otu); data(smp_meta)
smp_meta$readct <- apply(smp_otu, 2, sum)
smp_meta$subjID = as.character(smp_meta$subjID)

# Rarefy to constant sampling depth 
set.seed(1) 
smp_otu_rarefy <- rrarefy(t(smp_otu), sample=min(smp_meta$readct))

# Temporal change metadata 
smp_changemeta <- temporalSubsampleMeta(smp_otu_rarefy, smp_meta, desired_spacing = 7, 
                                       window_width = 1, taxaAreRows = FALSE) 

# Calculate microbiome changes at different time gaps 
smp_changemats <- calcMicrobiomeChanges(smp_otu_rarefy, smp_changemeta, 
                                        Ds = NULL, taxaAreRows = FALSE)

# Long format - all changes  
smp_longchange <- longMicrobiomeChanges(mbchanges = smp_changemats)
smp_change_readct <- merge(smp_longchange, smp_meta[, c("subjID", "sampID", "readct")], 
                          by.x=c("subjID", "sampID1"), by.y=c("subjID", "sampID")) %>% 
  rename(samp1_readct = readct) %>% 
  merge(., smp_meta[, c("subjID", "sampID", "readct")], 
        by.x=c("subjID", "sampID2"), by.y=c("subjID", "sampID")) %>% 
  rename(samp2_readct = readct) %>% 
  mutate(avg_readct = (samp1_readct + samp2_readct)/2)

# Calculate volatility measures 
bks <- c(0, quantile(smp_change_readct$AvgTaxAbund, c(0.5, 0.75, 1)))
smp_change_readct <- smp_change_readct %>% 
  mutate(TaxAbundCuts = cut(AvgTaxAbund, breaks=bks)) %>% 
  mutate(TaxAbundCat = case_when(TaxAbundCuts == "(0,1.86e-06]" ~ "[0, 0.5]", 
                                 TaxAbundCuts == "(1.86e-06,5.03e-05]" ~ "(0.5, 0.75]", 
                                 TaxAbundCuts == "(5.03e-05,0.24]" ~ "(0.75, 1]"))
head(smp_change_readct)




## Gajer
data(gaj_otu); data(gaj_meta)
gaj_meta$readct <- apply(gaj_otu, 1, sum)
gaj_meta$subjID = as.character(gaj_meta$subjID)

# Rarefy to constant sampling depth 
set.seed(1) 
gaj_otu_rarefy <- rrarefy(gaj_otu, sample=min(gaj_meta$readct))

# Temporal change metadata 
gaj_changemeta <- temporalSubsampleMeta(gaj_otu_rarefy, gaj_meta, desired_spacing = 7, 
                                        window_width = 1, taxaAreRows = FALSE) 

# Calculate microbiome changes at different time gaps 
gaj_changemats <- calcMicrobiomeChanges(gaj_otu_rarefy, gaj_changemeta, 
                                        Ds = NULL, taxaAreRows = FALSE)

# Long format - all changes  
gaj_longchange <- longMicrobiomeChanges(mbchanges = gaj_changemats)
gaj_change_readct <- merge(gaj_longchange, gaj_meta[, c("subjID", "sampID", "readct")], 
                           by.x=c("subjID", "sampID1"), by.y=c("subjID", "sampID")) %>% 
  rename(samp1_readct = readct) %>% 
  merge(., gaj_meta[, c("subjID", "sampID", "readct")], 
        by.x=c("subjID", "sampID2"), by.y=c("subjID", "sampID")) %>% 
  rename(samp2_readct = readct) %>% 
  mutate(avg_readct = (samp1_readct + samp2_readct)/2)

# Calculate volatility measures 
bks <- c(0, quantile(gaj_change_readct$AvgTaxAbund, c(0.5, 0.75, 1)))
gaj_change_readct <- gaj_change_readct %>% 
  mutate(TaxAbundCuts = cut(AvgTaxAbund, breaks=bks)) %>% 
  mutate(TaxAbundCat = case_when(TaxAbundCuts == "(0,3.43e-05]" ~ "[0, 0.5]", 
                                 TaxAbundCuts == "(3.43e-05,0.000329]" ~ "(0.5, 0.75]", 
                                 TaxAbundCuts == "(0.000329,0.357]" ~ "(0.75, 1]"))
head(gaj_change_readct)



## Ravel
data(rav_otu); data(rav_meta)
rav_meta$readct <- apply(rav_otu, 1, sum)
rav_meta$subjID = as.character(rav_meta$subjID)

# Rarefy to constant sampling depth 
set.seed(1) 
rav_otu_rarefy <- rrarefy(rav_otu, sample=min(rav_meta$readct))

# Temporal change metadata 
rav_changemeta <- temporalSubsampleMeta(rav_otu_rarefy, rav_meta, desired_spacing = 7, 
                                        window_width = 1, taxaAreRows = FALSE) 

# Calculate microbiome changes at different time gaps 
rav_changemats <- calcMicrobiomeChanges(rav_otu_rarefy, rav_changemeta, 
                                        Ds = NULL, taxaAreRows = FALSE)

# Long format - all changes  
rav_longchange <- longMicrobiomeChanges(mbchanges = rav_changemats)
rav_change_readct <- merge(rav_longchange, rav_meta[, c("subjID", "sampID", "readct")], 
                           by.x=c("subjID", "sampID1"), by.y=c("subjID", "sampID")) %>% 
  rename(samp1_readct = readct) %>% 
  merge(., rav_meta[, c("subjID", "sampID", "readct")], 
        by.x=c("subjID", "sampID2"), by.y=c("subjID", "sampID")) %>% 
  rename(samp2_readct = readct) %>% 
  mutate(avg_readct = (samp1_readct + samp2_readct)/2)

# Calculate volatility measures 
bks <- c(0, quantile(rav_change_readct$AvgTaxAbund, c(0.5, 0.75, 1)))
rav_change_readct <- rav_change_readct %>% 
  mutate(TaxAbundCuts = cut(AvgTaxAbund, breaks=bks))  %>% 
  mutate(TaxAbundCat = case_when(TaxAbundCuts == "(0,0.00108]" ~ "[0, 0.5]", 
                                 TaxAbundCuts == "(0.00108,0.00431]" ~ "(0.5, 0.75]", 
                                 TaxAbundCuts == "(0.00431,0.427]" ~ "(0.75, 1]"))


all_change_readct <- rbind(cbind(Study = "Moving Pictures (Gut)", mp_change_readct), 
                           cbind(Study = "Student Microbiome Project (Gut)", smp_change_readct), 
                           cbind(Study = "Gajer (Vaginal)", gaj_change_readct), 
                           cbind(Study = "Ravel (Vaginal)", rav_change_readct)) 

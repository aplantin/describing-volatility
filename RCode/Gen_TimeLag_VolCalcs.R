#############
### Setup ###
#############

library(tidyverse)
library(vegan)
library(GUniFrac)
#devtools::install_github("aplantin/describing-volatility", subdir="MBVolDescrip")
library(MBVolDescrip)
pre <- "~/AnnaDocuments/Research/Methods_indiv/2022_DescribingVolatility/"



############################
### Moving pictures data ###
############################

# Note: in each case data have been pre-processed to have convenient column names 
# sampID, subjID, and time; ordered by subjID and time; and order of samples 
# has been confirmed to match between metadata and OTU table 

## Moving pictures 
## Original sampling frequency: daily 
# taxa are rows, samples are columns
data(mp_otu); data(mp_meta); data(mp_tree)
mp_meta$readct <- apply(mp_otu, 2, sum)
mp_meta$subjID = as.character(mp_meta$subjID)

# Rarefy to constant sampling depth 
set.seed(1) 
mp_otu_rarefy <- rrarefy(t(mp_otu), sample=min(mp_meta$readct))

# Temporal change metadata 
mp_changemeta <- lapply(c(1,3,7,28), FUN = function(gap) {
  temporalSubsampleMeta(mp_otu_rarefy, mp_meta, desired_spacing = gap, 
                        window_width = floor(gap/7), taxaAreRows = FALSE) 
})

# Calculate distance matrices 
mp_unifracs <- GUniFrac(mp_otu_rarefy, mp_tree, alpha=c(0, 0.5, 1))$unifracs
mp_Ds <- list(unweightedUF = mp_unifracs[,,"d_UW"], 
              genUF0.5 = mp_unifracs[,,"d_0.5"], 
              weightedUF = mp_unifracs[,,"d_1"], 
              BrayCurtis = as.matrix(vegdist(mp_otu_rarefy, method="bray"))) 

# Calculate microbiome changes at different time gaps 
mp_changemats <- lapply(mp_changemeta, FUN = function(cmet) {
  calcMicrobiomeChanges(mp_otu_rarefy, cmet, Ds = mp_Ds, taxaAreRows = FALSE)
})

# Long format - all changes  
mp_longchange <- lapply(mp_changemats, FUN = function(cmat) {
  longMicrobiomeChanges(mbchanges = cmat)
})

# Calculate volatility summaries 
mp_vol <- lapply(mp_changemats, FUN = function(cmat) {
  summMicrobiomeVolatility(mbchanges = cmat)
})

# Merge different time scales (for summaries by taxon & subject)
mp_vol_1_3 <- merge(mp_vol[[1]], mp_vol[[2]], by=c("subjID", "taxID"), 
                    suffixes = c("_day1", "_day3", all=T))
mp_vol_7_28 <- merge(mp_vol[[3]], mp_vol[[4]], by=c("subjID", "taxID"), 
                     suffixes = c("_day7", "_day28"), all=T)
mp_vol_all <- merge(mp_vol_1_3, mp_vol_7_28, by=c("subjID", "taxID"), all=T)

# Merge different time scales (for summaries by subject only -- distance based)
mp_distvol <- do.call(rbind, lapply(mp_changemats, FUN = function(l) l$ChangeMeta)) %>% 
  mutate(GapSize = case_when(abs(time2 - time1) == 1 ~ 1, 
                             abs(time2 - time1) == 3 ~ 3, 
                             abs(time2 - time1) %in% 6:8 ~ 7, 
                             abs(time2 - time1) %in% 24:32 ~ 28)) %>% 
  select(subjID, changeID, GapSize, unweightedUF:BrayCurtis)




#######################################
### Student Microbiome Project data ###
#######################################

## Student Microbiome Project
## Original sampling frequency: weekly 
# taxa are rows, samples are columns 
data(smp_otu); data(smp_meta)
smp_meta$readct <- apply(smp_otu, 2, sum)
smp_meta$subjID = as.character(smp_meta$subjID)

# Rarefy to constant sampling depth 
set.seed(1) 
smp_otu_rarefy <- rrarefy(t(smp_otu), sample=min(smp_meta$readct))

# Temporal change metadata 
smp_changemeta <- lapply(c(7,28), FUN = function(gap) {
  temporalSubsampleMeta(smp_otu_rarefy, smp_meta, desired_spacing = gap, 
                        window_width = floor(gap/7), taxaAreRows = FALSE) 
})

## Calculate distance matrices 
smp_D_BC <- as.matrix(vegdist(smp_otu_rarefy, method="bray"))

# Calculate microbiome changes at different time gaps 
smp_changemats <- lapply(smp_changemeta, FUN = function(cmet) {
  calcMicrobiomeChanges(smp_otu_rarefy, cmet, Ds = smp_D_BC, taxaAreRows = FALSE)
})

# Long format - all changes  
smp_longchange <- lapply(smp_changemats, FUN = function(cmat) {
  longMicrobiomeChanges(mbchanges = cmat)
})

# Calculate volatility summaries 
smp_vol <- lapply(smp_changemats, FUN = function(cmat) {
  summMicrobiomeVolatility(mbchanges = cmat)
})

# Merge different time scales (for summaries by taxon & subject))
smp_vol_all <- merge(smp_vol[[1]], smp_vol[[2]], by=c("subjID", "taxID"), 
                     suffixes = c("_day7", "_day28"), all=T)

# Merge different time scales (for summaries by subject only -- distance based)
smp_distvol <- do.call(rbind, lapply(smp_changemats, FUN = function(l) l$ChangeMeta)) %>% 
  mutate(GapSize = case_when(abs(time2 - time1) == 1 ~ 1, 
                             abs(time2 - time1) == 3 ~ 3, 
                             abs(time2 - time1) %in% 6:8 ~ 7, 
                             abs(time2 - time1) %in% 24:32 ~ 28)) %>% 
  rename(BrayCurtis=DistanceMetric) %>% 
  select(subjID, changeID, GapSize, BrayCurtis)




#####################################
### Gajer Vaginal Microbiome data ###
#####################################

## Gajer Vaginal Microbiome 
## Original sampling frequency: twice-weekly for up to 16 weeks
# samples are rows, taxa are columns  
data(gaj_otu); data(gaj_meta)
gaj_meta$readct <- apply(gaj_otu, 1, sum)
gaj_meta$subjID = as.character(gaj_meta$subjID)

# Rarefy to constant sampling depth 
set.seed(1) 
gaj_otu_rarefy <- rrarefy(gaj_otu, sample=min(gaj_meta$readct))

# Temporal change metadata 
gaj_changemeta <- lapply(c(3,7,28), FUN = function(gap) {
  temporalSubsampleMeta(gaj_otu_rarefy, gaj_meta, desired_spacing = gap, 
                        window_width = floor(gap/7), taxaAreRows = FALSE) 
})

## Calculate distance matrices 
gaj_D_BC <- as.matrix(vegdist(gaj_otu_rarefy, method="bray"))

# Calculate microbiome changes at different time gaps 
gaj_changemats <- lapply(gaj_changemeta, FUN = function(cmet) {
  calcMicrobiomeChanges(gaj_otu_rarefy, cmet, Ds = gaj_D_BC, taxaAreRows = FALSE)
})

# Long format - all changes  
gaj_longchange <- lapply(gaj_changemats, FUN = function(cmat) {
  longMicrobiomeChanges(mbchanges = cmat)
})

# Calculate volatility summaries 
gaj_vol <- lapply(gaj_changemats, FUN = function(cmat) {
  summMicrobiomeVolatility(mbchanges = cmat)
})

# Merge different time scales (for summaries by taxon & subject)
gaj_vol_3_7 <- merge(gaj_vol[[1]], gaj_vol[[2]], by=c("subjID", "taxID"), 
                     suffixes = c("_day3", "_day7"), all=T)
gaj_vol_7_28 <- merge(gaj_vol[[2]], gaj_vol[[3]], by=c("subjID", "taxID"), 
                      suffixes = c("_day7", "_day28"), all=T)
gaj_vol_all <- merge(gaj_vol_3_7, gaj_vol_7_28, all=T)

# Merge different time scales (for summaries by subject only -- distance based)
gaj_distvol <- do.call(rbind, lapply(gaj_changemats, FUN = function(l) l$ChangeMeta)) %>% 
  mutate(GapSize = case_when(abs(time2 - time1) == 1 ~ 1, 
                             abs(time2 - time1) == 3 ~ 3, 
                             abs(time2 - time1) %in% 6:8 ~ 7, 
                             abs(time2 - time1) %in% 24:32 ~ 28)) %>% 
  rename(BrayCurtis=DistanceMetric) %>% 
  select(subjID, changeID, GapSize, BrayCurtis)




#####################################
### Ravel Vaginal Microbiome data ###
#####################################

## Ravel Vaginal Microbiome 
## Original sampling frequency: daily for up to 1.5 months 
# samples are rows, taxa are columns  
data(rav_otu); data(rav_meta)
rav_meta$readct <- apply(rav_otu, 1, sum)
rav_meta$subjID = as.character(rav_meta$subjID)

# Rarefy to constant sampling depth 
set.seed(1) 
rav_otu_rarefy <- rrarefy(rav_otu, sample=500) ## sample=min(rav_meta$readct))

# Temporal change metadata 
rav_changemeta <- lapply(c(1,3,7), FUN = function(gap) {
  temporalSubsampleMeta(rav_otu_rarefy, rav_meta, desired_spacing = gap, 
                        window_width = floor(gap/7), taxaAreRows = FALSE) 
})

## Calculate distance matrices 
rav_D_BC <- as.matrix(vegdist(rav_otu_rarefy, method="bray"))

# Calculate microbiome changes at different time gaps 
rav_changemats <- lapply(rav_changemeta, FUN = function(cmet) {
  calcMicrobiomeChanges(rav_otu_rarefy, cmet, Ds = rav_D_BC, taxaAreRows = FALSE)
})

# Long format - all changes  
rav_longchange <- lapply(rav_changemats, FUN = function(cmat) {
  longMicrobiomeChanges(mbchanges = cmat)
})

# Calculate volatility summaries 
rav_vol <- lapply(rav_changemats, FUN = function(cmat) {
  summMicrobiomeVolatility(mbchanges = cmat)
})

# Merge different time scales (for summaries by taxon & subject)
rav_vol_1_3 <- merge(rav_vol[[1]], rav_vol[[2]], by=c("subjID", "taxID"), 
                     suffixes = c("_day1", "_day3"), all=T)
rav_vol_3_7 <- merge(rav_vol[[2]], rav_vol[[3]], by=c("subjID", "taxID"), 
                     suffixes = c("_day3", "_day7"), all=T)
rav_vol_all <- merge(rav_vol_1_3, rav_vol_3_7, all=T)

# Merge different time scales (for summaries by subject only -- distance based)
rav_distvol <- do.call(rbind, lapply(rav_changemats, FUN = function(l) l$ChangeMeta)) %>% 
  mutate(GapSize = case_when(abs(time2 - time1) == 1 ~ 1, 
                             abs(time2 - time1) == 3 ~ 3, 
                             abs(time2 - time1) %in% 6:8 ~ 7, 
                             abs(time2 - time1) %in% 24:32 ~ 28)) %>% 
  rename(BrayCurtis=DistanceMetric) %>% 
  select(subjID, changeID, GapSize, BrayCurtis)




#### 
saveRDS(mp_vol_all, "VolSumms/mp_vol_all_rarefyP100.rds")
saveRDS(mp_longchange, "VolSumms/mp_longchg_rarefyP100.rds")
saveRDS(mp_distvol, "VolSumms/mp_distvol_rarefyP100.rds")

saveRDS(smp_vol_all, "VolSumms/smp_vol_all_rarefyP100.rds")
saveRDS(smp_longchange, "VolSumms/smp_longchg_rarefyP100.rds")
saveRDS(smp_distvol, "VolSumms/smp_distvol_rarefyP100.rds")

saveRDS(gaj_vol_all, "VolSumms/gaj_vol_all_rarefyP100.rds")
saveRDS(gaj_longchange, "VolSumms/gaj_longchg_rarefyP100.rds")
saveRDS(gaj_distvol, "VolSumms/gaj_distvol_rarefyP100.rds")

saveRDS(rav_vol_all, "VolSumms/rav_vol_all_rarefyP100.rds")
saveRDS(rav_longchange, "VolSumms/rav_longchg_rarefyP100.rds")
saveRDS(rav_distvol, "VolSumms/rav_distvol_rarefyP100.rds")



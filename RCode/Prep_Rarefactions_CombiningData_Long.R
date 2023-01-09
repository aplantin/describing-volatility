## Moving Pictures 
# get set up with no rarefaction dataset 
mp_vol_norare <- readRDS("VolSumms/mp_longchg_rarefyNone.rds")
mp_whichday7 <- which(unlist(lapply(mp_vol_norare, 
                                    FUN = function(x) all(x$NominalLag == 7))))
mp_vol_norare <- mp_vol_norare[[mp_whichday7]] %>% 
  mutate(Rarefy = "None")

# add other three rarefactions 
mp_vol_rp100 <- readRDS("VolSumms/mp_longchg_rarefyP100.rds")[[mp_whichday7]] %>% 
  mutate(Rarefy = "P100")
mp_vol_rp80 <- readRDS("VolSumms/mp_longchg_rarefyP80.rds")[[mp_whichday7]] %>% 
  mutate(Rarefy = "P80")
mp_vol_rp60 <- readRDS("VolSumms/mp_longchg_rarefyP60.rds")[[mp_whichday7]] %>% 
  mutate(Rarefy = "P60")

# append into one long dataset 
mp_longchg_all <- mp_vol_norare %>% 
  rbind(., mp_vol_rp100) %>% 
  rbind(., mp_vol_rp80) %>% 
  rbind(., mp_vol_rp60) %>% 
  mutate(Study = "Moving Pictures (Gut)")



## Student Microbiome Project
# get set up with no rarefaction dataset 
smp_vol_norare <- readRDS("VolSumms/smp_longchg_rarefyNone.rds")
smp_whichday7 <- which(unlist(lapply(smp_vol_norare, 
                                    FUN = function(x) all(x$NominalLag == 7))))
smp_vol_norare <- smp_vol_norare[[smp_whichday7]] %>% 
  mutate(Rarefy = "None")

# add other three rarefactions 
smp_vol_rp100 <- readRDS("VolSumms/smp_longchg_rarefyP100.rds")[[smp_whichday7]] %>% 
  mutate(Rarefy = "P100")
smp_vol_rp80 <- readRDS("VolSumms/smp_longchg_rarefyP80.rds")[[smp_whichday7]] %>% 
  mutate(Rarefy = "P80")
smp_vol_rp60 <- readRDS("VolSumms/smp_longchg_rarefyP60.rds")[[smp_whichday7]] %>% 
  mutate(Rarefy = "P60")

# append into one long dataset 
smp_longchg_all <- smp_vol_norare %>% 
  rbind(., smp_vol_rp100) %>% 
  rbind(., smp_vol_rp80) %>% 
  rbind(., smp_vol_rp60) %>% 
  mutate(Study = "Student Microbiome Project (Gut)")



## Ravel
# get set up with no rarefaction dataset 
rav_vol_norare <- readRDS("VolSumms/rav_longchg_rarefyNone.rds")
rav_whichday7 <- which(unlist(lapply(rav_vol_norare, 
                                    FUN = function(x) all(x$NominalLag == 7))))
rav_vol_norare <- rav_vol_norare[[rav_whichday7]] %>% 
  mutate(Rarefy = "None")

# add other three rarefactions 
rav_vol_rp100 <- readRDS("VolSumms/rav_longchg_rarefyP100.rds")[[rav_whichday7]] %>% 
  mutate(Rarefy = "P100")
rav_vol_rp80 <- readRDS("VolSumms/rav_longchg_rarefyP80.rds")[[rav_whichday7]] %>% 
  mutate(Rarefy = "P80")
rav_vol_rp60 <- readRDS("VolSumms/rav_longchg_rarefyP60.rds")[[rav_whichday7]] %>% 
  mutate(Rarefy = "P60")

# append into one long dataset 
rav_longchg_all <- rav_vol_norare %>% 
  rbind(., rav_vol_rp100) %>% 
  rbind(., rav_vol_rp80) %>% 
  rbind(., rav_vol_rp60) %>% 
  mutate(Study = "Ravel (Vaginal)")




## Gajer
# get set up with no rarefaction dataset 
gaj_vol_norare <- readRDS("VolSumms/gaj_longchg_rarefyNone.rds")
gaj_whichday7 <- which(unlist(lapply(gaj_vol_norare, 
                                     FUN = function(x) all(x$NominalLag == 7))))
gaj_vol_norare <- gaj_vol_norare[[gaj_whichday7]] %>% 
  mutate(Rarefy = "None")

# add other three rarefactions 
gaj_vol_rp100 <- readRDS("VolSumms/gaj_longchg_rarefyP100.rds")[[gaj_whichday7]] %>% 
  mutate(Rarefy = "P100")
gaj_vol_rp80 <- readRDS("VolSumms/gaj_longchg_rarefyP80.rds")[[gaj_whichday7]] %>% 
  mutate(Rarefy = "P80")
gaj_vol_rp60 <- readRDS("VolSumms/gaj_longchg_rarefyP60.rds")[[gaj_whichday7]] %>% 
  mutate(Rarefy = "P60")

# append into one long dataset 
gaj_longchg_all <- gaj_vol_norare %>% 
  rbind(., gaj_vol_rp100) %>% 
  rbind(., gaj_vol_rp80) %>% 
  rbind(., gaj_vol_rp60) %>% 
  mutate(Study = "Gajer (Vaginal)")




## Combine all 
vol_all_rarefy <- rbind(mp_longchg_all, smp_longchg_all, 
                        rav_longchg_all, gaj_longchg_all)

## Save 
saveRDS(vol_all_rarefy, file="VolSumms/vol_long_day7_allstudy_allrarefy.rds")


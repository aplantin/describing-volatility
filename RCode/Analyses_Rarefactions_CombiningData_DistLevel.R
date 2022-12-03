
## No rarefaction 
mp_vol_norare <- readRDS("VolSumms/mp_distvol_rarefyNone.rds") 
smp_vol_norare <- readRDS("VolSumms/smp_distvol_rarefyNone.rds") 
gaj_vol_norare <- readRDS("VolSumms/gaj_distvol_rarefyNone.rds") 
rav_vol_norare <- readRDS("VolSumms/rav_distvol_rarefyNone.rds") 

## 100% rarefaction 
mp_vol_rp100 <- readRDS("VolSumms/mp_distvol_rarefyP100.rds")
smp_vol_rp100 <- readRDS("VolSumms/smp_distvol_rarefyP100.rds") 
gaj_vol_rp100 <- readRDS("VolSumms/gaj_distvol_rarefyP100.rds") 
rav_vol_rp100 <- readRDS("VolSumms/rav_distvol_rarefyP100.rds")

## 100% rarefaction 
mp_vol_rp80 <- readRDS("VolSumms/mp_distvol_rarefyP80.rds") 
smp_vol_rp80 <- readRDS("VolSumms/smp_distvol_rarefyP80.rds") 
gaj_vol_rp80 <- readRDS("VolSumms/gaj_distvol_rarefyP80.rds")
rav_vol_rp80 <- readRDS("VolSumms/rav_distvol_rarefyP80.rds") 

## 100% rarefaction 
mp_vol_rp60 <- readRDS("VolSumms/mp_distvol_rarefyP60.rds") 
smp_vol_rp60 <- readRDS("VolSumms/smp_distvol_rarefyP60.rds")
gaj_vol_rp60 <- readRDS("VolSumms/gaj_distvol_rarefyP60.rds") 
rav_vol_rp60 <- readRDS("VolSumms/rav_distvol_rarefyP60.rds") 


## Combine all MP data (all distances) 
distvol_mp_rarefy <- rbind(cbind(Study = "Moving Pictures (Gut)", Rarefy = "None", mp_vol_norare), 
                           cbind(Study = "Moving Pictures (Gut)", Rarefy = "P100", mp_vol_rp100), 
                           cbind(Study = "Moving Pictures (Gut)", Rarefy = "P80", mp_vol_rp80), 
                           cbind(Study = "Moving Pictures (Gut)", Rarefy = "P60", mp_vol_rp60))

## Combine all studies (BrayCurtis only)
distvol_rarefy <- rbind(distvol_mp_rarefy[, c("Study", "Rarefy", "subjID", "changeID", "GapSize", "BrayCurtis")], 
                        cbind(Study = "Student Microbiome Project (Gut)", Rarefy = "None", smp_vol_norare),
                        cbind(Study = "Student Microbiome Project (Gut)", Rarefy = "P100", smp_vol_rp100), 
                        cbind(Study = "Student Microbiome Project (Gut)", Rarefy = "P80", smp_vol_rp80), 
                        cbind(Study = "Student Microbiome Project (Gut)", Rarefy = "P60", smp_vol_rp60), 
                        cbind(Study = "Gajer (Vaginal)", Rarefy = "None", gaj_vol_norare), 
                        cbind(Study = "Gajer (Vaginal)", Rarefy = "P100", gaj_vol_rp100), 
                        cbind(Study = "Gajer (Vaginal)", Rarefy = "P80", gaj_vol_rp80), 
                        cbind(Study = "Gajer (Vaginal)", Rarefy = "P60", gaj_vol_rp60), 
                        cbind(Study = "Ravel (Vaginal)", Rarefy = "None", rav_vol_norare), 
                        cbind(Study = "Ravel (Vaginal)", Rarefy = "P100", rav_vol_rp100), 
                        cbind(Study = "Ravel (Vaginal)", Rarefy = "P80", rav_vol_rp80), 
                        cbind(Study = "Ravel (Vaginal)", Rarefy = "P60", rav_vol_rp60))  

## Save 
saveRDS(distvol_rarefy, file="VolSumms/distvol_rarefy.rds")
saveRDS(distvol_mp_rarefy, file="VolSumms/distvol_rarefy_mp.rds")



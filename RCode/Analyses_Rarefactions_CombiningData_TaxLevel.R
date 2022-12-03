
## No rarefaction 
mp_vol_norare <- readRDS("VolSumms/mp_vol_all_rarefyNone.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))
smp_vol_norare <- readRDS("VolSumms/smp_vol_all_rarefyNone.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))
gaj_vol_norare <- readRDS("VolSumms/gaj_vol_all_rarefyNone.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))
rav_vol_norare <- readRDS("VolSumms/rav_vol_all_rarefyNone.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))

## 100% rarefaction 
mp_vol_rp100 <- readRDS("VolSumms/mp_vol_all_rarefyP100.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))
smp_vol_rp100 <- readRDS("VolSumms/smp_vol_all_rarefyP100.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))
gaj_vol_rp100 <- readRDS("VolSumms/gaj_vol_all_rarefyP100.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))
rav_vol_rp100 <- readRDS("VolSumms/rav_vol_all_rarefyP100.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))

## 100% rarefaction 
mp_vol_rp80 <- readRDS("VolSumms/mp_vol_all_rarefyP80.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))
smp_vol_rp80 <- readRDS("VolSumms/smp_vol_all_rarefyP80.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))
gaj_vol_rp80 <- readRDS("VolSumms/gaj_vol_all_rarefyP80.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))
rav_vol_rp80 <- readRDS("VolSumms/rav_vol_all_rarefyP80.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))

## 100% rarefaction 
mp_vol_rp60 <- readRDS("VolSumms/mp_vol_all_rarefyP60.rds") %>% 
  select(subjID, taxID, ends_with("_day7")) 
smp_vol_rp60 <- readRDS("VolSumms/smp_vol_all_rarefyP60.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))
gaj_vol_rp60 <- readRDS("VolSumms/gaj_vol_all_rarefyP60.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))
rav_vol_rp60 <- readRDS("VolSumms/rav_vol_all_rarefyP60.rds") %>% 
  select(subjID, taxID, ends_with("_day7"))


## Combine all 
vol_all_rarefy <- rbind(cbind(Study = "Moving Pictures (Gut)", Rarefy = "None", mp_vol_norare), 
                        cbind(Study = "Moving Pictures (Gut)", Rarefy = "P100", mp_vol_rp100), 
                        cbind(Study = "Moving Pictures (Gut)", Rarefy = "P80", mp_vol_rp80), 
                        cbind(Study = "Moving Pictures (Gut)", Rarefy = "P60", mp_vol_rp60), 
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
saveRDS(vol_all_rarefy, file="VolSumms/vol_all_rarefy.rds")




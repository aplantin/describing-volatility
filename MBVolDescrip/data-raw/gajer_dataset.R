## code to process Gajer healthy vaginal microbiome dataset
library(tidyverse) 

filepath.metadata <- system.file("extdata", "gajer_metadata.csv", package="MBVolDescrip")
filepath.otus <- system.file("extdata", "gajer_OTUs.csv", package="MBVolDescrip")

## All data in Table S2 of Gajer article 
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3722878/
## doi: https://doi.org/10.1126/scitranslmed.3003605 

# Metadata 
# Interpretable values for previously-coded variables (e.g., race) 
gv_meta_raw <- read.csv(filepath.metadata) 
gaj_meta <- gv_meta_raw %>% 
  mutate(Race = case_when(Racea == 1 ~ "White", 
                          Racea == 0 ~ "Black", 
                          Racea == 5 ~ "Hispanic", 
                          Racea == 4 ~ "Others"), 
         subjID = paste0("subj_", Subject.ID)) %>% 
  dplyr::rename(sampID = Sample.ID, 
                time = Time.in.study, 
                NugentScore = Nugent.Score, 
                NugentCat = Nugent.Categoryb, 
                CST = Community.State.Typec, 
                ReadCount = Total.Read.Countsd) %>% 
  arrange(subjID, time)


# Load proportions (from article) and check match of sample IDs 
gv_otus_props <- read.csv(filepath.otus, row.names = 1)
gv_otus_props <- gv_otus_props[gaj_meta$sampID, ]
all(rownames(gv_otus_props) == gaj_meta$sampID) 

# Convert to counts -- rounded to whole numbers 
gaj_otu <- gv_otus_props 
for (i in 1:nrow(gaj_otu)) {
  gaj_otu[i, ] <- round(gv_otus_props[i, ] * gaj_meta$ReadCount/100)
}

# Summary stats: number of samples per subject
table(table(gaj_meta$subjID))
median(table(gaj_meta$subjID))

# Summary stats: number of OTUs per sample 
# Taxa are columns 
summary(gaj_meta$ReadCount)
dim(gaj_otu)

# save as package data 
usethis::use_data(gaj_meta, overwrite = TRUE)
usethis::use_data(gaj_otu, overwrite = TRUE)

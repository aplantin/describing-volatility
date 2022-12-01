## code to prepare Moving Pictures dataset 
library(tidyverse) 
library(MicrobeDS)
data("MovingPictures")

# load all components of Moving Pictures dataset 
mp_meta_orig <- data.frame(sample_data(MovingPictures))
mp_otu_orig <- data.frame(otu_table(MovingPictures))
mp_tax <- data.frame(tax_table(MovingPictures))
mp_tree <- phy_tree(MovingPictures)

# filter to gut microbiome time series;  sort by subject ID and time point 
mp_meta <- mp_meta_orig %>% 
  filter(body_site =="UBERON:feces") %>% 
  select(X.SampleID, age, collection_timestamp, host_subject_id) %>% 
  arrange(host_subject_id, age)
# count number of time points per subject 
table(mp_meta$host_subject_id)

# change column names in preparation for later application of functions 
colnames(mp_meta) <- c("sampID", "age", "collectiontime", "subjID")
mp_meta <- mp_meta %>% 
  mutate(sampID = paste0("X", as.character(sampID))) %>% 
  mutate(day = as.Date(as.character(mp_meta$collectiontime), '%m/%d/%y %H:%M')) %>% 
  group_by(subjID) %>% 
  mutate(firstday = min(day)) %>% 
  mutate(time = as.numeric(day - firstday)) %>% 
  arrange(subjID, time) 

# rows are OTUs  
# remove OTUs that don't appear in any subject 
mp_otu <- mp_otu_orig[, mp_meta$sampID]
all(colnames(mp_otu) == mp_meta$sampID)
if (any(apply(mp_otu, 1, FUN = function(x) sum(x != 0)) == 0)) {
  zotus <- which(apply(mp_otu, 1, FUN = function(x) sum(x != 0)) == 0)
  mp_otu <- mp_otu[-zotus, ]
}
# summarize read counts and unique OTUs 
summary(apply(mp_otu, 2, sum))
dim(mp_otu)

# save as package data 
usethis::use_data(mp_meta, overwrite = TRUE)
usethis::use_data(mp_otu, overwrite = TRUE)
usethis::use_data(mp_tree, overwrite = TRUE)
usethis::use_data(mp_tax, overwrite = TRUE)

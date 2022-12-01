## code to process Ravel vaginal microbiome dataset
library(tidyverse) 

filepath.metadata <- system.file("extdata", "ravel_metadata_AdditionalFile1.csv", package="MBVolDescrip")
filepath.otus <- system.file("extdata", "ravel_otu_props_AdditionalFile3.csv", package="MBVolDescrip")

## Metadata in Additional File 1 
## OTU data in Additional File 3 
## Article link: https://microbiomejournal.biomedcentral.com/articles/10.1186/2049-2618-1-29 

# Load metadata 
rv_meta_raw <- read.csv(filepath.metadata)

# Exclude all observations after first BV occurrence 
# Exclude subjects observed fewer than 20 times 
excl_date_df <- rv_meta_raw %>% 
  arrange(subjectID, dayInStudy) %>% 
  group_by(subjectID) %>% 
  summarize(excl_date = ifelse(any(SBV == 1 | ABV == 1), 
                               min(dayInStudy[SBV == 1 | ABV == 1]), 
                               NA))
rv_meta <- rv_meta_raw %>% 
  merge(., excl_date_df, by="subjectID", all=T) %>% 
  mutate(excl_flag = ifelse(excl_date <=  dayInStudy, 1, 0)) %>% 
  filter(excl_flag == 0) %>% 
  filter(! subjectID %in%  c(3, 5, 130)) ## observed fewer than 20 times

# Load OTUs: columns are OTUs, rows are samples 
rv_otus_raw <- read.csv(filepath.otus, row.names=1)

# Filter to retained samples 
rv_meta <- rv_meta %>% 
  filter(sampleID %in% rownames(rv_otus_raw))
rv_otus <- rv_otus_raw[rv_meta$sampleID, ]

# Summarize subject level information: number of subjects, time points per subject 
table(rv_meta$subjectID)

# Summarize read counts 
all(rownames(rv_otus) == rv_meta$sampleID)
rv_meta$ReadCount = rv_otus$Total.Reads
summary(rv_meta$ReadCount)

# Summarize OTU level information 
rv_otus$Total.Reads <- NULL 

# Convert proportions to counts 
rv_otu_counts <- rv_otus 
for (i in 1:nrow(rv_otu_counts)) {
  rv_otu_counts[i, ] <- round(rv_otus[i,] * rv_meta$ReadCount[i]/100)
}

# Remove never-present OTUs  
if (any(apply(rv_otu_counts, 2, sum) == 0)) {
  rv_otu_counts <- rv_otu_counts[, -which(apply(rv_otu_counts, 2, sum) == 0)]
}

# Number of unique OTUs 
dim(rv_otu_counts)

# Rename metadata columns 
rav_meta <- rv_meta %>% 
  dplyr::rename(sampID = sampleID, 
                subjID = subjectID,
                time = dayInStudy) %>% 
  arrange(subjID, time) %>% 
  select(-excl_date, -excl_flag)
rav_otu <- rv_otu_counts[rav_meta$sampID, ]

all(rownames(rav_otu) == rav_meta$sampID)

# save as package data 
usethis::use_data(rav_meta, overwrite = TRUE)
usethis::use_data(rav_otu, overwrite = TRUE)

#############
### Setup ###
#############

library(tidyverse)
library(MicrobeDS) # for Moving Pictures dataset 
library(biomformat) # for SMP biom file 
#remotes::install_github("microbiome/microbiomeDataSets")
library(microbiomeDataSets) # for Blekhman baboon data 

# Processed data will be saved in this directory 
# (Metadata and OTU) 
if (!dir.exists("./FinalDatasets/")) {
  dir.create("./FinalDatasets/")
}




##################################################
### Dataset 1: Caporaso 2011 - Moving Pictures ###
##################################################

## Data available in many locations, including MicrobeDS R package 

# load data (from MicrobeDS package if not yet loaded and saved; otherwise from file)
if (!file.exists("./Caporaso_MovingPictures/metadata.rds")) {
  data("MovingPictures")
  
  # load all components of Moving Pictures dataset 
  mp_meta <- data.frame(sample_data(MovingPictures))
  mp_tax <- data.frame(tax_table(MovingPictures))
  mp_otu <- data.frame(otu_table(MovingPictures))
  
  # save separate objects 
  saveRDS(mp_meta, file="./Caporaso_MovingPictures/metadata.rds")
  saveRDS(mp_tax, file="Caporaso_MovingPictures/tax.rds")
  saveRDS(mp_otu, file="Caporaso_MovingPictures/asv.rds")
  
} else {
  
  # load previously saved objects 
  mp_meta_orig <- readRDS("./Caporaso_MovingPictures/metadata.rds")
  mp_tax_orig <- readRDS("./Caporaso_MovingPictures/tax.rds")
  mp_otu_orig <- readRDS("./Caporaso_MovingPictures/asv.rds")
  
}

# filter to gut microbiome time series;  sort by subject ID and time point 
mp_meta <- mp_meta_orig %>% 
  filter(body_site =="UBERON:feces") %>% 
  select(X.SampleID, age, collection_timestamp, host_subject_id) %>% 
  arrange(host_subject_id, age)
# count number of time points per subject 
table(mp_meta$host_subject_id)

# rows are OTUs  
# remove OTUs that don't appear in any subject 
mp_otu <- mp_otu_orig[, paste0("X", mp_meta$X.SampleID)]
if (any(apply(mp_otu, 1, FUN = function(x) sum(x != 0)) == 0)) {
  zotus <- which(apply(mp_otu, 1, FUN = function(x) sum(x != 0)) == 0)
  mp_otu <- mp_otu[-zotus, ]
}
# summarize read counts and unique OTUs 
summary(apply(mp_otu, 2, sum))
dim(mp_otu)

# Save datasets 
saveRDS(mp_otu, "./FinalDatasets/mp_otu.rds")
saveRDS(mp_meta, "./FinalDatasets/mp_meta.rds")




###########################################################
### Dataset 2: Flores 2014 - Student Microbiome Project ###
###########################################################

## Data available on GitHub 
## https://github.com/caporaso-lab/student-microbiome-project

# Metadata 
# filter to gut time series and no antibiotic use 
smp_metadata <- read.csv("./Flores_SMP/smp_metadata.csv") %>% 
  filter(GutTimeseries == "Yes") %>% 
  filter(PersonalIDAntibioticDisturbance == "No") %>% 
  arrange(PersonalID, WeeksSinceStart) %>% 
  select(SampleID, PersonalID, WeeksSinceStart, SequenceCountAbove10000, 
         Control, PD, ObservedOTUs, ShannonEvenness, Age, Gender, 
         BMI, BMIclass, SampleSicknessDisturbance, 
         SampleMenstruationDisturbance, 
         PersonalIDSicknessDisturbance, 
         PersonalIDMenstruationDisturbance) 
# count subjects and time points per subject 
length(unique(smp_metadata$PersonalID))
table(smp_metadata$PersonalID)

# OTU matrix -- taxa are rows 
smp_otus_biom <- read_biom("./Flores_SMP/otu_table.ts_only.biom")
tax_key <- observation_metadata(smp_otus_biom)

# keep samples to match metadata 
smp_otus <- as(biom_data(smp_otus_biom), "matrix")
smp_otus <- smp_otus[, smp_metadata$SampleID]

# exclude OTUs with zero counts for all retained samples 
if (any(apply(smp_otus, 1, FUN = function(x) sum(x != 0)) == 0)) {
  zotus <- which(apply(smp_otus, 1, FUN = function(x) sum(x != 0)) == 0)
  smp_otus <- smp_otus[-zotus, ]
}

# summarize characteristics: unique OTUs, read depth 
summary(apply(smp_otus, 2, sum))
dim(smp_otus) 

# Save datasets 
saveRDS(smp_otus, "./FinalDatasets/smp_otu.rds")
saveRDS(smp_metadata, "./FinalDatasets/smp_meta.rds")




##################################################
### Dataset 3: Gajer Vaginal Microbiome (2012) ###
##################################################

## All data in Table S2 of Gajer article 
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3722878/
## doi: https://doi.org/10.1126/scitranslmed.3003605 

# Metadata 
# Interpretable values for previously-coded variables (e.g., race) 
gv_meta_raw <- read.csv("./Gajer_Vaginal/metadata.csv") 
gv_meta <- gv_meta_raw %>% 
  mutate(Race = case_when(Racea == 1 ~ "White", 
                          Racea == 0 ~ "Black", 
                          Racea == 5 ~ "Hispanic", 
                          Racea == 4 ~ "Others")) %>% 
  dplyr::rename(SampleID = Sample.ID, 
                StudyTime = Time.in.study, 
                SubjectID = Subject.ID, 
                NugentScore = Nugent.Score, 
                NugentCat = Nugent.Categoryb, 
                CST = Community.State.Typec, 
                ReadCount = Total.Read.Countsd)

# Load proportions (from article) and check match of sample IDs 
gv_otus_props <- read.csv("./Gajer_Vaginal/OTUs.csv", row.names = 1)
all(rownames(gv_otus_props) == gv_meta$Sample.ID) 

# Convert to counts -- rounded to whole numbers 
gv_otus_counts <- gv_otus_props 
for (i in 1:nrow(gv_otus_counts)) {
  gv_otus_counts[i, ] <- round(gv_otus_props[i, ] * gv_meta$ReadCount/100)
}

# Summary stats: number of samples per subject
table(table(gv_meta$SubjectID))
median(table(gv_meta$SubjectID))

# Summary stats: number of OTUs per sample 
# Taxa are columns 
summary(gv_meta$ReadCount)
dim(gv_otus_counts)

# Save datasets 
saveRDS(gv_otus_counts, "./FinalDatasets/gajer_otu.rds")
saveRDS(gv_meta, "./FinalDatasets/gajer_meta.rds")


##################################################
### Dataset 3: Ravel Vaginal Microbiome (2013) ###
##################################################

## Metadata in Additional File 1 
## OTU data in Additional File 3 
## Article link: https://microbiomejournal.biomedcentral.com/articles/10.1186/2049-2618-1-29 

# Load metadata 
rv_meta_raw <- read.csv("./Ravel_Vaginal/ravel_metadata.csv")

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
rv_otus_raw <- read.csv("./Ravel_Vaginal/ravel_otu_props.csv", 
                        row.names=1)
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

# Save datasets 
saveRDS(rv_otu_counts, "./FinalDatasets/ravel_otu.rds")
saveRDS(rv_meta, "./FinalDatasets/ravel_meta.rds")

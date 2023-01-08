## code to process Student Microbiome Project dataset 
library(tidyverse) 
library(biomformat) 

filepath.metadata <- system.file("extdata", "smp_metadata.csv", package="MBVolDescrip")
filepath.otus <- system.file("extdata", "smp_otu_table_gut.ts_only.biom.gz", package="MBVolDescrip")

## Read and process
smp_metadata <- read.csv(filepath.metadata) %>% 
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

# change column names for later function use 
smp_meta <- smp_metadata %>% 
  dplyr::rename(sampID = SampleID, 
                subjID = PersonalID) %>% 
  mutate(time = as.numeric(WeeksSinceStart)*7) %>% 
  arrange(subjID, time) ## time in days 

# OTU matrix -- taxa are rows 
smp_otus_biom <- read_biom(filepath.otus)

# keep samples to match metadata 
smp_otu_species <- as(biom_data(smp_otus_biom), "matrix")
smp_otu_species <- smp_otu_species[, smp_meta$sampID]

# merge OTUs at genus level 
tax_key <- observation_metadata(smp_otus_biom)
tax_long <- data.frame(taxID = rownames(tax_key), 
                       taxonomy = apply(tax_key[, 1:6], 1, FUN = function(x) paste0(x, collapse="_", sep="")))
tax_long <- tax_long[rownames(smp_otu_species), ]
all(tax_long$taxID == rownames(smp_otu_species))

new_tax_key <- data.frame(taxID = NA, 
                          taxon = aggregate(smp_otu_species[,1], 
                                            by=list(tax_long$taxonomy), 
                                            FUN = function(y) sum(y))$Group)
new_tax_key$taxID <- paste("OTU", 1:nrow(new_tax_key), sep="_")

smp_otu_genus <- apply(smp_otu_species, 2, FUN = function(x) {
  aggregate(x, by=list(tax_long$taxonomy), FUN = function(y) sum(y))$x
})
rownames(smp_otu_genus) <- new_tax_key$taxID

# exclude OTUs with zero counts for all retained samples 
if (any(apply(smp_otu_genus, 1, FUN = function(x) sum(x != 0)) == 0)) {
  zotus <- which(apply(smp_otu_genus, 1, FUN = function(x) sum(x != 0)) == 0)
  smp_otu_genus <- smp_otu_genus[-zotus, ]
}

# summarize characteristics: unique OTUs, read depth 
summary(apply(smp_otu_genus, 2, sum))
dim(smp_otu_genus) 

# Double check match between samples in OTU table and metadata
all(colnames(smp_otu_genus) == smp_meta$sampID)



# save as package data 
smp_otu <- smp_otu_genus 
smp_tax_key <- new_tax_key 

usethis::use_data(smp_meta, overwrite = TRUE)
usethis::use_data(smp_otu, overwrite = TRUE)
usethis::use_data(smp_tax_key, overwrite = TRUE)

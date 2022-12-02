# Describing Microbiome Volatility 

This project aims to describe how measures of microbiome volatility depend on factors such as density of sampling (and relatedly, balanced or unbalanced data) and rarefaction. 

The key functions, raw datasets, data processing code, and final datasets are available in the R package MBVolDescrip. You can install this package via: 

```
library(devtools)
devtools::install_github("aplantin/describing-volatility", subdir="MBVolDescrip")
```

### Datasets 

Throughout these investigations, we use four publicly available microbiome time series datasets: 

1. Moving Pictures, Caporaso et al. (2011) 
    - gut microbiome 
    - 2 subjects
    - 131 and 336 daily samples 
2. Student Microbiome Project, Flores et al. (2014) 
    - gut microbiome
    - 58 subjects (with no antibiotic treatment during the study)
    - weekly samples over the course of 3 months; 7-10 time samples per subject
3. Healthy vaginal microbiome, Gajer et al. (2012) 
    - vaginal microbiome 
    - 32 subjects 
    - 25-33 samples per subject, twice-weekly 
4. Pre-bacterial vaginosis microbiome, Ravel et al. (2013) 
    - vaginal microbiome 
    - 6 subjects (with at least 20 daily samples prior to BV diagnosis) 
    - 23-38 daily samples 

Processed metadata and OTU or ASV tables are available in the R package: 

```
# Moving Pictures
data(mp_otu); data(mp_meta); data(mp_tree); data(mp_tax)

# Student Microbiome Project 
data(smp_otu); data(smp_meta); data(smp_tax_key) 

# Gajer vaginal microbiome 
data(gaj_otu); data(gaj_meta) 

# Ravel vaginal microbiome
data(rav_otu); data(rav_meta) 
```


### Function Usage 

The package provides three main functions. 

- `temporalSubsampleMeta()` identifies pairs of samples separated by the desired time spacing. 
- `calcMicrobiomeChanges()` calculates the additive, multiplicative, and qualitative (presence/absence) between these pairs of time points, as well as extracting the relevant distances from a distance matrix if provided. 
- `summMicrobiomeVolatility()` calculates several ad hoc measures of microbiome volatility. 
    
    
```
data(gaj_otu) 
data(gaj_meta) 

# Rarefying so that total read count is the same across time points and subjects 
min_readct <- min(apply(gaj_otu, 1, sum))
gaj_otu_rarefy <- rrarefy(gaj_otu, sample=min_readct)

# Temporal change metadata (one week spacing with one day wiggle room: samples are 6-8 days apart) 
gaj_changemeta_7 <- temporalSubsampleMeta(gaj_otu_rarefy, gaj_meta, desired_spacing = 7, 
                                          window_width = 1, taxaAreRows = FALSE) 

# Calculate distance matrix on full dataset for input to change
gaj_D_BC <- as.matrix(vegdist(gaj_otu_rarefy, method="bray"))

# Calculate microbiome changes at desired time gaps 
gaj_changematr_7 <- calcMicrobiomeChanges(gaj_otu_rarefy, gaj_changemeta_7, 
                                          Ds = gaj_D_BC, taxaAreRows = FALSE)

# Calculate volatility summaries 
gaj_vol_7 <- summMicrobiomeVolatility(mbchanges = gaj_changematr_7)
```

### Analyses 

The code to replicate the analyses in the manuscript is available in the R/ directory. 

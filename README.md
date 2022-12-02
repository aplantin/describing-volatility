# Describing Microbiome Volatility 

This project aims to describe how measures of microbiome volatility depend on factors such as density of sampling (and relatedly, balanced or unbalanced data) and rarefaction. 

The key functions, raw datasets, data processing code, and final datasets are available in the R package MBVolDescrip. You can istall this package via: 

```
library(devtools)
devtools::install_github("aplantin/describing-volatility", subdir="MBVolDescrip")
```

### Datasets 

To do this, we use four publicly available microbiome time series datasets: 

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

- Moving Pictures Gut Microbiome: mp_otu, mp_meta, mp_tree, mp_tax
- Student Microbiome Project Gut Microbiome: smp_otu, smp_meta, smp_tax_key 
- Gajer Vaginal Microbiome: gaj_otu, gaj_meta 
- Ravel Vaginal Microbiome: rav_otu, rav_meta 


### Analyses 


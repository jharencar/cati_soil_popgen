# Code and data for "Extremely fine-scale soil heterogeneity in a rare serpentine endemic plant shapes patterns of genetic diversity"
This repository contains scripts and input files used in the analysis for the manuscript titled, "Extremely fine-scale soil heterogeneity in a rare serpentine endemic plant shapes patterns of genetic diversity", which is on [biorxiv here](https://www.biorxiv.org/content/10.1101/2025.09.01.673272v1.full.pdf) and will be submitted to Evolution and Ecology shortly. 

## Contents:
### code:
`GDM_script.r` - code for the Generalized Dissimilarity Matrix analysis \
`Sample_map.qmd` - code for generating Figure 1, the sample map with example soil chemistry plots \
`mantel_tests.R` - code for running the Mantel tests \
`soil_covariance.qmd` - code for generating soil covariances and covariance matrix plots \

### data:
`CATIgind.Rdata` - genind file with RADseq data for the analyses saved as a `.Rdata` file \
`lat_long_r.csv` - coordinates of the sample collection locations for generating the sample map \
`soils_r.csv` - soil chemistry data for all sample collection locations (THIS is the one I use in `soil_covariance.qmd`; check with Joe about `soils_complete.csv`) \
`soils_complete.csv` - soil chemistry data for all sample collection locations \

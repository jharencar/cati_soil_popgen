#### CaTi genomics scripts #####

# Overview #
# This script tests for correlations between SNPs identified in a ddRAD
# Calochortus tiburonensis sequqncing run and patterns of soil chemistry and
# composition using individual mantel tests for each soil compound/measurement

# Load dependencies
## install.packages("vegan")
## install.packages("adegenet")
## install.packages("poppr")

library(vegan)
library(adegenet)
library(poppr)
library(distances)
library(readr)


# Plants are listed in genid object in numerical order based on their ID
load("C:/Users/PlantagoMacine/Documents/GitHub/cati_soil_popgen/CATIgind.Rdata")

# Load soils data. Convert to data frame because I am old
soils_r <- read_csv("soils_r.csv")
soils <- as.data.frame(soils_r)

# Sort soil data by plant ID to correspond with genetic data
soils <- soils[order(soils$plantID), ]

# trim empty rows
soils <- soils[1:24,]

# Produce matrix of genetic distances for CaTi
# We use euclidian distances because many alleles are NA for many individuals
# most estimates of genetic distance are also designed to compare populations and
# not individuals
eucl_dist <- dist(CATIgind$tab)


# Produce matrix of differnces in the abundance of soil components (Ba, Co, Mg, 
# etc. ), then perform mantel test for correlations with each and the matrix of 
# genetic distance. Output the r and p values for each comparison. This will be 
# done with a for loop that itterates through each chemical compound.


# Initialize empty vectors to store information from the loop
soil_compound <- c()
r_value <- c()
p_value <- c()

# The first three columns are not soil chemical data, so we start the loop with
# columm four

for (i in 4:ncol(soils)){
  
  # isolate one soil component
  this_compound <- soils[,i]
  
  # Produce distance matrix for the soil component
  chem_matrix <- dist(this_compound)
  
  # Perform mantel test
  mantel_out <- mantel(chem_matrix, eucl_dist,method = "pearson",
    permutations = 999,  strata = NULL, na.rm = FALSE)
  
  # Save information from the mantel test
  soil_compound <- c(soil_compound, colnames(soils_r[,i]))
  r_value <- c(r_value, mantel_out$statistic)
  p_value <- c(p_value, mantel_out$signif)
  
  
}

# Assemble dataframe of results
Mantel_results <- data.frame(soil_compound, r_value, p_value)








####### Analyses with non-redundant soil data ##############
# 14 of the plants sampled were within 1-5cm of each other. As a result, they
# occupied the same soil environment. Here we removed the second sampled plant
# and re-tested the relationships

# Create list of plant numbers to remove
to_remove <- c(2, 5, 8, 9, 14, 17, 23)

# remove double soil plants 
soils_mod <- soils[-which(soils$plantID %in% to_remove), ]

# remove plants from table of genetic data
cati_alleles <- CATIgind$tab
alleles_mod <- cati_alleles[-to_remove,]

## Produce distance matrix for modified alleles table
gen_dist_mod <- dist(alleles_mod)


# Initialize empty vectors to store information from the loop
soil_compound <- c()
r_value_mod <- c()
p_value_mod <- c()

# The first three columns are not soil chemical data, so we start the loop with
# columm four

for (i in 4:ncol(soils_mod)){
  
  # isolate one soil component
  this_compound <- soils_mod[,i]
  
  # Produce distance matrix for the soil component
  chem_matrix <- dist(this_compound)
  
  # Perform mantel test
  mantel_out <- mantel(chem_matrix, gen_dist_mod, method = "pearson",
                       permutations = 999,  strata = NULL, na.rm = FALSE)
  
  # Save information from the mantel test
  soil_compound <- c(soil_compound, colnames(soils_r[,i]))
  r_value_mod <- c(r_value_mod, mantel_out$statistic)
  p_value_mod <- c(p_value_mod, mantel_out$signif)
  
  
}

# Assemble dataframe of results
Mantel_results_mod <- data.frame(soil_compound, r_value_mod, p_value_mod)





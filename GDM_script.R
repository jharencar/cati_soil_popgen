#### CaTi genomics scripts #####

# Overview #
# This script tests for correlations between SNPs identified in a ddRAD
# Calochortus tiburonensis sequqncing run and patterns of soil chemistry and
# composition using generalized dissimilarity models


# Read file

# Load dependencies
## install.packages("vegan")
## install.packages("adegenet")
## install.packages("poppr")
## install.packages("gdm")

library(vegan)
library(adegenet)
library(poppr)
library(distances)
library(readr)
library(gdm)

# Plants are listed in genid object in numerical order based on their ID
load("C:/Users/PlantagoMacine/Documents/GitHub/cati_soil_popgen/CATIgind.Rdata")

# Load soils data. Convert to data frame because I am old
soils_r <- read_csv("soils_r.csv")
soils <- as.data.frame(soils_r)

# Sort soil data by plant ID to correspond with genetic data
soils <- soils[order(soils$plantID), ]

# trim empty rows
soils <- soils[1:24,]

# Load spatial coordinates
lat_long <- read_csv("lat_long_r.csv")
lat_long <- as.data.frame(lat_long)


# Combine lat long and soil data
descriptors <- cbind(lat_long, soils)

# Trim to only lat/long and soil variables. This will be the environmental data 
# for the GDM
descriptors_mod <- descriptors[c(2,3,4, 10:ncol(descriptors))]

# Remove plant 9, which does not have a corresponding soil sample
descriptors_mod <- descriptors_mod[-9,]

# Get euclidian estimates of genetic distance
cati_alleles <- CATIgind$tab

# Remove plant 9 from the data frame of alleles also
cati_alleles_mod <- cati_alleles[-9,]



# calculate euclidian genetic distance. This approach is used to create the
# square difference matrix and is capable of removing NAs from the calculation
# of distance
eucl_dist <- as.matrix( vegdist(cati_alleles_mod, method="euclidean", diag = TRUE, upper = TRUE, na.rm = TRUE))



## Distance values need to be scaled between 0 and 1. Divide all values in matrix by max value
eucl_dist <- as.matrix(eucl_dist/max(eucl_dist))


# Attach the plant IDs to the data frame of genetic distance
eucl_dist <- cbind(descriptors_mod_2$plantID, as.data.frame(eucl_dist))
colnames(eucl_dist)[1] <- c("plantID")




#Peform GDM
gdmTab.dis <- formatsitepair(bioData=eucl_dist, 
                             bioFormat=3, #diss matrix 
                             XColumn="Long_dd", 
                             YColumn="Lat_dd", 
                             predData=descriptors_mod, 
                             siteColumn="plantID")

gdm.1 <- gdm(data=gdmTab.dis, geo=TRUE)



gdm.1.splineDat <- isplineExtract(gdm.1)

plot(gdm.1.splineDat$x[,"Barium"], 
     gdm.1.splineDat$y[,"Barium"], 
     lwd=3,
     type="l", 
     xlab="Barium Content", 
     ylab="Partial genetic distance")
#### CaTi genomics scripts #####

# Overview #
# This script uses an RDA to test for correlations between genetic variation in
# a ddRAD Calochortus tiburonensis sequqncing run and soil components, notably
# Mg, Ca, and Ni. Additionaly, this script will perform a PCA of soil variation
# and then test for correlations between genetic variation and soil PC axes



# Load Dependencies
## install.packages("vegan")
## install.packages("adegenet")
## install.packages("viridis")
## install.packages("ggfortify")


library(vegan)
library(adegenet)
library(viridis)
library(ggfortify)


# Read in files
# Plants are listed in genid object in numerical order based on their ID
load("C:/Users/PlantagoMacine/Documents/GitHub/cati_soil_popgen/CATIgind.Rdata")

# Load soils data. Convert to data frame because I am old
soils_r <- read.csv("soils_r.csv")
soils <- as.data.frame(soils_r)

# Sort soil data by plant ID to correspond with genetic data
soils <- soils[order(soils$plantID), ]

# trim empty rows
soils <- soils[1:24,]

# Load spatial coordinates
lat_long <- read.csv("lat_long_r.csv")
lat_long <- as.data.frame(lat_long)


# Combine lat long and soil data
descriptors <- cbind(lat_long, soils)

# Remove last column from soil table, because it is Mg/Ca ratio and redundant
# with columns for Mg and Ca

descriptors <- descriptors[,-ncol(descriptors)]

# Trim to only lat/long and soil variables. This will be the environmental data 
# for the GDM
descriptors_mod <- descriptors[c(2,3,4,5, 10:ncol(descriptors))]

# Remove plant 9, which does not have a corresponding soil sample
descriptors_mod <- descriptors_mod[-9,]

# Get euclidian estimates of genetic distance
cati_alleles <- CATIgind$tab

# Remove plant 9 from the data frame of alleles also
cati_alleles_mod <- cati_alleles[-9,]




#### Create DF only with loci called in all individuals ###
#cati_no_na <- cati_alleles_mod[ , colSums(is.na(cati_alleles_mod))==0]




# RDA cannot be performed with NAs in the response variable. We thus need to 
# replace the NAs in our dataset somehow. Here we impute the most common allele 
# as a replacement for NAs. OTHER METHODS COULD BE MORE APPROPIRATE and it is 
# worth testing how the results change based on how NAs are replaced

# Replace NA with the most common allele for each locus
cati_allels_imp <- apply(cati_alleles_mod , 2, function(x) replace(x, is.na(x),
                                      as.numeric(names(which.max(table(x))))))


# Need to rename rows in allele data frame for graphing later
# split up plantIDs by 0
names1 <- unlist(strsplit(rownames(cati_allels_imp), "T"))
names1 <- names1[seq(2, length(names1), 2)]
# split up by "_"
names2 <- unlist(strsplit(names1, split= '_'))
names2 <- names2[seq(1, length(names2), 2)]

# replace CATI genetic data rownames
rownames(cati_allels_imp) <- names2




### Perform RDA ###
cati_dbrda <- dbrda(cati_allels_imp ~ Magnesium + Nickel + Calcium,
                data = descriptors_mod,
                distance = "euclidean",
                scale = FALSE,
                na.action = "na.exclude"
)
anova.cca(cati_dbrda)
anova.cca(cati_dbrda, by="axis")

# Full model is not significant, but dbRDA axis 1 is...
# full - F = 1.0639, P = 0.127
# RDA axis1: F = 1.2072, P = 0.006
# RDA axis1: F = 1.0934, P = 0.573
# RDA axis1: F = 0.8911, P = 0.781

#check dbRDA output
dbrda_summary <- summary(cati_dbrda)
dbrda_summary$cont

# dbRDA axis 1 explains 5.44% of the total variation, similar to the first
# unconstrained axis1 which explains 6.15% 

# Calculate adjusted r squared for the RDA
RsquareAdj(cati_dbrda)

dbrda_summary$biplot
# the first dbRDA axis is most strongly correlated with Mg (-0.9874), 
# Ni = 0.3497, Ca = 0.3334



#Plot RDA 
#set up colors for 6 populations with viridis package
levels(descriptors_mod$population) <- c("FernGulch", "SanQuentin","PetroglyphRock", "Calamagrostis","Taylor", "WestwardRidge")
pops <- descriptors_mod$population
v_cols <- viridis(6)
  

# plot RDA output
plot(cati_dbrda, type ="n", scaling = 3)
text(cati_dbrda, scaling=3, display="bp", col="#0868ac", cex=1)  
points(cati_dbrda, display="sites", pch=21, cex=1.3,
       col="gray32", scaling=3, bg = v_cols)
legend("bottomleft", legend=levels(pops), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)



#######################################################################
#### RDA with PCA of soil chemistry ####
######################################################################

# It is possible that we can account for total soil differences in the RDA 
# framework by using the princpal component loadings of soil variation

# To perform a PCA of soil variation, we must use the entire soils dataset. 
# Plant loadings can still be extracted from this more complete representation
# of soil variation.

# Load full soils data
soils_complete <- read.csv("soils_complete.csv")

# Load soils data. Convert to data frame because I am old
soils_r <- read.csv("soils_r.csv")
soils <- as.data.frame(soils_r)

# Sort soil data by plant ID to correspond with genetic data
soils <- soils[order(soils$plantID), ]

# Remove first three columns of soils data frame, which contain plant ID, 
# patch name, plot id, moisture, pH, and organic matter
soils_complete_mod <- as.data.frame(soils_complete[,-c(1:5, ncol(soils_complete))])


#retain soil plot IDs
plot_IDs <- as.data.frame(soils_complete[,2])

# Perform the PCA
soils_pca <- prcomp(soils_complete_mod, scale = TRUE)

# Extract points for the first four PC axes. The first four axes explain 
# 72.6% of the variation (0.329, 0.1869, 0.1172, 0.09285). Combine this data 
# plot IDs

pca_values <- data.frame(plot_IDs, soils_pca$x[,1:4])
colnames(pca_values)[1] <- c("plotID")


# Make a new dataframe by pulling PCA values via patching plant IDs in the 'soils'
# object to the PCA values via "Plot ID"
pca_descriptors <- merge(soils[,c(1:3)], pca_values, by.x = c("plotID"))

# remove redundant assignments and plant 9
pca_descriptors <- pca_descriptors[ -c(22, 24, 25), ]

# order the PCA table to match the genetic data
pca_descriptors <- pca_descriptors[order(pca_descriptors$plantID), ]


### Perform dbRDA ###

cati_rda_pca <- dbrda(cati_allels_imp ~ PC1 + PC2 + PC3,
                data = pca_descriptors,
                scale = FALSE,
                na.action = "na.exclude",
                subset = NULL,
)


# Calculate adjusted r squared for the RDA
RsquareAdj(cati_rda_pca)

# Use an anova to test for significant effects of predictors in RDA
anova.cca(cati_rda_pca)

# PCA model is significant F = 1.1091, P = 0.004
# Test for significance of factors
# Need to remove some large objects to allocate memory




rm(cati_alleles, cati_alleles_mod, CATIgind)
memory.limit(size=56000)

axis_out <- anova.cca(cati_rda_pca, by="axis")

# The first RDA axis was significant (F = 1.2042, P = 0.016, similar to the 
# full model). Axes 2 and 3 were not significant (F = 1.1198, p = 0.349; 
# F = 1.0033, P = 0.540 respectively)


#Plot RDA 
#set up colors for 6 populations with viridis package
levels(pca_descriptors$population) <- c("FernGulch","SanQuentin","PetroglyphRock",
                                        "Calamagrostis","Taylor", "WestwardRidge")
pops <- pca_descriptors$population
v_cols <- viridis(6)


# plot RDA output
plot(cati_rda_pca, type ="n", scaling = 3)
text(cati_rda_pca, scaling=3, display="bp", col="#0868ac", cex=1)  
points(cati_rda_pca, display="sites", pch=21, cex=1.3,
       col="gray32", scaling=3, bg = v_cols)
legend("bottomright", legend=levels(pops), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)




########## PCA of genetic diversity ###########
# The PCA is here because it requires the interpolated data from 

cati_alleles <- CATIgind$tab
cati_allels_imp <- apply(cati_alleles , 2, function(x) replace(x, is.na(x),
                   as.numeric(names(which.max(table(x))))))


snp_pca <- prcomp(cati_allels_imp )

plantID <- c(1:24)
cati_data <- cbind(plantID , cati_allels_imp)

#plot PCA
### Need to associate population assignments with PCA output
#Start with the modified file from soil PCA
to_mod <- pca_descriptors[,c(2:3)]
to_mod[24,] <- c(0, "ex") 
to_mod[c(10:24),] <- to_mod[c(9:23),]
to_mod[9,] <- c(9, "PetroglyphRock")

cati_gen_pops <- cbind(to_mod, cati_allels_imp)



# Make dataframe that has population and color assignments
v_cols <- cbind(viridis(6), c("FernGulch","SanQuentin","PetroglyphRock",
                            "Calamagrostis","Taylor", "WestwardRidge"))
colnames(v_cols) <- c("color", "population")

# Attach colors to genetic dataframe
cati_gen_colors <- merge(cati_gen_pops, v_cols, by.x = c("population"))


p <- autoplot(snp_pca, data = cati_gen_colors, colour = 'colors')


ggplotly(p)
#### CaTi genomics scripts #####

# Overview #
# This script uses NMDS to test for clustering in soil structure across
# Calochortus tiburonensis sampling sites
# 


# Load dependencies
## install.packages("vegan")
## install.packages("adegenet")
## install.packages("poppr")
## install.packages("gdm")

library(vegan)
library(ggplot2)
library(viridis)

# Load soils data. Convert to data frame because I am old
soils_r <- read.csv("soils_r.csv")
soils <- as.data.frame(soils_r)


# Remove row 18, which is a double measurement for two nearby plants
soils <- soils[-18,]

# Save patch labels
patch <- soils [,1:2]



# Remove all identifiers, redundant measurements (Ca/Mg), and dubious measurements (soil moisture, pH)
soil_comp <- soils [,6:30]


######## ORDINATION START######################
# Chose random number to set seed
#runif(1, 1, 1000) #139
set.seed(139)

#Run NMDS
soil.ord <- metaMDS(soil_comp , dist = "bray", wascores = TRUE, trymax = 999)
#stress = 0.15


# Create NMDS plot
# set up colors
levels(patch$population ) <- c("FernGulch", "SanQuentin","PetroglyphRock", "Calamagrostis","Taylor", "WestwardRidge")
groups <- patch$population
v_cols <- viridis(6)


#plot NMDS output
mds.fig.soil <- ordiplot(soil.ord, type = "none",  xlim=c(-0.08, 0.07), ylim=c(-0.08,0.06))
points(mds.fig.soil, "sites", pch = 19, col = v_cols[1], select = patch$population == c("FernGulch"))
points(mds.fig.soil, "sites", pch = 19, col = v_cols[2], select = patch$population == c("SanQuentin"))
points(mds.fig.soil, "sites", pch = 19, col = v_cols[3], select = patch$population == c("PetroglyphRock")) 
points(mds.fig.soil, "sites", pch = 19, col = v_cols[4], select = patch$population == c("Calamagrostis")) 
points(mds.fig.soil, "sites", pch = 19, col = v_cols[5], select = patch$population == c("Taylor")) 
points(mds.fig.soil, "sites", pch = 19, col = v_cols[6], select = patch$population == c("WestwardRidge"))



############### NMDS for plant fitness compounds only ####################3
soil_tox <- soil_comp [,c(6,10,12,13,15,16,17,20)]

set.seed(140)

#Run NMDS
soil.tox.ord <- metaMDS(soil_tox , dist = "bray", wascores = TRUE, trymax = 999)

#plot NMDS output
mds.fig.tox <- ordiplot(soil.tox.ord, type = "none",  xlim=c(-0.06, 0.11), ylim=c(-0.06,0.09))
points(mds.fig.tox, "sites", pch = 19, col = v_cols[1], select = patch$population == c("FernGulch"))
points(mds.fig.tox, "sites", pch = 19, col = v_cols[2], select = patch$population == c("SanQuentin"))
points(mds.fig.tox, "sites", pch = 19, col = v_cols[3], select = patch$population == c("PetroglyphRock")) 
points(mds.fig.tox, "sites", pch = 19, col = v_cols[4], select = patch$population == c("Calamagrostis")) 
points(mds.fig.tox, "sites", pch = 19, col = v_cols[5], select = patch$population == c("Taylor")) 
points(mds.fig.tox, "sites", pch = 19, col = v_cols[6], select = patch$population == c("WestwardRidge"))
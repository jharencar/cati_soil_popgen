#### CaTi genomics scripts #####

# Overview #
# This script tests for correlations between SNPs identified in a ddRAD
# Calochortus tiburonensis sequqncing run and patterns of soil chemistry and
# composition using generalized dissimilarity models


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
library(ggplot2)

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

# Remove last column from soil table, because it is Mg/Ca ratio and redundant
# with columns for Mg and Ca

descriptors <- descriptors[,-ncol(descriptors)]

# Trim to only lat/long and soil variables. This will be the environmental data 
# for the GDM
descriptors_mod <- descriptors[c(3,4,5, 10:ncol(descriptors))]

# Remove plant 9, which does not have a corresponding soil sample
descriptors_mod <- descriptors_mod[-9,]

# Get euclidean estimates of genetic distance
cati_alleles <- CATIgind$tab

# Remove plant 9 from the data frame of alleles also
cati_alleles_mod <- cati_alleles[-9,]



# calculate euclidian genetic distance. This approach is used to create the
# square difference matrix and is capable of removing NAs from the calculation
# of distance
eucl_dist <- as.matrix( vegdist(cati_alleles_mod, method="euclidean", diag = FALSE, upper = FALSE, na.rm = TRUE))



## Distance values need to be scaled between 0 and 1. Divide all values in matrix by max value
eucl_dist <- as.matrix(eucl_dist/max(eucl_dist))


# Attach the plant IDs to the data frame of genetic distance
eucl_dist <- cbind(descriptors_mod$plantID, as.data.frame(eucl_dist))
colnames(eucl_dist)[1] <- c("plantID")



set.seed(200)
#Peform GDM
gdmTab.dis <- formatsitepair(bioData=eucl_dist, 
                             bioFormat=3, #diss matrix 
                             XColumn="Long_dd", 
                             YColumn="Lat_dd", 
                             predData=descriptors_mod, 
                             siteColumn="plantID")

gdm.1 <- gdm(data=gdmTab.dis, geo=TRUE)



gdm.1.splineDat <- isplineExtract(gdm.1)

plot(gdm.1.splineDat$x[,"Nickel"], 
     gdm.1.splineDat$y[,"Nickel"], 
     lwd=3,
     type="l", 
     xlab="Nickel Abundance (mg/Kg)", 
     ylab="Partial genetic turnover")

## Plot relative contributions of all variables
# extract isocline y max values

isocline_y <- gdm.1.splineDat$y
cont_vector <- isocline_y[nrow(isocline_y),]
cont_matrix <- as.data.frame(cont_vector)
max_isoclines <- data.frame(rownames(cont_matrix), cont_matrix$cont_vector)
colnames(max_isoclines) <- c("variable", "contribution")
max_isoclines <- max_isoclines[order(-max_isoclines$contribution),]
max_isoclines$contribution <- max_isoclines$contribution + 0.0001


ggplot( data = max_isoclines,aes( x = factor(variable, 
            level = c(max_isoclines[order(max_isoclines$contribution),1])),
            y = contribution)) +
  geom_bar(stat = "identity", aes(fill = contribution))+
  coord_flip()+
  theme_gray()+
  theme(legend.position="none")+
  labs(x="Soil Element", y = "Relative Contribution")


# Test for significance 
modTest <- gdm.varImp(gdmTab.dis, geo=T, nPerm=500, parallel=T, cores=10, predSelect=F)
# No individual soil component was significant, but the entire model was
# Model assessment`
#All predictors
#Model deviance                     17.605
#Percent deviance explained         41.869
#Model p-value                       0.002
#Fitted permutations               500.000



######## Loop for multiple runs of GDM


for (i in c(1:500)){
  
  #Peform GDM
  gdmTab.dis <- formatsitepair(bioData=eucl_dist, 
                               bioFormat=3, #diss matrix 
                               XColumn="Long_dd", 
                               YColumn="Lat_dd", 
                               predData=descriptors_mod, 
                               siteColumn="plantID")
  
  gdm.1 <- gdm(data=gdmTab.dis, geo=TRUE)
  
  
  
  gdm.1.splineDat <- isplineExtract(gdm.1)
  
  plot(gdm.1.splineDat$x[,"Sodium"], 
       gdm.1.splineDat$y[,"Sodium"], 
       lwd=3,
       type="l", 
       xlab="Sodium Abundance", 
       ylab="Partial genetic turnover")
  
  ## Plot relative contributions of all variables
  # extract isocline y max values

if( i == 1){  
    
  isocline_y <- gdm.1.splineDat$y
  cont_vector <- isocline_y[nrow(isocline_y),]
  cont_matrix <- as.data.frame(cont_vector)
  max_isoclines <- data.frame(rownames(cont_matrix), cont_matrix$cont_vector)
  colnames(max_isoclines) <- c("variable", "contribution")

}
  else {
    isocline_y <- gdm.1.splineDat$y
    cont_vector <- isocline_y[nrow(isocline_y),]
    cont_matrix <- as.data.frame(cont_vector)
    max_isoclines_2 <- data.frame(rownames(cont_matrix), cont_matrix$cont_vector)
    
    max_isoclines[ , ncol(max_isoclines) +1] <- max_isoclines_2[,2]
    colnames(max_isoclines)[ncol(max_isoclines)] <- paste0("new", i)
    
    print(i)
  }
  
}

# Check cumulative values for each element
test <- rowSums(max_isoclines[,-1])
out <- cbind.data.frame(max_isoclines$variable ,test )

## Many elements were not significant
## Remove them from the data set and re-run

descriptors_mod <- descriptors_mod[, -which(names(descriptors_mod) %in% 
                          c("Aluminum", "Arsenic", "Boron", "Calcium",
                            "Cadmium", "Chromium", "Potassium", "Manganese",
                            "Lead", "Sulfur", "Selenium", "Silicon", "NH4-N",
                             "Strontium"))]




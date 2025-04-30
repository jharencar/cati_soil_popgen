### These are scripts to perform DAPC with the CATI data

library(adegenet)
library(factoextra)


## load file
load("C:/Users/PlantagoMacine/Documents/GitHub/cati_soil_popgen/CATIgind.Rdata")


#### Create DF only with loci called in all individuals ####
cati_alleles <- CATIgind$tab


# PCA cannot be performed with NAs in the response variable. We thus need to 
# replace the NAs in our dataset somehow. Here we impute the most common allele 
# as a replacement for NAs. OTHER METHODS COULD BE MORE APPROPIRATE and it is 
# worth testing how the results change based on how NAs are replaced

# Replace NA with the most common allele for each locus
cati_allels_imp <- apply(cati_alleles, 2, function(x) replace(x, is.na(x),
                                                                   as.numeric(names(which.max(table(x))))))



#### Find best number of PCs and Clusters ####

Cluster.out <- find.clusters(CATIgind, max.n.clust=15) # Chose 10 as an estimate above the expected number of clusters based on sampling
  
  # All PC axes retain important information, keep 23). 
  # Min BIC is 1 (BIC = 233)



#Transpose matrix so individuals are columns
gen.curated <- t(cati_allels_imp)

# Do the PCA
PCA.out <- prcomp(gen.curated)

#Plot eigenvector values
fviz_eig(pca1, n=23)


#19 eigenvectors explain over 90% of the variance
sum(get_eig(PCA.out)[1:19,2])







#Perform PCA
pca1 <- dudi.pca(cati_allels_imp)


#Select colors for plots
colors <- brewer.pal(length(unique(pop(CATIgind))), name = "Dark2")
colors <- viridis(6)

levels(CATIgind$pop) <- c("Fern Gulch","San Quentin","Petroglyph Rock",
                                        "Calamagrostis","Taylor", "Westward Ridge")
pops <- CATIgind$pop



#Plot PC1 and PC2 no labels for paper
fviz_pca_ind(pca1, axes = c(1, 2),
             mean.point = FALSE,
             col.ind = pops, 
             geom = c("point"), pointsize = 5,alpha.ind = 0.6,
             palette = colors,
             legend.title = "Site Locations",
             repel = TRUE,
             ggtheme =theme_classic())+
  scale_shape_manual(values=c(19,19,19,19,19,19))



#Plot PC1 and PC2 with labels for interpretation
fviz_pca_ind(pca1, axes = c(1, 2),
             mean.point = FALSE,
             col.ind = pops, 
             geom = c("point", "text"), pointsize = 5,alpha.ind = 0.6,
             palette = colors,
             legend.title = "Site Locations",
             repel = TRUE,
             ggtheme =theme_classic())+
  scale_shape_manual(values=c(19,19,19,19,19,19))






#Plot PC3 and PC4 no labels for paper
fviz_pca_ind(pca1, axes = c(3, 4),
             mean.point = FALSE,
             col.ind = pops, 
             geom = c("point"), pointsize = 5,alpha.ind = 0.6,
             palette = colors,
             legend.title = "Site Locations",
             repel = TRUE,
             ggtheme =theme_classic())+
  scale_shape_manual(values=c(19,19,19,19,19,19))



#Plot PC3 and PC4 with labels for interpretation
fviz_pca_ind(pca1, axes = c(3, 4),
             mean.point = FALSE,
             col.ind = pops, 
             geom = c("point", "text"), pointsize = 5,alpha.ind = 0.6,
             palette = colors,
             legend.title = "Site Locations",
             repel = TRUE,
             ggtheme =theme_classic())+
  scale_shape_manual(values=c(19,19,19,19,19,19))





######################################################
### PCAs with no NAs
######################################################




#Transpose matrix so individuals are columns
No_NA_t <- t(cati_no_na)

# Do the PCA
PCA.out <- prcomp(No_NA_t)

#Plot eigenvector values
fviz_eig(PCA.out, n=23)


#19 eigenvectors explain over 90% of the variance
sum(get_eig(PCA.out)[1:19,2])


fviz_pca_ind(PCA.out,  
             col.ind = CADI.str$pop, 
             #           geom = "point",
             palette = brewer.pal(length(unique(pop(CADI.str))), name = "Dark2"),
             legend.title = "Site Locations",
             repel = TRUE)+
  scale_shape_manual(values=c(19,19,19,19,19,19))


#graph PCA, no labels, axes 1 and 2
fviz(PCA.out, element="ind", geom= "point",  pointsize = 4, alpha = 0.65,
     habillage=CADI.str$pop,
     palette = brewer.pal(length(unique(pop(CADI.str))), name = "Dark2"), invisible="quali")+
  scale_shape_manual(values=c(19,19,19,19,19,19))



#graph PCA, with labels, axes 1 and 2
fviz(PCA.out, element="ind", geom= c("point", "text"),  pointsize = 2, alpha = 0.65,
     habillage=CADI.str$pop, repel = TRUE, 
     palette = brewer.pal(length(unique(pop(CADI.str))), name = "Dark2"), invisible="quali")+
  scale_shape_manual(values=c(19,19,19,19,19,19))




#graph PCA, no labels, axes 1 and 2
fviz(PCA.out, element="ind", geom= "point",  pointsize = 4, alpha = 0.65,
     habillage=CADI.str$pop, axes = c(3, 4),
     palette = brewer.pal(length(unique(pop(CADI.str))), name = "Dark2"), invisible="quali")+
  scale_shape_manual(values=c(19,19,19,19,19,19))



#graph PCA, with labels, axes 1 and 2
fviz(PCA.out, element="ind", geom= c("point", "text"),  pointsize = 2, alpha = 0.65,
     habillage=CADI.str$pop, repel = TRUE, axes = c(3, 4),
     palette = brewer.pal(length(unique(pop(CADI.str))), name = "Dark2"), invisible="quali")+
  scale_shape_manual(values=c(19,19,19,19,19,19))





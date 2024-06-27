# quick DAPC 

library(adegenet)

# load gind 
load("CATIgind.RData")

# data(CATIgind) 
# class(CATIgind)
# names(CATIgind)

grp <- find.clusters(CATIgind, max.n.clust=7, n.pca = 50,)

table(pop(CATIgind), grp$grp)

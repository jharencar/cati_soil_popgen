## Cati Neighbor Joining Tree
# Julia Harencar
# 2020 June 15

#set directory
setwd("/Users/juliaharencar/Google Drive (jharenca@ucsc.edu)/Cati_work")

#package for importing and converting .stru to geneind
library(adegenet)
# import data and convert to geneind (24 genotypes, 16384 loci, ea. genotype on two rows, individual names in first column, pop info in second column, header row with loci names, do not askfor otional info about the dataset)
CATIgind <- read.structure('CATI.structure.stru', n.ind = 24, n.loc = 128405, onerowperind = FALSE, col.lab = 1, col.pop = 2, row.marknames = 1, NA.char = 0)


# saving the genind object for future andlysis (bc it takes forever to load)
# save(CATIgind, file="CATIgind.RData")
# reloading genind object
load("CATIgind.RData")
library("poppr")
library("ape")

# Calculating dengrogram with bootstrap support (no strata, bootstrapping not within "population")
Cati_phylo <- aboot(CATIgind, strata = NULL, tree = "nj", distance = "nei.dist", root = FALSE)
plot.phylo(Cati_phylo, type = 'unrooted')
plot.phylo(Cati_phylo)

quartz.save("Cati_NJT.png", type="png", dpi = 1000)

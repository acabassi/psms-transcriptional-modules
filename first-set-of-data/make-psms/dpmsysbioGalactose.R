rm(list=ls())

library(klic)

galactose_output <- read.csv("GalactoseData_Results_Chain1.csv")
colnames(galactose_output)

galactose_cluster_samples <- galactose_output[2001:4000, 2:206]

psm <- matrix(0,
              dim(galactose_cluster_samples)[2],
              dim(galactose_cluster_samples)[2])
rownames(psm) <- colnames(psm) <- colnames(galactose_cluster_samples)
for(i in 1:dim(galactose_cluster_samples)[1]){
  for(j in 1:max(galactose_cluster_samples[i,])){
    psm <- psm + crossprod(galactose_cluster_samples[i,]==j)
  }
}

psm <- psm/(dim(galactose_cluster_samples)[1])
coph_corr <- copheneticCorrelation(psm)

heatmap(psm, scale = "none")

save(psm, coph_corr, file = "../dpmsysbio/psm_Galactose.RData")

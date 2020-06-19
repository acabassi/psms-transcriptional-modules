rm(list=ls())

library(klic)

harbison_output <- read.csv("HarbisonData_Results_Chain1.csv")
colnames(harbison_output)

harbison_cluster_samples <- harbison_output[2001:4000, 2:206]

psm <- matrix(0,
              dim(harbison_cluster_samples)[2],
              dim(harbison_cluster_samples)[2])
rownames(psm) <- colnames(psm) <- colnames(harbison_cluster_samples)
for(i in 1:dim(harbison_cluster_samples)[1]){
  for(j in 1:max(harbison_cluster_samples[i,])){
    psm <- psm + crossprod(harbison_cluster_samples[i,]==j)
  }
}

psm <- psm/(dim(harbison_cluster_samples)[1])
coph_corr <- copheneticCorrelation(psm)

heatmap(psm, scale = "none")

save(psm, coph_corr, file = "../dpmsysbio/psm_Harbison.RData")

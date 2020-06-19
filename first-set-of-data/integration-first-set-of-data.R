####################### Transcriptional module discovery #######################
########################### Unsupervised integration ###########################

rm(list = ls())

set.seed(1)

library(ComplexHeatmap)
library(colorspace)
library(rcartocolor)

# library(devtools)
# install_github("acabassi/coca")
library(coca)
# install_github("acabassi/klic")
library(klic)

############################ Load PSM of ChIP data #############################

# load("premium/harbison_chain1_psm_exclude_y_var_sel.RData")
load("dpmsysbio/psm_Harbison.RData")
ccH <- psm; rm(psm)
cophCorr_ccH <- coph_corr; rm(coph_corr)

########################## Load PSM of Ideker data #############################

# load("premium/galactose_chain1_psm_exclude_y_var_sel.RData")
load("dpmsysbio/psm_Galactose.RData")
ccM <- psm; rm(psm)
cophCorr_ccM <- coph_corr; rm(coph_corr)

############## Check that statistical units are in the same order ##############

sum(1-(rownames(ccM) == rownames(ccH))) # Should be zero

################################ Spectrum shift ################################

# Cophenetic correlation before spectrum shift
cophCorr_ccH_before <- copheneticCorrelation(ccH)
cophCorr_ccM_before <- copheneticCorrelation(ccM)

ccH <- spectrumShift(ccH,    verbose = TRUE, coeff = 1.01)
ccM <- spectrumShift(ccM, verbose = TRUE, coeff = 1.01)

cophCorr_ccH_after <- copheneticCorrelation(ccH)
cophCorr_ccM_after <- copheneticCorrelation(ccM)

######################### Create array of all matrices #########################

n_datasets <- 2
datasetNames <- c("Harbison", "Ideker")
N <- dim(ccH)[1]
CM <- array(NA, c(N, N, n_datasets))
dimnames(CM) <- list(rownames(ccH), rownames(ccH), datasetNames)

CM[rownames(ccH), colnames(ccH), "Harbison"] <- ccH
CM[rownames(ccM), colnames(ccM), "Ideker"]   <- ccM

################################# Integration ##################################
# Use localised multiple kernel k-means to integrate the datasets

# Initialise array of kernel matrices
maxK <- 13
# KM <- array(0, c(N, N, maxK - 1))
# clLabels <- array(NA, c(maxK - 1, N))
# theta <- array(NA, c(maxK - 1, N, n_datasets))
# 
# parameters <- list()
# parameters$iteration_count <-
#     100 # set the maximum number of iterations
# 
# for(i in 2:maxK) {
#     # Use kernel k-means with K=i to find weights and cluster labels
#     parameters$cluster_count <- i # set the number of clusters K
#     lmkkm <- lmkkmeans(CM, parameters, verbose = TRUE)
# 
#     # Compute weighted matrix
#     for (j in 1:dim(CM)[3]) {
#         KM[, , i - 1] <-
#             KM[, , i - 1] +
#             (lmkkm$Theta[, j] %*% t(lmkkm$Theta[, j])) * CM[, , j]
#     }
# 
#     # Save cluster labels
#     clLabels[i - 1, ] <- lmkkm$clustering
#     theta[i - 1, , ] <- lmkkm$Theta
# }

############################### Maximise silhouette ############################

# maxSil <- maximiseSilhouette(KM, clLabels, maxK = maxK)

# save(
#     KM,
#     clLabels,
#     theta,
#     maxSil,
#     file = "data/unsupervised_output_dpmsysbio.RData"
# )

load("data/unsupervised_output_dpmsysbio.RData")

png(paste0("figures/first_set_of_data_silhouette.png"),
    height = 550,
    width = 550)
plot(2:maxK, maxSil$silhouette,
     type = "b", xlab = "Clusters", ylab = "Average silhouette",cex.lab=1.5)
dev.off()

dimnames(theta) <- list(paste(as.character(2:maxK), "_clusters"),
                        rownames(CM),
                        c("ChIP", "Expression")
)

################################## Plot output #################################

bestK <- 5#maxSil$K
inds <- rownames(CM)
clustersBestK <- clLabels[bestK - 1, ]
weightsBestK <- theta[bestK-1,,]
names(clustersBestK) <- rownames(CM)

# Write clusters to csv files:
# write.table(
#     as.data.frame(clustersBestK),
#     paste(
#         "goto-scores/unsupervised_",
#         bestK,
#         "clusters_dpmsysbio.csv",
#         sep = ""
#     ),
#     col.names = FALSE,
#     quote = TRUE
# )
# 
# 
# shuffledClustersBestK <- sample(clustersBestK)
# table(clustersBestK)
# table(shuffledClustersBestK)
# names(shuffledClustersBestK) <- names(clustersBestK)
# sum(1-(clustersBestK==shuffledClustersBestK))
# write.table(
#     as.data.frame(shuffledClustersBestK),
#     paste(
#         "goto-scores/unsupervised_",
#         bestK,
#         "clusters_shuffled_dpmsysbio.csv",
#         sep = ""
#     ),
#     col.names = FALSE,
#     quote = TRUE
# )

# Let's just plot the clusters:
table(clustersBestK)
sortedClusters <- sort(clustersBestK, index.return = T)
myBlues <-
    colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)

my_anno_legend_param <- function(name, nrow=1) {
    return(list(
        title = name,
        labels_gp = gpar(fontsize = 18),
        title_gp = gpar(fontsize = 22),
        nrow = nrow
    ))
}

my_heatmap_legend_param <- list(labels_gp = gpar(fontsize = 18),
                                title_gp = gpar(fontsize = 22),
                                direction = "horizontal",
                                nrow = 1,
                                legend_width = unit(4.5, "cm"))
col_bw = c("white", "black")

palette25 = c("#B7BF10", # Light green
              "#4E5B31", # Core green
              "#115E67", # Core Cambridge blue
              "#85b09A", # Light Cambridge blue
              "#0072ce", # Core blue
              "#6CACE4", # Light blue
              "#E89CAE", # Light pink
              "#af95a6", # Light purple
              "#8C77A3", # Modified core purple
              "#D50032", # Core red
              "#E87722", # Core orange
              "#F1BE48", # Light yellow
              "#edf373", # Light green
              "#adbe86", # Core green
              "#57d4e3", # Core Cambridge blue
              "#c2d7cc", # Light Cambridge blue
              "#66bbff", # Core blue
              "#b5d5f1", # Light blue
              "#f3cdd6", # Light pink
              "#d7cad2", # Light purple
              "#c5bbd1", # Modified core purple
              "#ff6a8d", # Core red
              "#f3bb90", # Core orange
              "#f7dea3", # Light yellow
              "#000000") # Black

shuffled_palette25 <- palette25[sample(1:25,bestK)]
names(shuffled_palette25) <- as.character(1:bestK)


Hclusters <-
    rowAnnotation(
        finalClusters = as.factor(sortedClusters$x),
        name = "Final clusters",
        annotation_legend_param = my_anno_legend_param("Clusters", nrow=1),
        show_annotation_name = FALSE,
        col = list(finalClusters =  shuffled_palette25))
HwMatrix <-
    Heatmap(
        KM[sortedClusters$ix, sortedClusters$ix, bestK - 1],
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_column_names = FALSE,
        name = "Weighted kernel",
        col = myBlues,
        # row_split = sortedClusters$x,
        # column_split = sortedClusters$x,
        right_annotation = Hclusters,
        heatmap_legend_param = my_heatmap_legend_param,
        width = unit(20, "cm"),
        height = unit(20, "cm"))

png(
    paste0("figures/first_set_of_data_weightedKernel_", bestK, "clusters.png"),
    height = 660,
    width = 600
)
draw(
    HwMatrix,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

save(clustersBestK,
     shuffled_palette25,
     weightsBestK,
     file = "results/combined_clusters.RData")

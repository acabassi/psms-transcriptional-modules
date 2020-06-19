####################### Transcriptional module discovery #######################
############### Clustering the dataset of Harbison et al. (2004) ###############

rm(list = ls())

library(circlize)
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

################################ Spectrum shift ################################

cophCorr_ccH_before <- copheneticCorrelation(ccH)
ccH <- spectrumShift(ccH, verbose = TRUE, coeff = 1.01)
cophCorr_ccH_after <- copheneticCorrelation(ccH)

################################# Integration ##################################
# Use localised multiple kernel k-means to integrate the datasets

N <- dim(ccH)[1]

# Initialise array of kernel matrices
maxK <- 30
# clLabels <- array(NA, c(maxK - 1, N))
# ccH_rep <- array(NA, c(N, N, maxK - 1))
# 
# parameters <- list()
# 
# for(i in 2:maxK) {
#     # Use kernel k-means with K=i to find weights and cluster labels
#     parameters$cluster_count <- i # set the number of clusters K
#     set.seed(1)
#     kkm <- kkmeans(ccH, parameters)
#     # Save cluster labels
#     clLabels[i - 1, ] <- kkm$clustering
#     ccH_rep[,,i-1] <- ccH
# }

############################### Maximise silhouette ############################
# maxSil <- maximiseSilhouette(ccH_rep, clLabels, maxK = maxK)
# 
# save(clLabels,
#     maxSil,
#     file = "data/harbison_output_dpmsysbio.RData"
# )

load("data/harbison_output_dpmsysbio.RData")

################################## Plot output #################################

bestK <- 5# maxSil$K
inds <- rownames(ccH)
clustersBestK <- clLabels[bestK - 1, ]
names(clustersBestK) <- rownames(ccH)

# Write clusters to csv files:
# write.table(
#     as.data.frame(clustersBestK),
#     paste("goto-scores/harbison_", bestK, "clusters_dpmsysbio.csv", sep = ""),
#     col.names = FALSE,
#     quote = TRUE
# )

# shuffledClustersBestK <- sample(clustersBestK)
# table(clustersBestK)
# table(shuffledClustersBestK)
# names(shuffledClustersBestK) <- names(clustersBestK)
# sum(1-(clustersBestK==shuffledClustersBestK))
# write.table(
#     as.data.frame(shuffledClustersBestK),
#     paste("goto-scores/harbison_", bestK, "clusters_shuffled_dpmsysbio.csv", sep = ""),
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
                                at = c(0,1),
                                legend_width = unit(2.8, "cm"))
col_bw = c("white", "black")

harbisonDataNumeric <- as.matrix(read.csv("data/HarbisonData.csv",
                                          row.names = 1))
harbisonDataNumeric <- harbisonDataNumeric[rownames(ccH),]
harbCl <- as.factor(as.character(clustersBestK))

bigPalette <- c("#6CACE4", # Light blue
                "#E89CAE", # Light pink
                "#F1BE48", # Light yellow
                "#B7BF10", # Light green
                "#85b09A", # Light Cambridge blue
                "#af95a6", # Light purple
                "#D50032", # Core red
                "#0072ce", # Core blue
                "#E87722", # Core orange
                "#4E5B31", # Core green
                "#3F2A56", # Core purple
                "#115E67", # Core Cambridge blue
                "#003C71", # Dark blue
                "#8A1538") # Dark red
bigPalette <- bigPalette[7:12]
names(bigPalette) <- as.character(1:6)

HharClusters <- rowAnnotation(
    harbisonClusters = as.factor(sortedClusters$x),
    name = "ChIP Clusters",
    annotation_legend_param = my_anno_legend_param("ChIP clusters", nrow = 1),
    show_annotation_name = FALSE,
    col = list(
        harbisonClusters = bigPalette
    )
)

HharbisonData <-
    Heatmap(
        as.matrix(harbisonDataNumeric[sortedClusters$ix,]),
        cluster_rows = F,
        cluster_columns = T,
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        name = "ChIP data",
        right_annotation = HharClusters,
        col = col_bw,
        heatmap_legend_param = my_heatmap_legend_param
    )

png(
    paste0("figures/first_set_of_data_harbison_",
           bestK,
           "_clusters_dpmsysbio.png"),
    height = 440,
    width = 400
)
draw(
    HharbisonData,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

############################# Plot PSM and weights #############################

load("results/combined_clusters.RData")

HharbisonWeights <- 
    Heatmap(
        as.matrix(weightsBestK[sortedClusters$ix,"ChIP"]),
        cluster_rows = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        name = "Weight",
        heatmap_legend_param = my_heatmap_legend_param,
        col =  colorRamp2(c(0, 1), c("white", "#8A1538")),
        width = unit(0.5, "cm")
    )

HharbisonPSM <- 
    Heatmap(
        as.matrix(ccH[sortedClusters$ix,sortedClusters$ix]),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        name = "Similarity",
        right_annotation = HharClusters,
        col = myBlues,
        heatmap_legend_param = my_heatmap_legend_param
    )

png(
    paste0("figures/first_set_of_data_harbison_",
           bestK,
           "_clusters_psm_dpmsysbio.png"),
    height = 440,
    width = 400
)
draw(
    HharbisonWeights + HharbisonPSM,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

######################### Plot with combined clusters ##########################

load("results/combined_clusters.RData")
# This contains
# - weightsBestK
# - clustersBestK
# - shuffled_palette25 (only has 5 elements)

table(clustersBestK)
sortedClusters <- sort(clustersBestK, index.return = T)

HfinalClusters <-
    rowAnnotation(
        HarbisonClusters = as.factor(sortedClusters$x),
        name = "ChIP clusters",
        annotation_legend_param = my_anno_legend_param("ChIP clusters"),
        show_annotation_name = FALSE,
        col = list(HarbisonClusters = shuffled_palette25)
    )

HharbisonData <-
    Heatmap(
        as.matrix(harbisonDataNumeric[sortedClusters$ix,]),
        cluster_rows = F,
        cluster_columns = T,
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        name = "ChIP data",
        right_annotation = HfinalClusters,
        col = col_bw,
        heatmap_legend_param = my_heatmap_legend_param
    )


png("figures/first_set_of_data_harbison_final_clusters.png",
    height = 440,
    width = 400
)
draw(
    HharbisonData,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

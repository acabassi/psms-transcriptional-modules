####################### Transcriptional module discovery #######################
############### Clustering the dataset of Harbison et al. (2004) ###############

rm(list = ls())

library(ComplexHeatmap)
library(colorspace)
library(rcartocolor)

# library(devtools)
# install_github("acabassi/coca")
library(coca)
# install_github("acabassi/klic")
library(klic)

############################ Load PSM of ChIP data #############################

load("GTA_example_codes/Harbison_psm.RDa")
ccH <- Harbison_psm
rm(Harbison_psm)

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
#     file = "data/harbison_output.RData"
# )

load("data/harbison_output.RData")

################################## Plot output #################################

# bestK <- maxSil$K 
bestK <- 25
inds <- rownames(ccH)
clustersBestK <- clLabels[bestK - 1, ]
names(clustersBestK) <- rownames(ccH)

# Write clusters to csv files:
# write.table(
#     as.data.frame(clustersBestK),
#     paste("goto-scores/harbison_", bestK, "clusters.csv", sep = ""),
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
#     paste("goto-scores/harbison_", bestK, "clusters_shuffled.csv", sep = ""),
#     col.names = FALSE,
#     quote = TRUE
# )


############## Plot with clusters found based on the Harbison data #############
table(clustersBestK)
sortedClusters <- sort(clustersBestK, index.return = T)
myBlues <-
    colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)


my_anno_legend_param <- function(name, nrow=1) {
    return(list(
        title = name,
        labels_gp = gpar(fontsize = 18),
        title_gp = gpar(fontsize = 22),
        nrow = nrow,
        direction = "horizontal",
        grid_height = unit(5, "mm"),
        grid_width = unit(5, "mm")
    ))
}

my_heatmap_legend_param <- list(labels_gp = gpar(fontsize = 18),
                                title_gp = gpar(fontsize = 22),
                                direction = "horizontal",
                                nrow = 1,
                                at = c(0, 1),
                                legend_width = unit(5, "cm"))
col_bw = c("white", "black")

load("data/harbisonData.RData")
harbCl <- as.factor(as.character(clustersBestK))
names(harbCl) <- names(clustersBestK)
inds <- match(rownames(harbisonDataNumeric), names(harbCl))
harbCl <- harbCl[inds]

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

shuffled_palette25 <- palette25[sample(1:25,25)]
names(shuffled_palette25) <- as.character(1:25)

HharClusters <- rowAnnotation(
    harbisonClusters = as.factor(sortedClusters$x),
    name = "ChIP Clusters",
    annotation_legend_param = my_anno_legend_param("ChIP clusters", nrow = 3),
    show_annotation_name = FALSE,
    col = list(harbisonClusters = shuffled_palette25))

HharbisonData <-
    Heatmap(
        as.matrix(harbisonDataNumeric[sortedClusters$ix,]),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        # row_split = sortedClusters$x,
        # row_gap = unit(4, "mm"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        name = "ChIP data",
        right_annotation = HharClusters,
        col = col_bw,
        heatmap_legend_param = my_heatmap_legend_param,
        # width = unit(8, "cm")
    )
HharbisonData

png(
    paste0("figures/second_set_of_data_harbison_", bestK, "clusters.png"),
    height = 660,
    width = 600
)
draw(
    HharbisonData,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()
  

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
  paste0("figures/second_set_of_data_harbison_",
         bestK,
         "_clusters_psm_dpmsysbio.png"),
  height = 660,
  width = 600
)
draw(
  HharbisonPSM,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)
dev.off()

############################### Plot with weights ##############################

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
  paste0("figures/second_set_of_data_harbison_",
         bestK,
         "_clusters_psm_dpmsysbio.png"),
  height = 500,
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

table(clustersBestK)
sortedClusters <- sort(clustersBestK, index.return = T)
myBlues <-
    colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)


my_anno_legend_param <- function(name, nrow=1) {
    return(list(
        title = name,
        labels_gp = gpar(fontsize = 18),
        title_gp = gpar(fontsize = 22),
        nrow = nrow,
        direction = "horizontal"
    ))
}

my_heatmap_legend_param <- list(labels_gp = gpar(fontsize = 18),
                                title_gp = gpar(fontsize = 22),
                                direction = "horizontal",
                                nrow = 1,
                                legend_width = unit(5, "cm"))
col_bw = c("white", "black")

combCl <- as.factor(as.character(clustersBestK))
names(combCl) <- names(clustersBestK)
inds <- match(rownames(harbisonDataNumeric), names(combCl))
combCl <- combCl[inds]

HfinalClusters <- rowAnnotation(
    harbisonClusters = as.factor(sortedClusters$x),
    name = "Final Clusters",
    annotation_legend_param = my_anno_legend_param("Final clusters", nrow = 3),
    show_annotation_name = FALSE,
    col = list(harbisonClusters = shuffled_palette25))

HharbisonData <-
    Heatmap(
        as.matrix(harbisonDataNumeric[sortedClusters$ix,]),
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        # row_split = sortedClusters$x,
        # row_gap = unit(4, "mm"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        name = "ChIP data",
        right_annotation = HfinalClusters,
        col = col_bw,
        heatmap_legend_param = my_heatmap_legend_param,
        # width = unit(8, "cm")
    )
HharbisonData


png(
    paste0("figures/second_set_of_data_harbison_final_clusters.png"),
    height = 660,
    width = 600
)
draw(
    HharbisonData,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

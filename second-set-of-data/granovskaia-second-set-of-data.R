####################### Transcriptional module discovery #######################
############### Clustering the data of Granovskaia et al. (2010) ###############

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

load("GTA_example_codes/Harbison_psm.RDa")
ccH <- Harbison_psm
rm(Harbison_psm)

############################## Load Marina data ################################
load("GTA_example_codes/Marina_psm.RDa")

sum(1 - (rownames(ccH) == rownames(Marina_psm)))
sum(1 - (rownames(Marina_psm) == rownames(ccH)))

ccM <- Marina_psm
rm(Marina_psm)

################################ Spectrum shift ################################

cophCorr_ccM_before <- copheneticCorrelation(ccM)
ccM <- spectrumShift(ccM, verbose = TRUE, coeff = 1.01)
cophCorr_ccM_after <- copheneticCorrelation(ccM)

################################# Integration ##################################
# Use localised multiple kernel k-means to integrate the datasets
N <- dim(ccM)[1]

# Initialise array of kernel matrices
maxK <- 30
# clLabels <- array(NA, c(maxK - 1, N))
# ccM_rep <- array(0, c(N, N, maxK - 1))
# 
# parameters <- list()
# 
# for(i in 2:maxK) {
#     # Use kernel k-means with K=i to find weights and cluster labels
#     parameters$cluster_count <- i # set the number of clusters K
#     
#     set.seed(1)
#     kkm <- kkmeans(ccM, parameters)
#     
#     # Save cluster labels
#     clLabels[i - 1, ] <- kkm$clustering
#     
#     ccM_rep[,,i-1] <- ccM
# }

############################### Maximise silhouette ############################

# maxSil <- maximiseSilhouette(ccM_rep, clLabels, maxK = maxK)
# 
# save(clLabels,
#     maxSil,
#     file = "data/granovskaia_output.RData"
# )

load("data/granovskaia_output.RData")

################################## Plot output #################################

bestK <- 25# maxSil$K
clustersBestK <- clLabels[bestK - 1, ]
names(clustersBestK) <- rownames(ccM)

# Write clusters to csv files:
# write.table(
#     as.data.frame(clustersBestK),
#     paste0("goto-scores/granovskaia_", bestK, "clusters.csv"),
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
#     paste0("goto-scores/granovskaia_", bestK, "clusters_shuffled.csv"),
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
                                legend_width = unit(5, "cm"))
col_bw = c("white", "black")

load("data/marinaData.RData")
mariCl <- as.factor(as.character(clustersBestK))
names(mariCl) <- names(clustersBestK)
inds <- match(rownames(marinaData), names(mariCl))
mariCl <- mariCl[inds]

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

HmarClusters <- rowAnnotation(
    marinaClusters = as.factor(sortedClusters$x),
    name = "Expr. clusters",
    annotation_legend_param = my_anno_legend_param("Expr. clusters", nrow = 3),
    show_annotation_name = FALSE,
    col = list(marinaClusters = shuffled_palette25))

HmarinaData <-
    Heatmap(
        as.matrix(marinaData[sortedClusters$ix,]),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        # row_split = sortedClusters$x,
        # row_gap = unit(4, "mm"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        name = "Expr. data",
        right_annotation = HmarClusters,
        col = colorRamp2(c(min(marinaData), 0, max(marinaData)),
                         c("#8A1538", "white", "#003C71")),
        heatmap_legend_param = my_heatmap_legend_param,
        # width = unit(8, "cm")
    )

png(
    paste0("figures/second_set_of_data_marina_", bestK, "clusters.png"),
    height = 660,
    width = 600
)
draw(
    HmarinaData,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

HmarinaPSM <- 
    Heatmap(
        as.matrix(ccM[sortedClusters$ix,sortedClusters$ix]),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        name = "Similarity",
        right_annotation = HmarClusters,
        col = myBlues,
        heatmap_legend_param = my_heatmap_legend_param
    )

png(
    paste0("figures/second_set_of_data_granovskaia_",
           bestK,
           "_clusters_psm_dpmsysbio.png"),
    height = 660,
    width = 600
)
draw(
    HmarinaPSM,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

############################### Plot with weights ##############################

load("results/combined_clusters.RData")

my_heatmap_legend_param <- list(labels_gp = gpar(fontsize = 18),
                                title_gp = gpar(fontsize = 22),
                                direction = "horizontal",
                                nrow = 1,
                                at = c(0, 1),
                                legend_width = unit(5, "cm"))

HmarinaWeights <- 
    Heatmap(
        as.matrix(weightsBestK[sortedClusters$ix, "Expression"]),
        cluster_rows = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        name = "Weight",
        heatmap_legend_param = my_heatmap_legend_param,
        col =  colorRamp2(c(0, 1), c("white", "#8A1538")),
        width = unit(0.5, "cm")
    )

HmarinaPSM <- 
    Heatmap(
        as.matrix(ccM[sortedClusters$ix,sortedClusters$ix]),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        name = "Similarity",
        right_annotation = HmarClusters,
        col = myBlues,
        heatmap_legend_param = my_heatmap_legend_param
    )

png(
    paste0("figures/second_set_of_data_granovskaia_",
           bestK,
           "_clusters_psm_dpmsysbio.png"),
    height = 500,
    width = 400
)
draw(
    HmarinaWeights + HmarinaPSM,
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

HfinalClusters <- rowAnnotation(
    harbisonClusters = as.factor(sortedClusters$x),
    name = "Final Clusters",
    annotation_legend_param = my_anno_legend_param("Final clusters", nrow = 3),
    show_annotation_name = FALSE,
    col = list(harbisonClusters = shuffled_palette25))

HmarinaData <-
    Heatmap(
        as.matrix(marinaData[sortedClusters$ix,]),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        # row_split = sortedClusters$x,
        # row_gap = unit(4, "mm"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        name = "Expression data",
        right_annotation = HfinalClusters,
        col = colorRamp2(c(min(marinaData), 0, max(marinaData)),
                         c("#8A1538", "white", "#003C71")),
        heatmap_legend_param = my_heatmap_legend_param,
        # width = unit(8, "cm")
    )
HmarinaData

png(
    "figures/second_set_of_data_marina_final_clusters.png",
    height = 660,
    width = 600
)
draw(
    HmarinaData,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

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
#     set.seed(1)
#     kkm <- kkmeans(ccM, parameters)
#     # Save cluster labels
#     clLabels[i - 1, ] <- kkm$clustering
#     ccM_rep[,,i-1] <- ccM
# }

############################### Maximise silhouette ############################
 
# maxSil <- maximiseSilhouette(ccM_rep, clLabels, maxK = maxK)
# 
# save(clLabels,
#     maxSil,
#     file = "data/ideker_output_dpmsysbio.RData"
# )

load("data/ideker_output_dpmsysbio.RData")

################################## Plot output #################################

bestK <- 5#maxSil$K
clustersBestK <- clLabels[bestK - 1, ]
names(clustersBestK) <- rownames(ccM)

# Write clusters to csv files:
# write.table(
#     as.data.frame(clustersBestK),
#     paste0("goto-scores/ideker_", bestK, "clusters_dpmsysbio.csv"),
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
#     paste0("goto-scores/ideker_", bestK, "clusters_shuffled_dpmsysbio.csv"),
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
                                legend_width = unit(5, "cm"))
col_bw = c("white", "black")

bigPalette <- c("#B7BF10", # Light green
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

ideClusters_palette <- bigPalette[c(9,1,7,4,8,6)]
names(ideClusters_palette) <- as.character(1:6)

idekerData <- as.matrix(read.csv("data/GalactoseData.csv", row.names = 1))
HideClusters <-
    rowAnnotation(
        IdekerClusters = as.factor(sortedClusters$x),
        name = "Expr. clusters",
        annotation_legend_param = my_anno_legend_param("Expr. clusters"),
        show_annotation_name = FALSE,
        col = list(IdekerClusters = ideClusters_palette)
    )

palette <- colorRampPalette(c("#8A1538", "white", "#003C71"))(8)

HidekerData <- Heatmap(
    as.matrix(idekerData[sortedClusters$ix,]),
    cluster_rows = F,
    cluster_columns = T,
    # row_split = sortedClusters$x,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_column_dend = FALSE,
    name = "Expr. data",
    right_annotation = HideClusters,
    heatmap_legend_param = my_heatmap_legend_param,
    col = c("1"=palette[2], "2"="white", "3"=palette[7])
)

png(
    paste0("figures/first_set_of_data_ideker_",
           bestK,
           "_clusters_dpmsysbio.png"),
    height = 440,
    width = 400
)
draw(
    HidekerData,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

############################# Plot PSM and weights #############################

load("results/combined_clusters.RData")
# This contains
# - weightsBestK
# - clustersBestK
# - shuffled_palette25 (only has 5 elements)

my_heatmap_legend_param <- list(labels_gp = gpar(fontsize = 18),
                                title_gp = gpar(fontsize = 22),
                                direction = "horizontal",
                                nrow = 1,
                                at = c(0,1),
                                legend_width = unit(2.8, "cm"))
HidekerWeights <- 
    Heatmap(
        as.matrix(weightsBestK[sortedClusters$ix,"Expression"]),
        cluster_rows = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        name = "Weight",
        heatmap_legend_param = my_heatmap_legend_param,
        col =  colorRamp2(c(0, 1), c("white", "#8A1538")),
        width = unit(0.5, "cm")
    )

HidekerPSM <-
    Heatmap(
        as.matrix(ccM[sortedClusters$ix,sortedClusters$ix]),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        name = "Similarity",
        right_annotation = HideClusters,
        col = myBlues,
        heatmap_legend_param = my_heatmap_legend_param
    )

png(
    paste0("figures/first_set_of_data_ideker_",
           bestK,
           "_clusters_psm_dpmsysbio.png"),
    height = 440,
    width = 400
)
draw(
    HidekerWeights + HidekerPSM,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

######################### Plot with combined clusters ##########################

# Start using combined clusters (uploaded above) from now on
table(clustersBestK)
sortedClusters <- sort(clustersBestK, index.return = T)

my_heatmap_legend_param <- list(labels_gp = gpar(fontsize = 18),
                                title_gp = gpar(fontsize = 22),
                                direction = "horizontal",
                                nrow = 1,
                                legend_width = unit(5, "cm"))

HfinalClusters <-
    rowAnnotation(
        IdekerClusters = as.factor(sortedClusters$x),
        name = "Expr. clusters",
        annotation_legend_param = my_anno_legend_param("Expr. clusters"),
        show_annotation_name = FALSE,
        col = list(IdekerClusters = shuffled_palette25)
    )

HideData <-
    Heatmap(
        as.matrix(idekerData[sortedClusters$ix,]),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        name = "Expression data",
        right_annotation = HfinalClusters,
        heatmap_legend_param = my_heatmap_legend_param,
        col = c("1"=palette[2], "2"="white", "3"=palette[7])
    )


png("figures/first_set_of_data_ideker_final_clusters.png",
    height = 440,
    width = 400
)
draw(
    HideData,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()

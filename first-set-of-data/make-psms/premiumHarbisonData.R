rm(list = ls())

library(gplots)
library(klic)
library(PReMiuM)

## HOW TO USE PREMIUM (UNSUPERVISED CASE)

# 1: Read in the data

myData                  <- read.csv("../data/HarbisonData.csv", row.names = 1)
covNames                <- names(myData)
# Even in the unsupervised case, we require a (dummy) column in the data matrix,
# for the response
outcome                 <- c(seq(0,0,length = ceiling(nrow(myData)/2)),
                             seq(1,1,length = floor(nrow(myData)/2)))
myData                  <- cbind(outcome, myData) # This is the data file 
                                            # (includes response and covariates)

# Visualise the data:
#heatmap(as.matrix(myData), Rowv= NA, Colv = NA)

# 2: Put data and modelling options into an "inputs" data frame
inputs <- # This is just to initialise the "inputs" variable, so that it has all
          # of the right fieldnames
  generateSampleDataFile(clusSummaryBernoulliDiscrete()) 
inputs$inputData <- myData   #Put our data in the "inputData" slot
inputs$covNames <- covNames
inputs$nCovariates <- length(covNames)
inputs$fixedEffectNames <- NULL
chain <- 1

# 3: Run PReMiuM!
# For the unsupervised case, we set "excludeY = TRUE"
prof_regr <- profRegr(
  yModel = inputs$yModel,
  xModel = inputs$xModel,
  nSweeps = 100000,
  nClusInit = 20,
  nBurn = 1,
  data = inputs$inputData,
  output = paste0("../premium/harbison_chain", chain),
  covNames = inputs$covNames,
  reportBurnIn = TRUE,
  excludeY = TRUE,
  seed = chain
)

save(prof_regr,
     file = paste0("../premium/harbison_",
                   chain,"output_prof_regr_exclude_y_var_sel.RData"))

dissimObj <- PReMiuM::calcDissimilarityMatrix(prof_regr)
dissMat <- PReMiuM::vec2mat(dissimObj$disSimMat, nrow = length(outcome))
psm <- 1-dissMat
coph_corr <- copheneticCorrelation(psm)
rownames(psm) <- colnames(psm) <- rownames(myData)

save(psm, coph_corr,
     file = paste0(
       "../premium/harbison_chain", chain,"_psm_exclude_y_var_sel.RData"))

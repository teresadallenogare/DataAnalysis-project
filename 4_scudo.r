# ============================ SCUDO CLASSIFICATION ==============================
# Author: Teresa Dalle Nogare
# Version : 06-04-2023
# - scudo classification
# ==============================================================================
# ====================== Install packages  ===================================== 
install.packages("BiocManager")
install.packages("igraph") 
install.packages("ggplot2")
install.packages('caret')
install.packages("glmnet") 
install.packages('Matrix')
install.packages('igraph')
BiocManager::install("GEOquery")
BiocManager::install("useful")
BiocManager::install("genefilter")
BiocManager::install("rScudo")
BiocManager::install('clusterProfiler')
BiocManager::install("AnnotationDbi") 
BiocManager::install("org.Hs.eg.db") 
BiocManager::install("HsAgilentDesign026652.db")

library("GEOquery")
library("useful")
library(ggplot2)
library("RColorBrewer")
library("genefilter")
library("MASS")
library("gplots")
library('ggplot2')
library('lattice')
library("caret")
library("glmnet")
library("rScudo")
library('igraph')
library('clusterProfiler')
library("AnnotationDbi") 
library("org.Hs.eg.db")
library('HsAgilentDesign026652.db') # agilent conversion id probes
# ============================= SCUDO ==================================
load(file = "~/Desktop/DA-project/0_data/Y.Rdata")
set.seed(123)
Y <- Y[, c(rep(1:3), rep(10:12), rep(19:21), rep(28:30), rep(4:9), rep(13:18), rep(22:27), rep(31:35))]
# 0 : FBS , 1 : GCRa-COMB
group01 <- c(rep(0,12), rep(1,23))
y <- c(group01)
f <- factor(y, labels = c('FBS', 'GCRA-COMB'))
#Y1 <- Y[,c(rep(1:18))]
#Y2 <- Y[,c(rep(19:35))]
n.train <- 15 # Y2_small_clear : 18 (2/3)
n.test <- 10# Y3 : 17 (1/3)

# take index at random
idx.train <- sample(c(rep(1:6), rep(13:24)), n.train) # random sample 12 indexes from 1 to 18 (train dataset in first 18 raws of dat_tot)
idx.test <- sample(c(rep(7:12), rep(25:35)),n.test)   # random sample 6 indexes from 19 to 35 (test dataset in last 17 rows of dat_tot)
idx.train
idx.test

Y_train <- Y[,idx.train]
Y_test <- Y[,idx.test]

trainRes <- scudoTrain(Y_train, groups = f[idx.train], nTop = 50, nBottom = 50, alpha = 0.05)
show(trainRes)
upS <- upSignatures(trainRes)
cons <- consensusUpSignatures(trainRes)
distMat <- distMatrix(trainRes)
distMat
# inspect signatures 
#upSignatures(trainRes)[1:10,1:10] 
#consensusUpSignatures(trainRes)[1:5, ]
# generate and plot map of training samples 
trainNet <- scudoNetwork(trainRes, N = 0.45) 
scudoPlot(trainNet, vertex.label = NA)

# perform validation using testing samples
testRes <- scudoTest(trainRes, Y_test, f[idx.test],
                     nTop = 50, nBottom = 50, foldChange = TRUE) 
testRes
testNet <- scudoNetwork(testRes, N = 0.4)
scudoPlot(testNet, vertex.label = NA)

# identify clusters on map
testClust <- igraph::cluster_spinglass(testNet, spins = 2) 
plot(testClust, testNet, vertex.label = NA)

# perform classification
classRes <- scudoClassify(Y_train, Y_test, N = 0.4,
                          nTop = 50, nBottom = 50,
                          trainGroups = f[idx.train], alpha = 0.05) 
classRes
caret::confusionMatrix(classRes$predicted, f[idx.test])

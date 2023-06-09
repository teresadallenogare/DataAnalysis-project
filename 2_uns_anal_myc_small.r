# ============================ UNSUPERVISED METHODS ============================
# Author: Teresa Dalle Nogare
# Version : 27-04-2023
# - hierarchical clustering
# - k-means clusering
# ==============================================================================
# ====================== Install packages  ===================================== 
install.packages("BiocManager")
#install.packages('ggfortify')
BiocManager::install("GEOquery")
BiocManager::install("useful")
library("GEOquery")
library("useful")
library(ggplot2)
library(ggrepel)
#library('ggfortify')


# ------------------------------- GPL21185 -------------------------------------
load(file = "./0_data/Y3.Rdata")

# 1. Hierarchical clustering
dist_matrix3 <- dist(t(Y3)) # pass the transpose of the data matrix 
hc_result3 <- hclust(dist_matrix3, method = 'ave')
k <- 6
grp3 <- cutree(hc_result3, k=k)
table(grp3)
par(mar=c(5.1, 4.1, 4.1, 2.1))
plot(hc_result3, hang <- -1, labels=grp3)
rect.hclust(hc_result3, k = 6, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups

# 2. K-means
k <- 6 # expect 6 clusters since 3 treatments for each cell type
kmeans_result3 <- kmeans(t(Y3), k)
table(kmeans_result3$cluster)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(kmeans_result3, data=t(Y3)) + geom_text_repel(aes(label = colnames(Y3)), size = 3)

# ---------------------------- GPL21185 - clear --------------------------------
load(file = "./0_data/Y3_clear.Rdata")

# 1. Hierarchical clustering
dist_matrix3_c <- dist(t(Y3_clear)) # pass the transpose of the data matrix 
hc_result3_c <- hclust(dist_matrix3_c, method = 'ave')
k <- 6
grp3_c <- cutree(hc_result3_c, k=k)
table(grp3_c)
par(mar=c(5.1, 4.1, 4.1, 2.1))
plot(hc_result3_c, hang <- -1, labels=grp3_c)
rect.hclust(hc_result3_c, k = 6, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups
# same dendrogram as the one with no cleaned data

# 2. K-means
k <- 6 # expect 6 clusters since 3 treatments for each cell type
kmeans_result3_c <- kmeans(t(Y3_clear), k)
table(kmeans_result3_c$cluster)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(kmeans_result3_c, data=t(Y3_clear)) + geom_text_repel(aes(label = colnames(Y3_clear)), size = 3)


# ---------------------------- GPL17077 - myc ----------------------------------
load(file = "./0_data/Y2_myc.Rdata")

# 1. Hierarchical clustering
dist_matrix2_myc <- dist(t(Y2_myc)) # pass the transpose of the data matrix 
hc_result2_myc <- hclust(dist_matrix2_myc, method = 'ave')
k <- 6
grp2_myc <- cutree(hc_result2_myc, k=k)
table(grp2_myc)
par(mar=c(5.1, 4.1, 4.1, 2.1))
plot(hc_result2_myc, hang <- -1, labels=grp2_myc)
rect.hclust(hc_result2_myc, k = 6, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups

# 2. K-means
k <- 6 # expect 6 clusters since 3 treatments for each cell type
kmeans_result2_myc <- kmeans(t(Y2_myc), k)
table(kmeans_result2_myc$cluster)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(kmeans_result2_myc, data=t(Y2_myc)) + geom_text_repel(aes(label = colnames(Y2_myc)), size = 3)

# --------------------------- GPL17077 - small ---------------------------------
load(file = "./0_data/Y2_small.Rdata")

# 1. Hierarchical clustering
dist_matrix2_small <- dist(t(Y2_small)) # pass the transpose of the data matrix 
hc_result2_small <- hclust(dist_matrix2_small, method = 'ave')
k <- 6
grp2_small <- cutree(hc_result2_small, k=k)
table(grp2_small)
par(mar=c(5.1, 4.1, 4.1, 2.1))
plot(hc_result2_small, hang <- -1, labels=grp2_small)
rect.hclust(hc_result2_small, k = 6, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups

# 2. K-means
k <- 6 # expect 6 clusters since 3 treatments for each cell type
kmeans_result2_small <- kmeans(t(Y2_small), k)
table(kmeans_result2_small$cluster)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(kmeans_result2_small, data=t(Y2_small)) + geom_text_repel(aes(label = colnames(Y2_small)), size = 3)

# ------------------------ GPL17077 - small clear ------------------------------
load(file = "./0_data/Y2_small_clear.Rdata")

# 1. Hierarchical clustering
dist_matrix2_c <- dist(t(Y2_small_clear)) # pass the transpose of the data matrix 
hc_result2_c <- hclust(dist_matrix2_c, method = 'ave')
k <- 6
grp2_c <- cutree(hc_result2_c, k=k)
table(grp2_c)
par(mar=c(5.1, 4.1, 4.1, 2.1))
plot(hc_result2_c, hang <- -1, labels=grp2_c)
rect.hclust(hc_result2_c, k = 6, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups
# same dendrogram as the one with no cleaned data

# 2. K-means
k <- 6 # expect 6 clusters since 3 treatments for each cell type
kmeans_result2_c <- kmeans(t(Y2_small_clear), k)
table(kmeans_result2_c$cluster)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(kmeans_result2_c, data=t(Y2_small_clear)) + geom_text_repel(aes(label = colnames(Y2_small_clear)), size = 3)


# ------------------------ GPL17077+GPL21185 cleaned ------------------------------
# ==============================================================================
# Note :
# I expect 4 big clusters, each of which should contain 3 sub-clusters
# The 4 clusters are the one due to SCLC and NSCLC for Y2_small_clear and Y3_clear
# The 3 subclusters are the three different treatments
# Here I suppose to find the 4 main clusters, that arise when I consider both Y2_small_clear
# and Y3_clear together
# ==============================================================================
load(file = "./0_data/Y.Rdata")
# 1. Hierarchical clustering
dist_matrix <- dist(t(Y)) # pass the transpose of the data matrix 
hc_result <- hclust(dist_matrix, method = 'ave')
k <- 4
grp <- cutree(hc_result, k=k)
table(grp)
par(mar=c(5.1, 4.1, 4.1, 2.1))
plot(hc_result, hang <- -1, labels=grp)
rect.hclust(hc_result, k = 4, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups

# 2. K-means
k <- 4
kmeans_result <- kmeans(t(Y), k)
table(kmeans_result$cluster)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(kmeans_result, data=t(Y)) + geom_text_repel(aes(label = colnames(Y)), size = 3)


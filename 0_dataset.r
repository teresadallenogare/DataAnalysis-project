# ======================= IMPORT PACKAGES AND DATASET DEFINITION ===============
# Author: Teresa Dalle Nogare
# Version : 27-04-2023
# - import dataset and cleans it
# ==============================================================================
# ====================== Install packages  ===================================== 
install.packages("BiocManager")
#install.packages('ggfortify')
BiocManager::install("GEOquery")
BiocManager::install("useful")
library("GEOquery")
library("useful")
library(ggplot2)
#library('ggfortify')

# ========================= Import GSE66245 ====================================
gse <- getGEO("GSE66245",destdir=".") 
# -------------------------------- GPL13534 ------------------------------------
# Methylation profile 
gse1 <- gse[[1]]
ex1 <- exprs(gse1)
Y1 <- scale(ex1)

boxplot(Y1, 
        xlab = 'Channel', ylab = ' Mean Signal', main = 'GPL13534',
        aces = F)
save(Y1, file = './0_data/Y1.Rdata')
# -------------------------------- GPL17077 ------------------------------------
# Expression profile 1 
# clear by column 
gse2 <- gse[[2]]
ex2 <- exprs(gse2)
ex2_lung <- ex2[,c(1:55)] # only lung source name : discard first 9 samples
ex2_myc <- ex2_lung[,c(rep(1:36), rep(55))] # only MYC cell line
ex2_small <- ex2_myc[,c(rep(10:18), rep(28:36))] # only shBRG1 samples

log_ex2 <- log2(ex2)
log_ex2_lung <- log2(ex2_lung)
log_ex2_myc <- log2(ex2_myc)
log_ex2_small <- log2(ex2_small)

Y2 <- scale(log_ex2)
Y2_lung <- scale(log_ex2_lung)
Y2_myc <- scale(log_ex2_myc)
Y2_small <- scale(log_ex2_small)
boxplot(Y2_small, 
        xlab = 'Channel', ylab = 'Log Mean Signal', main = 'GPL17077 - small',
        aces = F)

save(Y2, file = './0_data/Y2.Rdata')
save(Y2_lung, file ='./0_data/Y2_lung.Rdata')
save(Y2_myc, file = './0_data/Y2_myc.Rdata')
save(Y2_small, file = './0_data/Y2_small.Rdata')
# -------------------------------- GPL21185 ------------------------------------
# Expression profile 2
gse3 <- gse[[3]]
ex3 <- exprs(gse3)
log_ex3 <- log2(ex3)
Y3 <- scale(log_ex3)
boxplot(Y3, 
        xlab = 'Channel', ylab = 'Log Mean Signal', main = 'GPL21185',
        aces = F)
save(Y3, file = './0_data/Y3.Rdata')
# ------------------------ GPL17077-small + GPL21185 ---------------------------
# clear by row
ex <- merge(ex2_small,ex3, by = 'row.names', all = TRUE)
ex <- na.omit(ex) 
rownames(ex) <- ex[,1] # add rownames as header
ex <- ex[,-1]
log_ex <- log2(ex)
Y <- scale(log_ex)

boxplot(Y, 
        xlab = 'Channel', ylab = 'Log Mean Signal', main = 'GPL17077 + GPL21185',
        aces = F)
save(Y, file = './0_data/Y.Rdata')

# separate in two datasetss
ex2_small_clear <- ex[,c(rep(1:18))] 
ex3_clear <- ex[,c(rep(19:35))]

log_ex2_small_clear <- log2(ex2_small_clear)
log_ex3_clear <- log2(ex3_clear)
Y2_small_clear <- scale(log_ex2_small_clear)
Y3_clear <- scale(log_ex3_clear)
boxplot(Y2_small_clear, 
        xlab = 'Channel', ylab = 'Log Mean Signal', main = 'GPL17077 ',
        aces = F)

boxplot(Y3_clear, 
        xlab = 'Channel', ylab = 'Log Mean Signal', main = 'GPL21185 ',
        aces = F)

save(Y2_small_clear, file = "./0_data/Y2_small_clear.Rdata")
save(Y3_clear, file = "./0_data/Y3_clear.Rdata")



# ============================== PCA ===========================================
# Author: Teresa Dalle Nogare
# Version : 27-04-2023
# - PCA of whole dataset and smaller ones
# ====================== Install packages  ===================================== 
install.packages("BiocManager")
#install.packages('ggfortify')
BiocManager::install("GEOquery")
BiocManager::install("useful")
library("GEOquery")
library("useful")
library(ggplot2)
#library('ggfortify')

# -------------------------------- PCA GPL21185 --------------------------------
load(file = "./0_data/Y3.Rdata")
pca3 <- prcomp(t(Y3)) # t(): transposes
head(pca3)
summary(pca3)
var_explained = pca3$sdev^2 / sum(pca3$sdev^2)
#create scree plot
qplot(c(1:11),var_explained[1:11])+
  geom_line() + 
  xlab("Principal Component index") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot : GPL21185") +
  ylim(0, 0.7)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
pch = rep(c(21,22,24), each = 3)
grpcol <- c(rep('chocolate1',9), rep('mediumpurple1',8))
plot(pca3$x[,1], pca3$x[,3], xlab='PC1', ylab='PC3', main= 'PCA for components 1&3', type='p', pch=pch, col='black', bg = grpcol, lwd = 0.9)
legend("topright", inset=c(-0.3,0), legend=c( 'FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col = c('black'), pt.bg = c('chocolate1'),title = 'SCLC')
legend("topright", inset=c(-0.3,0.3), legend=c('FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col= c('black'), pt.bg = c('mediumpurple1'), title = 'NSCLC')
text(pca3$x[,1][1:3],pca3$x[,3][1:3], labels = rownames(pca3$x)[1:3], cex=0.5, adj = c(-0.4,0.3))
text(pca3$x[,1][4],pca3$x[,3][4], labels = rownames(pca3$x)[4], cex=0.5, adj = c(-0.4,0))
text(pca3$x[,1][5],pca3$x[,3][5], labels = rownames(pca3$x)[5], cex=0.5, adj = c(-0.4,-0.5))
text(pca3$x[,1][6],pca3$x[,3][6], labels = rownames(pca3$x)[6], cex=0.5, adj = c(-0.4,-1))
text(pca3$x[,1][7],pca3$x[,3][7], labels = rownames(pca3$x)[7], cex=0.5, adj = c(-0.4,0))
text(pca3$x[,1][8],pca3$x[,3][8], labels = rownames(pca3$x)[8], cex=0.5, adj = c(-0.4,0.8))
text(pca3$x[,1][9],pca3$x[,3][9], labels = rownames(pca3$x)[9], cex=0.5, adj = c(-0.4,1))
text(pca3$x[,1][10],pca3$x[,3][10], labels = rownames(pca3$x)[10], cex=0.5, adj = c(1.4,0.5))
text(pca3$x[,1][11:12],pca3$x[,3][11:12], labels = rownames(pca3$x)[11:12], cex=0.5, adj = c(1.4,0))
text(pca3$x[,1][13:14],pca3$x[,3][13:14], labels = rownames(pca3$x)[13:14], cex=0.5, adj = c(1.4,1))
text(pca3$x[,1][15],pca3$x[,3][15], labels = rownames(pca3$x)[15], cex=0.5, adj = c(1.4,0))
text(pca3$x[,1][16:17],pca3$x[,3][16:17], labels = rownames(pca3$x)[16:17], cex=0.5, adj = c(1.4,0.3))


plot(pca3$x[,1], pca3$x[,2], xlab='PCA1', ylab='PCA2', main= 'PCA for components 1&2', type='p', pch=pch, col='black', bg = grpcol, lwd = 0.9)
legend("topright", inset=c(-0.3,0), legend=c( 'FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col = c('black'), pt.bg = c('chocolate1'),title = 'SCLC')
legend("topright", inset=c(-0.3,0.3), legend=c('FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col= c('black'), pt.bg = c('mediumpurple1'), title = 'NSCLC')
plot3d(pca3$x,col=grp_trt )

# -------------------------------- PCA GPL21185-clear --------------------------------
load(file = './0_data/Y3_clear.Rdata')
pca3_clear <- prcomp(t(Y3_clear)) # t(): transposes
head(pca3_clear)
summary(pca3_clear)
var_explained = pca3_clear$sdev^2 / sum(pca3_clear$sdev^2)
#create scree plot
qplot(c(1:11),var_explained[1:11])+
  geom_line() + 
  xlab("Principal Component index") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot : GPL21185 - clear") +
  ylim(0, 0.9)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
pch = rep(c(21,22,24), each = 3)
grpcol <- c(rep('chocolate1',9), rep('mediumpurple1',8))
plot(pca3_clear$x[,1], pca3_clear$x[,3], xlab='PC1', ylab='PC3', main= 'PCA for components 1&3', type='p', pch=pch, col='black', bg = grpcol, lwd = 0.9)
legend("topright", inset=c(-0.3,0), legend=c( 'FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col = c('black'), pt.bg = c('chocolate1'),title = 'SCLC')
legend("topright", inset=c(-0.3,0.3), legend=c('FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col= c('black'), pt.bg = c('mediumpurple1'), title = 'NSCLC')
text(pca3_clear$x[,1][1:2],pca3_clear$x[,3][1:2], labels = rownames(pca3_clear$x)[1:2], cex=0.5, adj = c(-0.4,0))
text(pca3_clear$x[,1][3],pca3_clear$x[,3][3], labels = rownames(pca3_clear$x)[3], cex=0.5, adj = c(-0.4,0.5))
text(pca3_clear$x[,1][4],pca3_clear$x[,3][4], labels = rownames(pca3_clear$x)[4], cex=0.5, adj = c(-0.4,-0.2))
text(pca3_clear$x[,1][5],pca3_clear$x[,3][5], labels = rownames(pca3_clear$x)[5], cex=0.5, adj = c(-0.4,-0.3))
text(pca3_clear$x[,1][6],pca3_clear$x[,3][6], labels = rownames(pca3_clear$x)[6], cex=0.5, adj = c(-0.4,-0.5))
text(pca3_clear$x[,1][7],pca3_clear$x[,3][7], labels = rownames(pca3_clear$x)[7], cex=0.5, adj = c(-0.4,0.5))
text(pca3_clear$x[,1][8],pca3_clear$x[,3][8], labels = rownames(pca3_clear$x)[8], cex=0.5, adj = c(-0.4,0))
text(pca3_clear$x[,1][9],pca3_clear$x[,3][9], labels = rownames(pca3_clear$x)[9], cex=0.5, adj = c(-0.4,0.5))
text(pca3_clear$x[,1][10],pca3_clear$x[,3][10], labels = rownames(pca3_clear$x)[10], cex=0.5, adj = c(1.4,0))
text(pca3_clear$x[,1][11],pca3_clear$x[,3][11], labels = rownames(pca3_clear$x)[11], cex=0.5, adj = c(1.4,0.5))
text(pca3_clear$x[,1][12],pca3_clear$x[,3][12], labels = rownames(pca3_clear$x)[12], cex=0.5, adj = c(1.4,0))
text(pca3_clear$x[,1][13:14],pca3_clear$x[,3][13:14], labels = rownames(pca3_clear$x)[13:14], cex=0.5, adj = c(1.4,0))
text(pca3_clear$x[,1][15],pca3_clear$x[,3][15], labels = rownames(pca3_clear$x)[15], cex=0.5, adj = c(1.4,1.4))
text(pca3_clear$x[,1][16:17],pca3_clear$x[,3][16:17], labels = rownames(pca3_clear$x)[16:17], cex=0.5, adj = c(1.4,0))

# ------------------------- PCA GPL17077 - myc -------------------------------
load(file = './0_data/Y2_myc.Rdata')
pca2_myc <- prcomp(t(Y2_myc)) # t(): transposes
head(pca2_myc)
summary(pca2_myc)
var_explained = pca2_myc$sdev^2 / sum(pca2_myc$sdev^2)
#create scree plot
qplot(c(1:11),var_explained[1:11])+
  geom_line() + 
  xlab("Principal Component index") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot : GPL17077-MYC") +
  ylim(0, 0.7)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
grpcol <- c(rep('cornflowerblue',18), rep('indianred1',19))
pch = rep(c(21,22,24), each = 3)
# 12
plot(pca2_myc$x[,1], pca2_myc$x[,2], xlab='PC1', ylab='PC2', main= 'PCA for components 1&2', type='p', pch=pch, col='black', bg = grpcol, lwd = 0.9)
legend("topright", inset=c(-0.3,0), legend=c( 'FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col=c('black'), pt.bg = c('cornflowerblue'), title = 'SCLC')
legend("topright", inset=c(-0.3,0.3), legend=c('FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col=c('black'), pt.bg = c('indianred1'), title = 'NSCLC')
# 13
plot(pca2_myc$x[,1], pca2_myc$x[,3], xlab='PC1', ylab='PC3', main= 'PCA for components 1&3', type='p', pch=pch, col='black', bg = grpcol, lwd = 0.9)
legend("topright", inset=c(-0.3,0), legend=c( 'FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col=c('black'), pt.bg = c('cornflowerblue'), title = 'SCLC')
legend("topright", inset=c(-0.3,0.3), legend=c('FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col=c('black'), pt.bg = c('indianred1'), title = 'NSCLC')

# ------------------------- PCA GPL17077 - small -------------------------------
load(file = './0_data/Y2_small.Rdata')
pca2_small <- prcomp(t(Y2_small)) # t(): transposes
head(pca2_small)
summary(pca2_small)
var_explained = pca2_small$sdev^2 / sum(pca2_small$sdev^2)
#create scree plot
qplot(c(1:11),var_explained[1:11])+
  geom_line() + 
  xlab("Principal Component index") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot : GPL17077") +
  ylim(0, 0.7)
grpcol <- c(rep('cornflowerblue',9), rep('indianred1',9))
plot(pca2_small$x[,1], pca2_small$x[,2], xlab='PC1', ylab='PC2', main= 'PCA for components 1&2', type='p', pch=pch, col=c('black'), bg = grpcol)
legend("topright", inset=c(-0.3,0), legend=c( 'FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col=c('black'), pt.bg = c('cornflowerblue'), title = 'SCLC')
legend("topright", inset=c(-0.3,0.3), legend=c('FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col=c('black'), pt.bg = c('indianred1'), title = 'NSCLC')
text(pca2_small$x[,1][1:9],pca2_small$x[,2][1:9], labels = rownames(pca2_small$x)[1:3], cex=0.5, adj = c(-0.4,0))
text(pca2_small$x[,1][10:14],pca2_small$x[,2][10:14], labels = rownames(pca2_small$x)[10:14], cex=0.5, adj = c(1.4,0))
text(pca2_small$x[,1][15],pca2_small$x[,2][15], labels = rownames(pca2_small$x)[15], cex=0.5, adj = c(1.4,1))
text(pca2_small$x[,1][16:18],pca2_small$x[,2][16:18], labels = rownames(pca2_small$x)[16:18], cex=0.5, adj = c(1.4,0))

plot(pca2_small$x[,1], pca2_small$x[,3], xlab='PC1', ylab='PC3', main= 'PCA for components 1&3', type='p', pch=pch, col=c('black'), bg = grpcol)
legend("topright", inset=c(-0.3,0), legend=c( 'FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col=c('black'), pt.bg = c('cornflowerblue'), title = 'SCLC')
legend("topright", inset=c(-0.3,0.3), legend=c('FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col=c('black'), pt.bg = c('indianred1'), title = 'NSCLC')


# --------------------- PCA GPL17077 - small clear -----------------------------
load(file = './0_data/Y2_small_clear.Rdata')
pca2_small_clear <- prcomp(t(Y2_small_clear)) # t(): transposes
head(pca2_small_clear)
summary(pca2_small_clear)
var_explained_ex2 = pca2_small_clear$sdev^2 / sum(pca2_small_clear$sdev^2)
#create scree plot
qplot(c(1:11),var_explained_ex2[1:11])+
  geom_line() + 
  xlab("Principal Component index") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot : GPL17077 - small clear") +
  ylim(0, 0.7)
grpcol <- c(rep('cornflowerblue',9), rep('indianred1',9))
plot(pca2_small_clear$x[,1], pca2_small_clear$x[,2], xlab='PC1', ylab='PC2', main= 'PCA for components 1&2', type='p', pch=pch, col=c('black'), bg = grpcol)
legend("topright", inset=c(-0.3,0), legend=c( 'FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col=c('black'), pt.bg = c('cornflowerblue'), title = 'SCLC')
legend("topright", inset=c(-0.3,0.3), legend=c('FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col=c('black'), pt.bg = c('indianred1'), title = 'NSCLC')
text(pca2_small_clear$x[,1][1:9],pca2_small_clear$x[,2][1:9], labels = rownames(pca2_small_clear$x)[1:9], cex=0.5, adj = c(-0.4,0))
text(pca2_small_clear$x[,1][10:14],pca2_small_clear$x[,2][10:14], labels = rownames(pca2_small_clear$x)[10:14], cex=0.5, adj = c(1.4,0))
text(pca2_small_clear$x[,1][15],pca2_small_clear$x[,2][15], labels = rownames(pca2_small_clear$x)[15], cex=0.5, adj = c(1.4,1))
text(pca2_small_clear$x[,1][16:18],pca2_small_clear$x[,2][16:18], labels = rownames(pca2_small_clear$x)[16:18], cex=0.5, adj = c(1.4,0))

plot(pca2_small_clear$x[,1], pca2_small_clear$x[,3], xlab='PC1', ylab='PC3', main= 'PCA for components 1&3', type='p', pch=pch, col=c('black'), bg = grpcol)
legend("topright", inset=c(-0.3,0), legend=c( 'FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col=c('black'), pt.bg = c('cornflowerblue'), title = 'SCLC')
legend("topright", inset=c(-0.3,0.3), legend=c('FBS', 'GC-RA', 'COMB' ), pch=rep(c(21,22,24)), col=c('black'), pt.bg = c('indianred1'), title = 'NSCLC')


# ----------------- cleaned data together --------------------------------------
load(file = './0_data/Y.Rdata')
pca_clear <- prcomp(t(Y)) # t(): transposes
head(pca_clear)
summary(pca_clear)
var_explained = pca_clear$sdev^2 / sum(pca_clear$sdev^2)
#create scree plot
qplot(c(1:11),var_explained[1:11])+
  geom_line() + 
  xlab("Principal Component index") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot ") +
  ylim(0, 0.7)
# [   Y2_small,  Y3*   ]
# [SCLC, NSCLC, SCLC, NSCLC]
grpcol <- c(rep('darkgoldenrod1',9), rep('darkolivegreen4',9), rep('cornflowerblue',9), rep('indianred1',8))
pch = rep(c(21,22,24), each = 3)
plot(pca_clear$x[,1], pca_clear$x[,2], xlab='PC1', ylab='PC2', main= 'PCA for components 1&2', type='p', pch=pch, col=c('black'), bg = grpcol)
legend("topright", inset=c(-0.3,0), legend=c( 'FBS', 'GC-RA', 'COMB', 'FBS*', 'GC-RA*', 'COMB*' ), pch=rep(c(21,22,24, 21, 22, 24)), col=c('black'), pt.bg = c('darkgoldenrod1', 'darkgoldenrod1', 'darkgoldenrod1', 'cornflowerblue', 'cornflowerblue', 'cornflowerblue'), title = 'SCLC')
legend("topright", inset=c(-0.3,0.5), legend=c('FBS', 'GC-RA', 'COMB' , 'FBS*', 'GC-RA*', 'COMB*' ), pch=rep(c(21,22,24, 21, 22, 24)), col=c('black'), pt.bg = c('darkolivegreen4', 'darkolivegreen4', 'darkolivegreen4', 'indianred1', 'indianred1', 'indianred1'), title = 'NSCLC')
# * indicates Y3.

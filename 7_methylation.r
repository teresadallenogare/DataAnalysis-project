# ======================= METHYLATION PROIFLE ===================================
# Author: Teresa Dalle Nogare
# Version : 23-05-2023

# ====================== Install packages  ===================================== 
install.packages("BiocManager")
#install.packages('ggfortify')
BiocManager::install("GEOquery")
BiocManager::install("useful")
BiocManager::install("mixOmics")
library("GEOquery")
library("genefilter")
library("MASS")
library(ggplot2)
library(mixOmics)
library("useful")
library("glmnet")
library("AnnotationDbi") 
library("org.Hs.eg.db")
library('HsAgilentDesign026652.db') # agilent conversion id probes

# methylation profile
load(file = "~/Desktop/DA-project/0_data/Y1.Rdata")
# expression profile of GPL21185
load(file = "~/Desktop/DA-project/0_data/Y3.Rdata")

# ------------------------------- prepare datasets ----------------------------------
# create expression and methylation profiles to compare - pay attention to have same order of samples
# expression profile
X <- Y3[, c(1,10, 4, 7, 13, 16)] 
groupX <- c(rep('FBS',2), rep('GCRA-COMB',4))
fX <- factor(groupX)
# row ttest
set.seed(1234)
ttX<- rowttests(X,fX)
keepers <- which(ttX$p.value<0.1) # no use p-adjust
length(keepers)
X_ttest <- X[keepers,]
tX_ttest <- t(X_ttest)
# methylation profile
Y <- Y1[, c(4, 10, 5, 6, 11, 12)]
groupY <- c(rep('FBS',2), rep('GCRA-COMB',4))
fY <- factor(groupY)
# row ttest
set.seed(1234)
ttY<- rowttests(Y,fY)
keepers <- which(ttY$p.value<0.1) # no use p-adjust
length(keepers)
Y_ttest <- Y[keepers,]
tY_ttest <- t(Y_ttest)

# ------------------------------- SPLS ----------------------------------
#group_total <- c(groupX, groupY)
#ftotal <- factor(group_total)

spls.lung <- spls(tX_ttest, tY_ttest, ncomp = 2, 
                  keepX = c(25,25),
                  keepY = c(10,10),
                  mode = "regression") 

# ------------------------------- Plots ----------------------------------
group_cancerType <- c(rep('SCLC', 1), rep('NSCLC', 1), rep('SCLC', 2), rep('NSCLC', 2) )
f_cancerType <- factor(group_cancerType)


plotIndiv(spls.lung, ind.names = TRUE, 
          rep.space = "X-variate", # plot in X-variate subspace
          group = groupX,
          pch = f_cancerType, 
          col.per.group = color.mixo(1:2), 
          ellipse = TRUE,
          legend = TRUE, legend.title = 'Treatment', legend.title.pch = 'Cancer type')

plotIndiv(spls.lung, ind.names = FALSE,
          rep.space = "Y-variate", # plot in Y-variate subspace
          group = groupY, # colour by time group
          pch = f_cancerType, 
          col.per.group = color.mixo(1:2), 
          ellipse = TRUE,
          legend = TRUE, legend.title = 'Treatment', legend.title.pch = 'Cancer type')

plotIndiv(spls.lung, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          group = groupX, # colour by time group
          pch = f_cancerType, # select symbol
          col.per.group = color.mixo(1:2),  
          ellipse = TRUE,
          legend = TRUE, legend.title = 'Treatment', legend.title.pch = 'Cancer type',
          title = 'Individuals plot')

plotArrow(spls.lung, ind.names = TRUE,
          method = 'shrinkage',
          group = groupX, # colour by time group
          col.per.group = color.mixo(1:2),
          legend.title = 'Treatment',
          cex = c(10,10))


plotVar(spls.lung, cex = c(3,3), legend = TRUE, style="ggplot2",var.names = c(FALSE, TRUE), abline = TRUE)
selectVar(spls.lung, comp = 1)$X$name 
selectVar(spls.lung, comp = 2)$X$name

selectVar(spls.lung, comp = 1)$Y$name
selectVar(spls.lung, comp = 2)$Y$name
plotLoadings(spls.lung, method = 'mean', contrib = 'max', subtitle = c('Expression', 'Methylation')) # depict weight assigned to each of these variables


loadX_df <- data.frame(spls.lung$loadings$X)
loadY_df <- data.frame(spls.lung$loadings$Y)
# remove all rows with both zeros
loadX_df <- loadX_df[rowSums(loadX_df[])>0,]
loadY_df <- loadY_df[rowSums(loadY_df[])>0,]

probe.namesX  <- rownames(loadX_df)
probe.namesX_df     <- data.frame(probe.namesX)

gene.entrez  <- lookUp(probe.namesX, "HsAgilentDesign026652", "ENTREZID")
gene.ensembl <- lookUp(probe.namesX, "HsAgilentDesign026652", "ENSEMBL")
gene.symbols <- lookUp(probe.namesX, "HsAgilentDesign026652", "SYMBOL")

probe.namesX_df["ENZYME"]  <- list(as.character(gene.entrez))
probe.namesX_df["SYMBOL"]  <- list(as.character(gene.symbols))
probe.namesX_df["ENSEMBL"] <- list(as.character(gene.ensembl))

# coefficiente di variabilità per fare su ogni riga (o varianza più alta) if variance greater than 5%
# siti di metilazione che sono posizioni su dna





# ============================ GSEA 2 ==============================
# Author: Teresa Dalle Nogare
# Version : 12-04-2023
# Differential gene expression analysis - to determine changes in expression 
# between different conditions - aka treated vs non-treated 
# ==============================================================================
# ====================== Install packages  ===================================== 
install.packages("BiocManager")
install.packages("igraph") 
install.packages("ggplot2")
install.packages('caret')
install.packages("glmnet") 
install.packages('Matrix')
install.packages('igraph')
install.packages('gprofiler2')
BiocManager::install("GEOquery")
BiocManager::install("useful")
BiocManager::install("genefilter")
BiocManager::install("rScudo")
BiocManager::install('DESeq2')
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
library('DESeq2')
library('gprofiler2')
library('limma')
library("XML")
library(annotate)

load(file = "~/Desktop/DA-project/0_data/Y.Rdata")
# order dataset based on treatment type and group
Y <- Y[, c(rep(1:3), rep(10:12), rep(19:21), rep(28:30), rep(4:9), rep(13:18), rep(22:27), rep(31:35))]
group <- c(rep('FBS',12), rep('GCRACOMB',23))
f <- factor(group)
design <- model.matrix(~ f+0 , data.frame(Y))
fit <- lmFit(Y, design)
cts <- paste('fFBS', 'fGCRACOMB', sep="-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", number = length(Y))
tT_df <- data.frame(tT)

# redefinition of group
group <- c(rep('FBS',12), rep('GCRA-COMB',23))
f <- factor(group)

# Random forest : selection of 100 most relevant genes
set.seed(1234) # deterministic algorithm
rf <- randomForest(x=t(Y), y=as.factor(group), ntree=1000)
probe.names  <- rownames(rf$importance)
top100index  <- order(rf$importance,decreasing = TRUE)[1:100]
top100       <- probe.names[order(rf$importance,decreasing = TRUE)[1:100]]
top100df     <- data.frame(top100)
# select only the 100 genes in top100
tT <- tT[row.names(tT) %in% top100, ]
# sort data in the same order of ttY
tT <- tT[match(top100, rownames(tT)), ]

gene.entrez  <- lookUp(top100, "HsAgilentDesign026652", "ENTREZID")
gene.ensembl <- lookUp(top100, "HsAgilentDesign026652", "ENSEMBL")
gene.symbols <- lookUp(top100, "HsAgilentDesign026652", "SYMBOL")

top100df["ENZYME"]  <- list(as.character(gene.entrez))
top100df["SYMBOL"]  <- list(as.character(gene.symbols))
top100df["ENSEMBL"] <- list(as.character(gene.ensembl))

# row ttest
set.seed(1234)
ttY <- rowttests(Y,f)

top100df["p.value"] <- ttY$p.value[top100index]
top100df["logFC"]   <- tT$logFC

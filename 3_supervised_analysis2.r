# ============================ SUPERVISED METHODS 2 ==============================
# Author: Teresa Dalle Nogare
# Version : 12-04-2023
# Note : work with cleared data only. 
# Classification based on the type of treatment. the type of treatment is contained in the
# field 'treatment:ch1'
#
# ==============================================================================
# ====================== Install packages  ===================================== 
install.packages("BiocManager")
install.packages("randomForest")
#install.packages('ggfortify')
install.packages("ggplot2")
install.packages("e1071")
install.packages('caret')
install.packages("glmnet") 
install.packages('Matrix')
BiocManager::install("GEOquery")
BiocManager::install("useful")
BiocManager::install(" randomForest")
BiocManager::install("pROC")
BiocManager::install("genefilter")
library("GEOquery")
library("useful")
library(ggplot2)
#library('ggfortify')
library("randomForest") 
library("RColorBrewer")
library("pROC")
library("genefilter")
library("MASS")
library("gplots")
library('ggplot2')
library("e1071")
library('lattice')
library("caret")
library("glmnet")
# ============================= Random forest ==================================
load(file = "~/Desktop/DA-project/0_data/Y.Rdata")
# order dataset based on treatment type and group
Y <- Y[, c(rep(1:3), rep(10:12), rep(19:21), rep(28:30), rep(4:9), rep(13:18), rep(22:27), rep(31:35))]
group <- c(rep('FBS',12), rep('GCRA-COMB',23))

set.seed(1234) # deterministic algorithm
rf <- randomForest(x=t(Y), y=as.factor(group), ntree=1000)
plot(rf, main = NULL)
rf
rf$votes
# importance 1 
plot(sort(rf$importance, decreasing=TRUE), xlab = 'Gene index', ylab = 'Importance', xlim= range(1:400))
sort(rf$importance, decreasing=TRUE)
probe.names <- rownames(rf$importance)
# importance 2
varImpPlot(rf, n.var=25, main='ALL Subset Results')
# Look at variable selection
#install.packages('varSelRF')
#library(varSelRF)
#set.seed(1234)

#rfsel <- varSelRF(t(Y),as.factor(group),
#                  ntree=1000, ntreeIterat=200, vars.drop.frac=0.02)
#rf.sig.gn <- rfsel$selected.vars 
#View(rf.sig.gn)

# top N genes
top30 <- probe.names[order(rf$importance, decreasing = TRUE)[1:30]]
top100 <- probe.names[order(rf$importance, decreasing=TRUE)[1:100]]
top200 <- probe.names[order(rf$importance, decreasing=TRUE)[1:200]]
#write.csv(top30, file = "~/Desktop/DA-project/gsea_data/FBS-GCRA/probes-top30_Y_whole.txt", quote=FALSE, row.names = FALSE, col.names=FALSE)
#write.csv(top100, file = "~/Desktop/DA-project/gsea_data/FBS-GCRA/probes-top100_Y_whole.txt", quote=FALSE, row.names = FALSE, col.names=FALSE)
#write.csv(top200, file = "~/Desktop/DA-project/gsea_data/FBS-GCRA/probes-top200_Y_whole.txt", quote=FALSE, row.names = FALSE, col.names=FALSE)

# heatmap
imp.temp <- abs(rf$importance[,])
t <- order(imp.temp,decreasing=TRUE) 
# Get subset of expression values for 25 most 'important' 
gn.imp <- names(imp.temp)[t]
gn.25 <- gn.imp[1:25] # vector of top 25 genes, in order 
t <- is.element(rownames(Y,),gn.25)
sig.Y <- Y[t,] # matrix of expression values,
hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256) 
colnames(sig.Y) <- group # This will label the heatmap 
csc <- rep(hmcol[50],35 )
csc[group=='GCRA-COMB'] <- hmcol[200]
par(oma=c(0,0,0,2.75))
heatmap.2 (sig.Y, 
           scale="row", 
           col=hmcol,
           ColSideColors=csc, 
           key = TRUE,
           trace = 'none',
           density = 'none' ,
           main = NULL,
           labRow = NULL,
           labCol = NULL)


# ============== Random forest with cross validation ===========================
f <- factor(group)
ttY <- rowttests(Y,f)
keepers <- which(p.adjust(ttY$p.value)<1e-1)
length(keepers)
Y_tt <- Y[keepers,]
tY_tt <- t(Y_tt) 
dat_tot <- cbind(as.data.frame(tY_tt), f)
colnames(dat_tot)[ncol(dat_tot)] <- 'TREAT'

# Run algorithm using 10-fold cross validation
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"
# CANCER_TYPE~. : left is dependent variable, right is dependent variable
fit.lda <- train(TREAT~., data=dat_tot, method="lda", metric=metric, trControl=control)
fit.rf <- train(TREAT~., data=dat_tot, method="rf", metric=metric, trControl=control)
results <- resamples(list(LDA=fit.lda, RF=fit.rf)) 
summary(results)
ggplot(results) + labs(y = "Accuracy")

fit.rf
# ================================= LDA ========================================
# consider the whole dataset
f <- factor(group)
ttY <- rowttests(Y,f) 
head(ttY)
keepers <- which(p.adjust(ttY$p.value)<1e-1)
length(keepers)
Y_tt <- Y[keepers,]
tY_tt <- t(Y_tt) 
# Add column with labels 
dat_tot <- cbind(as.data.frame(tY_tt), f)
colnames(dat_tot)[ncol(dat_tot)] <- 'TREAT'
n.train <- 15 # Y2_small_clear : 18 (2/3)
n.test <- 10  # Y3 : 17 (1/3)

# take index at random
idx.train <- sample(c(rep(1:6), rep(13:24)), n.train) # random sample 12 indexes from 1 to 18 (train dataset in first 18 raws of dat_tot)
idx.test <- sample(c(rep(7:12), rep(25:35)),n.test)   # random sample 6 indexes from 19 to 35 (test dataset in last 17 rows of dat_tot)
idx.train
idx.test

mod <- lda(TREAT~., data = dat_tot, prior = c(0.5,0.5), subset = idx.train )
plot(mod)
# TRAIN
mod.values <- predict(mod, dat_tot[idx.train,]) 
mod.values$class
plot(mod.values$x[,1], ylab='LDA Axis', main= 'LDA - training', type='p', pch=18, col=c('black') )
text(mod.values$x[,1], col=c(as.numeric(dat_tot[idx.train,"TREAT"])), cex=0.8, adj = c(1.6,0))
# TEST
preds<-predict(mod, dat_tot[idx.test,])
preds$class
plot(preds$x[,1], ylab=c("LDA Axis"), main= 'LDA - testing', type='p', pch=18, col=c('black') )
text(preds$x[,1],col=c(as.numeric(dat_tot[idx.test,"TREAT"])+18), cex=0.8, adj = c(1.6,0))

table(as.numeric(preds$class),
      as.numeric(dat_tot[idx.test, "TREAT"])+18)
roc_lda <- plot.roc(as.numeric(preds$class), as.numeric(dat_tot[idx.test, "TREAT"])+18 )

# ============================== LDA with CARET ================================
group <- c(rep('FBS',12), rep('GCRA-COMB',23))
f <- factor(group)
ttY <- rowttests(Y,f)
keepers <- which(p.adjust(ttY$p.value)<1e-1)
length(keepers)
Y_tt <- Y[keepers,]
tY_tt <- t(Y_tt) 
dat_tot <- cbind(as.data.frame(tY_tt), f)
colnames(dat_tot)[ncol(dat_tot)] <- 'TREAT'

# Run algorithm using 10-fold cross validation
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"
fit.lda <- train(TREAT~., data=dat_tot, method="lda", metric=metric, trControl=control)
fit.lda
summary(fit.lda)
# TRAIN
pred_lda_train <-  predict(fit.lda, dat_tot[idx.train,])
rownames(dat_tot[idx.train,])
pred_lda_train
# TEST
pred_lda_test = predict(fit.lda, dat_tot[idx.test,])
rownames(dat_tot[idx.test,])
pred_lda_test
# ============================== LASSO ================================
# '0' : FBS, '1' : GCRA-COMB
group01 <- c(rep(0,12), rep(1,23))
y <- c(group01)
f <- factor(y)
Y_t <- t(Y)
fit=glmnet(Y_t,y,standardize=FALSE,family="binomial") 
plot(fit, xvar = 'lambda', label=TRUE)
cfit=cv.glmnet(Y_t,y,standardize=FALSE,family="binomial")
plot(cfit)

coef(cfit, s=cfit$lambda.min)

# ------------- Repeat analysis with training and testing --------------------
n.train <- 15 # Y2_small_clear : 18 (2/3)
n.test <- 10  # Y3 : 17 (1/3)

# take index at random
idx.train <- sample(c(rep(1:6), rep(13:24)), n.train) # random sample 12 indexes from 1 to 18 (train dataset in first 18 raws of dat_tot)
idx.test <- sample(c(rep(7:12), rep(25:35)),n.test)   # random sample 6 indexes from 19 to 35 (test dataset in last 17 rows of dat_tot)
idx.train
idx.test

# fit on training subset
fit=glmnet(Y_t[idx.train,],y[idx.train],standardize=FALSE,family="binomial") 
plot(fit,xvar = 'lambda', label=TRUE)
cfit=cv.glmnet(Y_t[idx.train,],y[idx.train],standardize=FALSE,family="binomial") 
plot(cfit)
cfit$lambda.min

# TRAIN
pred_train <- predict(fit,Y_t[idx.train,], type="class", s= cfit$lambda.min)
pred_train
plot(pred_train, ylab=c('class value'), main= 'LASSO - training', type='p', pch=18, col=c('black') )
text(pred_train,col=c(as.numeric(dat_tot[idx.train,"TREAT"])+18), cex=0.8, adj = c(1.6,0))
# TEST
pred_test <- predict(fit,Y_t[idx.test,], type="response", s=cfit$lambda.min)
pred_test
plot(pred_test, ylab=c('class value'), main= 'LASSO - testing', type='p', pch=18, col=c('black') )
text(pred_test,col=c(as.numeric(dat_tot[idx.test,"TREAT"])+18), cex=0.8, adj = c(1.6,0))



# ============================== LASSO  with CARET ================================
group <- c(rep('FBS',12), rep('GCRA-COMB',23))
f <- factor(group)
ttY <- rowttests(Y,f)
keepers <- which(p.adjust(ttY$p.value)<1e-1)
length(keepers)
Y_tt <- Y[keepers,]
tY_tt <- t(Y_tt) 
dat_tot <- cbind(as.data.frame(tY_tt), f)
colnames(dat_tot)[ncol(dat_tot)] <- 'TREAT'
control <- trainControl(method="cv", number=10) 
metric <- "Accuracy"
fit.lasso <- train(TREAT~., data=dat_tot, method="glmnet", family = "binomial",
                   tuneGrid = expand.grid(alpha = 1, lambda = seq(0,1,by=0.05)),
                   trControl = control, metric = metric)
plot(fit.lasso, xlab = 'lambda')
fit.lasso$finalModel
fit.lasso$bestTune

# comparison with other classification methods 
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"
fit.rf <- train(TREAT~., data=dat_tot, method="rf",metric=metric, trControl=control)
fit.lda <- train(TREAT~., data=dat_tot, method="lda",metric=metric, trControl=control)
results <- resamples(list(RF=fit.rf, LDA=fit.lda, Lasso=fit.lasso)) 
summary(results)
ggplot(results) + labs(y = "Accuracy")

# Run algorithms using 10-fold cross validation, 10 times
control <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
fit.rf.2 <- train(TREAT~., data=dat_tot, method="rf", metric=metric, trControl=control)
fit.lda.2 <- train(TREAT~., data=dat_tot, method="lda", metric=metric, trControl=control)
results <- resamples(list(LDA=fit.lda.2, RF=fit.rf.2))
ggplot(results) + labs(y = "Accuracy")

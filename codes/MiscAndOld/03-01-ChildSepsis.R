# Focusing on the subsets of children sepsis data
# Generally, I use a variance cut-off filter along with a lasso feature
# selection for identifying a small set of pathways.

rm(list = ls())

library(pheatmap)
library(limma)
library(tidyverse)

child.mat <- readRDS("DATA/childMat.rds")
labs      <- readRDS("DATA/childTab2.rds")


png("images/ChildrenSubtypesGPL570.png")
pheatmap(child.mat, annotation_col = labs,show_rownames = F,
         clustering_method = "ward.D2" )
dev.off()


high.paths <- (apply(child.mat,1, var))

mean.cutoff <- mean(high.paths)
plotDensities(high.paths, main = "Pathway Variances in Children GPL570")
var.cutoff <- sd(as.matrix(child.mat))

#child.mat  <- child.mat[high.paths > (mean.cutoff + (1 * var.cutoff)), ]
child.mat  <- child.mat[high.paths >1, ]
#de.var <- rownames(child.mat) %in% rownames(tT.select)
#child.mat  <- child.mat[de.var, ]
ncol(child.mat)
nrow(child.mat)
rownames(child.mat)


#labs$age <- as.factor(round(labs$age))
write.csv(high.paths[high.paths>1],"Reports/ChildVaryingPathways.csv")
# Heatmap below is amazing. Does no show batch effects in big subsets

png("images/ChildrenVariablePathwaysSubtypes.png")
init.map <- pheatmap(child.mat, annotation_col = labs,show_rownames = F,
         clustering_method = "ward.D2")
dev.off()

clusts <- cutree(init.map$tree_col, k = 2)
clusts <- paste("G",clusts,sep = "")



labs2   <- cutree(init.map$tree_row, k = 2)



library(caret) 
library(pROC)

set.seed(2)
validation_index <- createDataPartition(clusts, p=0.7, list=FALSE)
gset.test        <-  child.mat[,-validation_index]
gene.exp.tot     <-  child.mat[,validation_index]        
#Using all the data
#gene.exp.tot     <- child.mat

train.set        <- as.data.frame(t(gene.exp.tot)) 
train.set$sml <- as.factor(clusts)[validation_index]
#train.set$sml    <- as.factor(clusts)


tests.set     <- as.data.frame(t(gset.test)) 
tests.set.sml <- as.factor(clusts)[-validation_index]
### Lasso on the train data for further dimention reduction
control <- trainControl(method="repeatedcv",number = 3,repeats = 10,
                        summaryFunction = twoClassSummary, classProbs = TRUE)
#control <- trainControl(method="LOOCV",
#                        summaryFunction = twoClassSummary, classProbs = TRUE)
#control <- trainControl(method="cv",number = 5, returnResamp = "all")


set.seed(7)
fit.lasso <- train(sml ~ .,train.set,
                   method = "glmnet",
                   trControl = control,
                   metric = "ROC",
                   preProc = c("center", "scale"),
                   tuneGrid = expand.grid(.alpha = c(1),
                                          .lambda = c((1:10)/100)))


probs.glm <- predict(fit.lasso,tests.set,type = "prob")
ROC.glm   <- roc(tests.set.sml, probs.glm[,"G2"], levels = levels(tests.set.sml))
plot(ROC.glm)
fit.lasso$results
fit.lasso$bestTune

#results <- resamples(list("glmNet" = fit.lasso))

aa <- varImp(fit.lasso, scale = FALSE)
plotDensities(unlist(aa$importance))
#names(aa$importance[aa$importance > 0])
lasso.feats <- rownames(aa$importance)
lasso.feats <- lasso.feats[which(aa$importance > 0.00)] 
lasso.feats <- gsub("`","",lasso.feats)
print(length(lasso.feats))  


write.csv(aa$importance, "Reports/LassoChildrenPathways.csv")

lasso.feats <- read.csv("Reports/S3_LassoChildrenPathways.csv")
lasso.feats <- lasso.feats$Pathway.Name
lasso.feats <- as.character(lasso.feats)

net.mat <- child.mat[lasso.feats,]


clusts <- cutree(init.map$tree_col, k = 2)

names(clusts) == rownames(labs)
labs$InitClusts <- clusts
png("images/FinalChildren.png")
pheatmap(net.mat,annotation_col = labs, clustering_method = "ward.D2")
dev.off()


PathExpTab  <- readRDS("DATA/pathGeneTable.rds")
temp.tab    <-PathExpTab %>% filter(., Pathway %in% lasso.feats) 

write.csv(unique(temp.tab$SYMBOL),file = "Reports/CheckThisGenes.csv")


# paths <- split(temp.tab$SYMBOL,temp.tab$Pathway)
# inds  <- sapply(paths, length)
# paths <- paths[inds > 0]


# 
# 
# 
# 
# cov.mat <- stats::cov(t(net.mat))
# invcov <- abs(round(solve(cov.mat),3))
# A1     <- 1*(invcov > 4)
# diag(A1) <-0
# 
# hist(invcov)
# 
# seps.network <- network::as.network.matrix(A1)
# 
# ?GGally::ggnet2
# set.seed(1)
# GGally::ggnet2(seps.network,node.size = 5, node.label = rownames(A1),
#        arrow.size = 4,label.size = 2,label.trim = T,arrow.gap = 0.015,
#        mode = "fruchtermanreingold") +
#     theme(legend.title=element_blank())



library("gplots")
venn(paths[4:6])


library(ggfortify)
### What's next is to find translation of the subgroups

net.pca <- prcomp(t(net.mat),center = TRUE, scale. = TRUE)

pdf("images/SelectPCAnormalizedChilds570.pdf")
autoplot(net.pca, data = labs, colour = "InitClusts")
dev.off()



### Some other type of clustering...
### Ignore the following unless you like to explore model based clustering

library(mclust)


BIC <- mclustBIC(t(child.mat))
plot(BIC)
summary(BIC)

BIC.r <- mclustBIC((child.mat))
plot(BIC.r)
summary(BIC.r)



mod.mRNA <- Mclust(t(child.mat), G = 2, modelName = "VVE")
mod.mRNA.r <- Mclust((child.mat), G = 3, modelName = "EII")

mod.mRNA$classification

clust.order <- unlist(tapply(1:ncol((child.mat)),
                             as.factor(mod.mRNA$classification),I)
                      ,use.names = F)
mRNA.mat    <- data.matrix(child.mat)[,clust.order]


clust.order <- unlist(tapply(1:ncol(t(child.mat)),
                             as.factor(mod.mRNA.r$classification),I)
                      ,use.names = F)
mRNA.mat    <- data.matrix(mRNA.mat)[clust.order,]


pheatmap(mRNA.mat, cluster_cols = F,cluster_rows = F, show_rownames  = F,annotation_col = labs)

# Just some color setting for visualization

my_group    <- as.numeric(as.factor(mod.mRNA$classification))
my_col1     <- brewer.pal(8, "Set1")[my_group]
my_col1     <- my_col1[clust.order]


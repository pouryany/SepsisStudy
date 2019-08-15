# Preliminary analysis using Pathway Summary Statistics (Rank Based)
# Other methodologies are investigated in 02-XX Codes. 

rm(list = ls())

gc()
library(pheatmap)
library(limma)
library(tidyverse)

run1.exp.mat <- readRDS("DATA/CorrectedSepsisGPL570.rds")
samp.tab     <- readxl::read_xlsx("DATA/TrueSepsisDataUsedforFinalAnalysisMay2015.xlsx")
left.out     <- samp.tab %>% filter(.,study != "GSE26378" & 
                                      platform == "GPL570") %>% 
                             select(.,arrayID) %>% pull()
PathExpTab   <- readRDS("DATA/pathGeneTable.rds")

MedianUp     <- function(zz){mean(zz[zz > median(zz,na.rm = T)],na.rm =T)}


run1.exp.mat <- run1.exp.mat  %>% mutate_if(.,is.numeric,rank) %>%
  mutate_if(.,is.numeric,function(X){X*X})
# Edit ranking scales above!

PathExpTab   <- left_join(PathExpTab,run1.exp.mat, by = c("SYMBOL" = "Gene"))
PathExpTab   <- PathExpTab %>% group_by(., Pathway) %>% 
  summarise_if(is.numeric,mean, na.rm=T)


tab.order    <- colnames(PathExpTab)[-c(1)]
tab.order    <- tab.order[tab.order %in% left.out]
description  <- samp.tab %>% filter(., arrayID %in% tab.order) %>% 
                        select(.,arrayID,stage,study,outcome,gender,age)
tab.order   <- left_join(as_tibble(tab.order),
                        description, by= c("value"= "arrayID"))

tab.order$age[tab.order$age > 12 ] <- NA
GPL570.tab  <- tab.order %>% group_by(.,study,stage) %>% 
               summarise(., n = n())

write_csv(GPL570.tab,"Reports/GPL570samples.csv")


PathExpMat           <- as.matrix(PathExpTab[,-c(1)])
PathExpMat           <- PathExpMat[,left.out]
rownames(PathExpMat) <- pull(PathExpTab[,1])
PathExpMat           <- as.data.frame(PathExpMat)

plotMA(PathExpMat)
plotDensities(PathExpMat,legend = F)


###### Just redoing things with a Z-normalization

PathExpMat <- apply(PathExpMat, 2, function(X){ (X - mean(X))/sd(X)})
PathExpMat <- as.data.frame(PathExpMat)


#plotMA(PathExpMat2[,3:5])

            design <- model.matrix(~ tab.order$stage + 0,PathExpMat)
            colnames(design) <- c("AD","CH")
            fit <- lmFit(PathExpMat, design)
            cont.matrix <- makeContrasts(AD-CH, levels=design)
            fit2 <- contrasts.fit(fit, cont.matrix)
            fit2 <- eBayes(fit2, 0.01)
            tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
            
            head(tT,5)
            hist(tT$adj.P.Val)
            tT.select <- tT[ tT$adj.P.Val < 0.05 & abs(tT$logFC) > 0.1,]
           # tT.select <-             head(tT,20)

            nrow(tT.select)
            rownames(tT.select)
          
            write.csv(tT.select,"Reports/NewReports/DEpathways.csv")
           
            
            
            
             
            labs <- as.data.frame(tab.order[,2:6])
            rownames(labs) <- tab.order$value
            
            ann_colors = list(
              outcome = c(NONSURVIVOR ="#E41A1C", 
                          SURVIVOR ="#377EB8",
                          UNKNOWN ="#4DAF4A"),
              stage =  c(adult ="#E41A1C", 
                         child ="#377EB8")
            )
            
            pheatmap((as.matrix(PathExpMat[rownames((tT.select)),])),
                     cutree_rows = 3,
                     cutree_cols = 3,
                     show_rownames = F,
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation",
                     show_colnames =  F,
                     annotation_col = labs,clustering_method = "ward.D2",
                     scale = "row",
                     annotation_colors = ann_colors)
    

            #
            #
            #
            #   DO NOT WORK WITH THE REST NOW
            #
            #
            #
            #
            #
            #
            #Subtypes exist in children. Not much in adults. 
            
    
 
            
            rownames(tT.select)              
            
            pathprints <- readxl::read_xlsx("DATA/PathPrintRes.xlsx")
            pathprints <- pathprints[complete.cases(pull(pathprints,2)),]
            pathprints <- pull(pathprints,2)
            
            dif.paths  <- rownames(tT.select)    
            
            
            intersect(dif.paths,pathprints)
            ?phyper
           
pdf("images/NormalizedPathwayIntesities.pdf")
plotDensities(PathExpMat,legend = F,
              main = "Z-normalized Pathway Intensity per Sample")
dev.off()



library(ggfortify)

# PCA of normalized samples
GPL570.pca <- prcomp(t(PathExpMat[,tab.order$value]), center = TRUE, scale. = TRUE)
plot(GPL570.pca,type = "l")

set.seed(1)
# Does not seem to be affected by study ... Yoohoo?
autoplot(GPL570.pca, data = tab.order, colour = "study",shape = "stage")


colnames(PathExpMat) == tab.order$value

#Batch effect for children is not that much of a problem seemingly

#saveRDS(labs,"DATA/childTab.rds")
#saveRDS(child.mat,"DATA/childMat.rds")

### Looking at adults. There not much subclasses here. We'll see what to do.


## Putting both studies together          
#clustering highly differentially expressed pathways

library(mclust)

dif.mat <- (PathExpMat[rownames(tT.select),])
BIC <- mclustBIC(t(dif.mat))
plot(BIC)
summary(BIC)





# Try to run some PCA down the road.

seq.pca <- prcomp(t(PathExpMat[,tab.order$value]), center = TRUE, scale. = TRUE)
plot(seq.pca,type = "l")


library(ggfortify)

set.seed(1)

# Does not seem to be affected by study ... Yoohoo?
autoplot(seq.pca, data = tab.order, colour = "study")

library(pca3d)

pca3d::pca3d(seq.pca,group = tab.order$stage , show.centroids = T, new = T)


## The following is differential pathway analysis

# design <- model.matrix(~ tab.order$stage + 0,PathExpMat)
# colnames(design) <- c("AD","CH")
# fit <- lmFit(PathExpMat, design)
# cont.matrix <- makeContrasts(AD-CH, levels=design)
# fit2 <- contrasts.fit(fit, cont.matrix)
# fit2 <- eBayes(fit2, 0.01)
# tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
# 
# 
# tT.select <- (head(tT,40))
# 






library(caret) 
library(pROC)

short.set   <- PathExpMat[rownames(tT.select),]


set.seed(1)
validation_index <- createDataPartition(tab.order$stage, p=1, list=FALSE)
gset.test      <-  short.set[,-validation_index]
gene.exp.tot   <-  short.set[,validation_index]        


train.set     <- as_data_frame(t(gene.exp.tot)) 
train.set$sml <- as.factor(tab.order$stage)[validation_index]


tests.set     <- as_data_frame(t(gset.test)) 
tests.set.sml <- as.factor(tab.order$stage)[-validation_index]
### Lasso on the train data for further dimention reduction
control <- trainControl(method="repeatedcv",number = 5,repeats = 5,
                        summaryFunction = twoClassSummary, classProbs = TRUE)
set.seed(7)
fit.lasso <- train(sml ~ .,train.set,
                   method = "glmnet",
                   trControl = control,
                   metric = "ROC",
                   preProc = c("center", "scale"),
                   tuneGrid = expand.grid(.alpha = seq(.05, 1,length = 15),
                                          .lambda = c((1:5)/10)))

probs.glm <- predict(fit.lasso,tests.set,type = "prob")
ROC.glm   <- roc(tests.set.sml, probs.glm[,"child"], levels = levels(tests.set.sml))
print(ROC.glm$auc)

fit.lasso$finalModel$beta


aa <- varImp(fit.lasso, scale = FALSE)
hist(unlist(aa$importance))

lasso.feats <- rownames(aa$importance)
lasso.feats <- lasso.feats[which(aa$importance > 0.06)] 
lasso.feats <- gsub("`","",lasso.feats)
print(length(lasso.feats))  



short.set   <- PathExpMat[rownames(tT.select)[rownames(tT.select) %in% lasso.feats],]


labs <- as.data.frame(tab.order[,2:4])
rownames(labs) <- tab.order$value 
pheatmap(as.matrix(short.set),show_rownames = F,
         annotation_col = labs,clustering_method = "ward.D")





# Just a separation of data based on the platform. 
# The current study mostly focuses on GPL570 Platform

rm(list= ls())

library(Biobase)
library(GEOquery)
library(tidyverse)

# Loading expressions and information from original sepsis paper
# Selecting GPL570 platform only

exp.sets <- readRDS("DATA/AllSepsisDatasets.rds")
samp.tab <- readxl::read_xlsx("DATA/TrueSepsisDataUsedforFinalAnalysisMay2015.xlsx")

datasets <- samp.tab %>% filter(., platform == "GPL570") %>% 
                          select(.,study) %>% unique() %>% pull

name1    <- gsub("\\_.*","",names(exp.sets))
names(exp.sets)     <- gsub("^.*\\.","",name1)

exp.sets <- exp.sets[datasets]

# Getting an additional dataset to boost the number of normals
new.sets            <- getGEO("GSE35571",GSEMatrix =TRUE)
name1               <- gsub("\\_.*","",names(new.sets))
names(new.sets)     <- gsub("^.*\\.","",name1)


news.tab <- pData(new.sets[[1]])
type.ok  <- which(news.tab$`disease status:ch1` == "non-asthmatic")
ages.ok  <- which(as.numeric(news.tab$`age (yr):ch1`)  < 11.5)
alls.ok  <- intersect(type.ok,ages.ok)
alls.ok  <- rownames(news.tab[alls.ok,])




# GSE26378 contains age. Disease state which contains normal control
# GSE26440 contains age. Disease state which contains normal control
norm.sets <- exp.sets[1:2]


# Cleaning GSE26378 to find normal control children above 5.

set1 <- (pData(norm.sets[[1]]))
names(norm.sets[1])
set1$characteristics_ch1.2 <- as.character(set1$characteristics_ch1.2)
set1$characteristics_ch1.2 <- str_sub(set1$characteristics_ch1.2,start = 14)
ind1 <- which(as.numeric(set1$characteristics_ch1.2) > 5)
ind2 <- grep("control",set1$characteristics_ch1.1)

normals <- rownames(set1)[intersect(ind1,ind2)]


# There's something weird with this data
    test.this.set <- exprs(norm.sets[[2]])
    pdf("../plots /MABad", pointsize = 20)
    plotDensities(log(test.this.set), legend = F)
    dev.off()
    MA.q <- normalizeBetweenArrays(test.this.set, method="quantile")
    MA.q[MA.q < 10] <- 0
    MA.q[MA.q > 10] <- 1
    rem.ind <- which(rowSums((MA.q )) > 0)
    names(rem.ind)
    MA.q <- test.this.set[-(rem.ind),]
    plotDensities(log(MA.q),legend = F)
    limma::plotMA((MA.q[,6:7]))
    boxplot(log(MA.q))


normals.exp1 <- exprs(norm.sets[[1]])[,normals]
norm.ages.1  <- set1$characteristics_ch1.2[intersect(ind1,ind2)]
normal.info1 <- data.frame("ArrayID"= normals, "Age" = norm.ages.1,
                           "study"= names(norm.sets[1]) )

# Cleaning GSE26440 to find normal control children above 5.

set1         <- (pData(norm.sets[[2]]))

pdf("../plots /MABAD", pointsize = 20)
plotMA(log(normals.exp1))
dev.off()


ind1         <- which(as.numeric(set1$`age (years):ch1`) > 5)
ind2         <- grep("control",set1$`disease state:ch1`)

normals      <- rownames(set1)[intersect(ind1,ind2)]
normals.exp2 <- exprs(norm.sets[[2]])[,normals]

norm.ages.2  <- set1$`age (years):ch1`[intersect(ind1,ind2)]
normal.info2 <- data.frame("ArrayID"= normals, "Age" = norm.ages.2,
                           "study"= names(norm.sets[2]) )



# There's something weird with Wong  data
test.this.set <- exprs(new.sets[[1]])
MA.q <- normalizeBetweenArrays(test.this.set, method="quantile")


plotDensities((MA.q),legend = F)
boxplot(log(test.this.set))

?normalizeBetweenArrays

normals.exp3 <- exprs(new.sets[[1]])[,alls.ok]
norm.ages.3  <- new.sets[[1]]$`age (yr):ch1`[intersect(type.ok,ages.ok)]

normal.info3 <- data.frame("ArrayID"= alls.ok, "Age" = norm.ages.3,
                           "study"= names(new.sets[1]) )

pdf("../plots /MAGoodLog", pointsize = 20)
plotMA(log(normals.exp3),legend = F)
dev.off()
plotMA(log(normals.exp3))

# unifying probes across studies
uni.genes    <- intersect(rownames(normals.exp1),rownames(normals.exp3))
normals.exp1 <- normals.exp1[uni.genes,]
normals.exp2 <- normals.exp2[uni.genes,]
normals.exp3 <- normals.exp3[uni.genes,]


normal.info  <- do.call(rbind,list(normal.info1,normal.info2,normal.info3))

#Just checking if the gene names match
sum(rownames(normals.exp1) != rownames(normals.exp3))



normalChildren <- do.call(cbind,list(normals.exp1,normals.exp2,normals.exp3))




fData.run    <-  fData(exp.sets[[1]])
sym.table    <-  fData.run[,c("ID","Gene Symbol")]
sym.table    <-  sym.table[uni.genes,]
#rm(exp.sets)


sum(rownames(normalChildren) != sym.table$ID)

normalChildren      <- as_tibble(normalChildren)
normalChildren$Gene <- sym.table$`Gene Symbol`
normalChildren      <- normalChildren %>% group_by(.,Gene) %>%
                       summarise_if(.,is.numeric,mean)
normalChildren      <- normalChildren %>% filter(.,Gene != "")



saveRDS(normal.info,"DATA/normalChildrenMetaData")
saveRDS(normalChildren,"DATA/normalChildrenExpression.rds")


# Looking at online datasets for whole blood samples
# GSE35571: Childern Asthma/Non-asthma with age variable. Many samples
# GSE58667: Control Samples with age
# GSE14844: Many 12 year olds here
# GSE72439: A lot of samples of patients with Growth hormone deficiency
# GSE6575: Normal samples of children but no age variable on first glance
# GSE11755: Contains children sepsis, time series, and normal, no age variable
# GSE8121: Hector Wong dataset
# GSE67684 not related but look at this dataset...
#
# GSE9692, GSE13904, and GSE4607 have median age below 5. The age data exists
# but is not available online. Maybe Contact for them is Hector Wong?



clin <- readRDS("clinical.rds")
gex <- readRDS("expression.rds")

idx <- intersect(colnames(gex), clin$sample)
clin_p <- clin[clin$sample_id %in% idx, ]
stage <- as.vector(clin_p$stage)
stage_p <- rep(0,length(stage))
for (i in 1:length(stage)){
  if(stage[i]=='stage i'){stage_p[i]=1}
  else if(stage[i]=='stage ii'){stage_p[i]=2}
  else if(stage[i]=='stage iii'){stage_p[i]=3}
  else if(stage[i]=='stage iv'){stage_p[i]=4}
  else if(stage[i]=='stage v'){stage_p[i]=5}
}

gex_p <- gex[, colnames(gex) %in% idx]
gex_p <- na.omit(gex_p)
library(preprocessCore)
temp <- normalize.quantiles(gex_p)
colnames(temp) <- colnames(gex_p)
rownames(temp) <- rownames(gex_p)
gex_p <- temp

normalize <- function(x)
{
  return((x- min(x))/(max(x)-min(x)))
}
gex_p <- apply(gex_p, MARGIN=1, normalize)

data <- data.frame(stage_p, gex_p)
#µ¿ÀÏ
sum(data$stage_p == 1);sum(data$stage_p == 2);sum(data$stage_p == 3);
set.seed(123)
stage2.sam <- sample(rownames(data)[data$stage_p == 2], size=min(sum(data$stage_p==1)))
stage3.sam <- sample(rownames(data)[data$stage_p == 3], size=min(sum(data$stage_p==1)))
data <- data[rownames(data) %in% c(rownames(data)[data$stage_p == 1], stage2.sam, stage3.sam), ]

feature.size <- 2
aov.gene.pval <- c()
for(i in 2:(ncol(data))){
  aov.gene.pval <- c(aov.gene.pval, summary(aov(data$stage_p ~ data[, i]))[[1]][["Pr(>F)"]][1])
}
aov.idx <- order(aov.gene.pval, decreasing=FALSE)[1:feature.size] 
data.aov <- data[, c(1, 1+aov.idx)] 

library(caret)
data.aov$stage_p <- as.factor(data.aov$stage_p) 
partition <- createDataPartition(data.aov$stage_p, p=0.8, list=FALSE) 
train <- data.aov[partition,]
test <- data.aov[-partition,]
k <- expand.grid(k = c(1:89))
set.seed(123)
ctrl <- trainControl(method="cv", number=5)
knn.cv <- train(stage_p~., train, method="knn", trControl=ctrl,
                preProcess=c("center", "scale"), tuneGrid=k)
knn.cv$results
plot(knn.cv)

#/////////////////////////HW1

feature.size <- 2
aov.gene.pval <- c()
for(i in 2:(ncol(data))){
  aov.gene.pval <- c(aov.gene.pval, summary(aov(data$stage_p ~ data[, i]))[[1]][["Pr(>F)"]][1])
}
aov.idx <- order(aov.gene.pval, decreasing=FALSE)[1:feature.size] 
data.aov <- data[, c(1, 1+aov.idx)] 

library(caret)
data.aov$stage_p <- as.factor(data.aov$stage_p) 
partition <- createDataPartition(data.aov$stage_p, p=0.8, list=FALSE) 
train <- data.aov[partition,]
test <- data.aov[-partition,]
k <- expand.grid(k = c(1:89))
set.seed(123)
ctrl <- trainControl(method="cv", number=5)
knn.cv <- train(stage_p~., train, method="knn", trControl=ctrl,
                preProcess=c("center", "scale"), tuneGrid=k)
knn.cv$results
plot(knn.cv)
#///hw2
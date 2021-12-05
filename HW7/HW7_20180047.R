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
temp <- normalize.quantiles(gex_p)
colnames(temp) <- colnames(gex_p)
rownames(temp) <- rownames(gex_p)
gex_p <- temp
data <- data.frame(stage_p, t(gex_p))

set.seed(1)
stage2.sam <- sample(rownames(data)[data$stage_p == 2], size=min(sum(data$stage_p==1)))
stage3.sam <- sample(rownames(data)[data$stage_p == 3], size=min(sum(data$stage_p==1)))
data <- data[rownames(data) %in% c(rownames(data)[data$stage_p == 1], stage2.sam, stage3.sam), ]
temp <- normalize.quantiles(gex_p)
colnames(temp) <- colnames(gex_p)
rownames(temp) <- rownames(gex_p)
gex_p <- temp
data <- data.frame(stage_p, t(gex_p))


set.seed(1)
stage2.sam <- sample(rownames(data)[data$stage_p == 2], size=min(sum(data$stage_p==1)))
stage3.sam <- sample(rownames(data)[data$stage_p == 3], size=min(sum(data$stage_p==1)))
data <- data[rownames(data) %in% c(rownames(data)[data$stage_p == 1], stage2.sam, stage3.sam), ]
feature.size <- 10
aov.gene.pval <- c()
for(i in 2:(ncol(data))){
  aov.gene.pval <- c(aov.gene.pval, summary(aov(data$stage_p ~ data[, i]))[[1]][["Pr(>F)"]][1])
}
aov.idx <- order(aov.gene.pval, decreasing=FALSE)[1:feature.size] 
data.aov <- data[, c(1, 1+aov.idx)]
data.aov$stage_p <- as.factor(data.aov$stage_p) # stage_p numeric to factor for classification
partition <- createDataPartition(data.aov$stage_p, p=0.8, list=FALSE) # Training:Test=0.8:0.2
train <- data.aov[partition,]
test <- data.aov[-partition,]
#전처리끝
linear.tune <- tune(svm, stage_p~., data = train, kernel = "linear",
                    ranges = list(cost = c(0.001, 0.01, 0.1, 1, 10)))
summary(linear.tune)

linear.model <- linear.tune$best.model
summary(linear.model)

linear.pred <- predict(linear.model, test)
confusionMatrix(linear.pred, test$stage_p)

mean(linear.pred == test$stage_p)
#linear kernel
poly.tune <- tune(svm, stage_p~., data = train, kernel = "polynomial",
                  ranges = list(cost = c(0.001, 0.01, 0.1, 1, 10), gamma = c(0.01,0.5,1,2)))
poly.model <- poly.tune$best.model
summary(poly.model)

poly.pred <- predict(poly.model, test)
confusionMatrix(poly.pred, test$stage_p)
mean(poly.pred == test$stage_p)
#poly
rad.tune <- tune(svm, stage_p~., data = train, kernel = "radial",
                 ranges = list(cost = c(0.001, 0.01, 0.1, 1, 10, 100), gamma = c(0,001,0.005,0.01,0.5)))
summary(rad.tune)

rad.model <- rad.tune$best.model
summary(rad.model)

rad.pred <- predict(rad.model, test)
confusionMatrix(rad.pred, test$stage_p)
mean(rad.pred == test$stage_p)
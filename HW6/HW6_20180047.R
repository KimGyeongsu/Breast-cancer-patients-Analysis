clin <- readRDS("clinical.rds")
gex <- readRDS("expression.rds")

clin <- clin[clin$survival_time!= 0, ]
idx <- intersect(colnames(gex), clin$sample)
clin_p <- clin[clin$sample_id %in% idx, ]
surv_time <- clin_p$survival_time

gex_p <- gex[, colnames(gex) %in% idx]
gex_p <- na.omit(gex_p)
temp <- normalize.quantiles(gex_p)
colnames(temp) <- colnames(gex_p)
rownames(temp) <- rownames(gex_p)
gex_p <- temp

cor.p <- cor.coef <- cor.idx <- c()
for(i in 1:dim(gex_p)[1]){
  cor.p <- c(cor.p, cor.test(clin_p$survival_time, gex_p[i, ], method="pearson")$p.value)
  cor.coef <- c(cor.coef, cor.test(clin_p$survival_time, gex_p[i, ], method="pearson")$estimate)
}
cor.idx <- order(cor.p, decreasing=FALSE)[1:100]
gex_p.cor <- t(gex_p[cor.idx, ])
data.cor <- data.frame(surv_time, gex_p.cor)
#여기까지 전처리 / forward stepwise selection 시작

regfit <- regsubsets(surv_time ~ ., data = data.cor, method="forward", nvmax = 100)
reg.sum <- summary(regfit)
#adjusted R^2 / plot그리기 /
par(mfrow=c(1,2))
plot(reg.sum$adjr2, type = 'l', col = "red",
     xlab = "Number of Predictors", ylab = "Adjusted RSq")
points(reg.sum$adjr2, col = "red", pch = 16)
which(reg.sum$adjr2==max(reg.sum$adjr2))

#Ridge Lasso
x <- gex_p.cor
y <- surv_time

grid <- 10^seq(10, -2, length = 100)
ridge <- glmnet(x, y, alpha = 0, lambda = grid, standardize = TRUE)
set.seed(1)
train <- sample(1:nrow(x), nrow(x)/2)
test <- (-train)
y.test <- y[test]
ridge.mod <- glmnet(x[train,], y[train], alpha=0, lambda = grid)
cv.ridge <- cv.glmnet(x[train, ], y[train], alpha=0, nfolds=10)
plot(cv.ridge)
bestlambda <- cv.ridge$lambda.1se
bestlambda
out <- glmnet(x, y, alpha = 0, lambda = grid)
predict(out, type = "coefficients", s = bestlambda)

lasso <- glmnet(x, y, alpha = 1, lambda = grid, standardize = TRUE)
lasso.mod <- glmnet(x[train,], y[train], alpha=1, lambda = grid)
cv.lasso <- cv.glmnet(x[train, ], y[train], alpha=1, nfolds=10)
plot(cv.lasso)
bestlambda <- cv.lasso$lambda.1se
bestlambda
out <- glmnet(x, y, alpha = 1, lambda = grid)
predict(out, type = "coefficients", s = bestlambda)
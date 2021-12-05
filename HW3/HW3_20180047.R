clinical<-readRDS("clinical.rds")
expression<-readRDS("expression.rds")
mutation<-readRDS("mutation.rds")
#데이터 뽑아오기
sample_id_filter<-intersect(colnames(expression),
                            intersect(unique(mutation$sample_id),
                                      clinical$sample_id[which(!is.na(clinical$subtype))]))
clinical_filter<-dplyr::filter(clinical,sample_id %in% sample_id_filter)
mutation_filter<-dplyr::filter(mutation, sample_id %in% sample_id_filter)
expression_filter<-expression[,sample_id_filter]
#NA subtype 제외.
cor.coef<-cor.pvalue<-c()
for(i in 1:dim(expression_filter)[1]){
  survival.time<-clinical_filter$survival_time
  expression<-expression_filter[i,]
  cor.result<-cor.test(survival.time,expression, method="pearson")
  cor.coef<-c(cor.coef, cor.result$estimate)
  cor.pvalue<-c(cor.pvalue, cor.result$p.value)
}
gene.correlation_sort<-row.names(expression_filter)[order(cor.coef,decreasing=TRUE)]
gene.correlation_sort_sort<-gene.correlation_sort[which(cor.pvalue<0.01)]
head(gene.correlation_sort_sort, 10)

colnames.data<-gene.correlation_sort_sort[c(1:10)]
data<-as.data.frame(t(expression_filter[colnames.data,]))
colnames(data)<-colnames.data
data$survival<-clinical_filter$survival_time
res.simple<-lm(survival~.,data=data)
summary(res.simple)
#C4orf7 제외(0.78174)
colnames.data<-gene.correlation_sort_sort[c(1,2,3,4,5,6,7,8,10)]
data<-as.data.frame(t(expression_filter[colnames.data,]))
colnames(data)<-colnames.data
data$survival<-clinical_filter$survival_time
res.simple<-lm(survival~.,data=data)
summary(res.simple)
#RIMS3 제외(0.76697)
colnames.data<-gene.correlation_sort_sort[c(1,2,3,4,5,6,8,10)]
data<-as.data.frame(t(expression_filter[colnames.data,]))
colnames(data)<-colnames.data
data$survival<-clinical_filter$survival_time
res.simple<-lm(survival~.,data=data)
summary(res.simple)
#STEAP3 제외(0.618799)
colnames.data<-gene.correlation_sort_sort[c(1,2,3,4,6,8,10)]
data<-as.data.frame(t(expression_filter[colnames.data,]))
colnames(data)<-colnames.data
data$survival<-clinical_filter$survival_time
res.simple<-lm(survival~.,data=data)
summary(res.simple)
#MALT1 제외(0.409867)
colnames.data<-gene.correlation_sort_sort[c(1,2,3,4,8,10)]
data<-as.data.frame(t(expression_filter[colnames.data,]))
colnames(data)<-colnames.data
data$survival<-clinical_filter$survival_time
res.simple<-lm(survival~.,data=data)
summary(res.simple)
#KRT5 제외(0.264202)
colnames.data<-gene.correlation_sort_sort[c(2,3,4,8,10)]
data<-as.data.frame(t(expression_filter[colnames.data,]))
colnames(data)<-colnames.data
data$survival<-clinical_filter$survival_time
res.simple<-lm(survival~.,data=data)
summary(res.simple)
#RBP5 제외(0.111280)
colnames.data<-gene.correlation_sort_sort[c(2,3,4,10)]
data<-as.data.frame(t(expression_filter[colnames.data,]))
colnames(data)<-colnames.data
data$survival<-clinical_filter$survival_time
res.simple<-lm(survival~.,data=data)
summary(res.simple)
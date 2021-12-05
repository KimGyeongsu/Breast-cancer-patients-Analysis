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
boxplot(survival_time~stage,data=clinical_filter,xlab="Clinical Stage",ylab="Survival Time (days)")
#그림그리기
t.test(clinical_filter$survival_time[which(clinical_filter$stage=="stage i")],
       clinical_filter$survival_time[which(clinical_filter$stage=="stage ii")],'greater')$p.value
#one-sided로 stage1,stage2 
t.test(clinical_filter$survival_time[which(clinical_filter$stage=="stage i")],
       clinical_filter$survival_time[which(clinical_filter$stage=="stage iii")],'greater')$p.value
#one-sided로 stage1,stage3
t.test(clinical_filter$survival_time[which(clinical_filter$stage=="stage ii")],
       clinical_filter$survival_time[which(clinical_filter$stage=="stage iii")],'greater')$p.value
#one-sdied로 stage2,stage3
stage1<-as.character(clinical_filter$sample_id[which(clinical_filter$stage=="stage i")])
stage2<-as.character(clinical_filter$sample_id[which(clinical_filter$stage=="stage ii")])
stage3<-as.character(clinical_filter$sample_id[which(clinical_filter$stage=="stage iii")])
expression.stage3<-expression_filter[,stage3]
expression.stage2<-expression_filter[,stage2]
expression.stage1<-expression_filter[,stage1]

pvalues<-sapply(c(1:dim(expression_filter)[1]),
                FUN=function(k){
                  pval<-t.test(expression.stage1[k,],expression.stage2[k,])$p.value
                  return(pval)
                })
pvalues.adj<-p.adjust(pvalues,method = "fdr")
length(pvalues[which(pvalues<0.05)])
length(pvalues.adj[which(pvalues.adj<0.05)])
#stage 1과 stage 2 // stage 2와 3 // stage 과 3에 대해서도 동일하게 진행
clinical_filter<-dplyr::filter(clinical, subtype!="NA")
clinical_filter<-dplyr::filter(clinical_filter, subtype!="Normal-like")
boxplot(survival_time~subtype,data=clinical_filter,xlab="Subtype",ylab="Survival Time (days)")

anova.result <- aov(survival_time~subtype,data=clinical_filter)
summary_ANOVA <- summary(anova.result)
summary_ANOVA
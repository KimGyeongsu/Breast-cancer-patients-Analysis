clinical<-readRDS("clinical.rds")
expression<-readRDS("expression.rds")
mutation<-readRDS("mutation.rds")
#데이터 뽑아오기
sample_id_filter<-intersect(colnames(expression),
                            intersect(unique(mutation$sample_id),
                                      clinical$sample_id[which(!is.na(clinical$subtype))]))
clinical_filter<-dplyr::filter(clinical,sample_id %in% sample_id_filter)
mutation_filter<-dplyr::filter(mutation, sample_id %in% sample_id_filter)
expression_filter<-expression[,as.character(clinical_filter$sample_id)]
clinical_filter$stage <- factor(clinical_filter$stage)
vital_status <- 1*(clinical_filter$vital_status == 0)
survival_time <- as.numeric(clinical_filter$survival_time)
su <- survival::Surv(survival_time, vital_status)
fit.survival<-survfit(su~stage, data=clinical_filter)
fit.survival
ggsurvplot(fit.survival, data=clinical_filter, pval=TRUE)
fit.cox<-coxph(su~stage, data=clinical_filter)
fit.cox
ggforest(fit.cox, data=clinical_filter)
clinical_filter$subtype <- factor(clinical_filter$subtype)
fit.survival<-survfit(su~subtype, data=clinical_filter)
fit.survival
ggsurvplot(fit.survival, data=clinical_filter, pval=TRUE)
fit.cox<-coxph(su~subtype, data=clinical_filter)
fit.cox
ggforest(fit.cox, data=clinical_filter)
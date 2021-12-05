clinical<-readRDS("clinical.rds")
mutation<-readRDS("mutation.rds")
clinical_filter<-dplyr::filter(clinical, subtype!="Normal-like")
subtype<-unique(clinical_filter$subtype)
key<-"Luminal A"
clinical_filter_key<-dplyr::filter(clinical_filter, subtype==key)
survived_patients<-length(which(clinical_filter_key$survival_time>1825))
total_patients<-length(which(clinical_filter_key$survival_time>=0))
survival_rate<-survived_patients/total_patients
print(paste("This is the result which subtype is ", key, ".  Total number of the survived patients are ", survived_patients, ". The number of all patients are ", total_patients, ". In conclusion, 5-year survival probability of breast cancer patients which subtype is ", key, "is ", round(survival_rate, 3), "."))
key<-"Basal-like"
clinical_filter_key<-dplyr::filter(clinical_filter, subtype==key)
survived_patients<-length(which(clinical_filter_key$survival_time>1825))
total_patients<-length(which(clinical_filter_key$survival_time>=0))
survival_rate<-survived_patients/total_patients
print(paste("This is the result which subtype is ", key, ".  Total number of the survived patients are ", survived_patients, ". The number of all patients are ", total_patients, ". In conclusion, 5-year survival probability of breast cancer patients which subtype is ", key, "is ", round(survival_rate, 3), "."))
key<-"Luminal B"
clinical_filter_key<-dplyr::filter(clinical_filter, subtype==key)
survived_patients<-length(which(clinical_filter_key$survival_time>1825))
total_patients<-length(which(clinical_filter_key$survival_time>=0))
survival_rate<-survived_patients/total_patients
print(paste("This is the result which subtype is ", key, ".  Total number of the survived patients are ", survived_patients, ". The number of all patients are ", total_patients, ". In conclusion, 5-year survival probability of breast cancer patients which subtype is ", key, "is ", round(survival_rate, 3), "."))
key<-"HER2-enriched"
clinical_filter_key<-dplyr::filter(clinical_filter, subtype==key)
survived_patients<-length(which(clinical_filter_key$survival_time>1825))
total_patients<-length(which(clinical_filter_key$survival_time>=0))
survival_rate<-survived_patients/total_patients
print(paste("This is the result which subtype is ", key, ".  Total number of the survived patients are ", survived_patients, ". The number of all patients are ", total_patients, ". In conclusion, 5-year survival probability of breast cancer patients which subtype is ", key, "is ", round(survival_rate, 3), "."))
key <- "Luminal A"
clinical_key<-dplyr::filter(clinical,subtype==key)         
sample_id_filter<-clinical_key$sample_id         
mutation_key<-dplyr::filter(mutation, sample_id %in% sample_id_filter)
gene_key_raw<-(mutation_key$Hugo_Symbol)
gene_key<-unique(gene_key_raw)
cnt<-rep(0,length(gene_key))
for(i in 1:length(gene_key_raw)){
  for(j in 1 : length(gene_key)){
    if(gene_key_raw[i]==gene_key[j]){
      cnt[j]<-cnt[j]+1
      break
    }
  }
}
print(paste("This is the result which subtype is ",key,". "))
for(i in 1:10){
  temp<-gene_key[which.max(cnt)]
  print(paste("Total number of patients is ",length(sample_id_filter),". The number of patients with ",temp," mutation is ",max(cnt),". In conclusion, conditional probability is ",round(max(cnt)/length(sample_id_filter),3),"."))
  cnt[which.max(cnt)]<-0
}
key <- "Basal-like"
clinical_key<-dplyr::filter(clinical,subtype==key)         
sample_id_filter<-clinical_key$sample_id         
mutation_key<-dplyr::filter(mutation, sample_id %in% sample_id_filter)
gene_key_raw<-(mutation_key$Hugo_Symbol)
gene_key<-unique(gene_key_raw)
cnt<-rep(0,length(gene_key))
for(i in 1:length(gene_key_raw)){
  for(j in 1 : length(gene_key)){
    if(gene_key_raw[i]==gene_key[j]){
      cnt[j]<-cnt[j]+1
      break
    }
  }
}
print(paste("This is the result which subtype is ",key,". "))
for(i in 1:10){
  temp<-gene_key[which.max(cnt)]
  print(paste("Total number of patients is ",length(sample_id_filter),". The number of patients with ",temp," mutation is ",max(cnt),". In conclusion, conditional probability is ",round(max(cnt)/length(sample_id_filter),3),"."))
  cnt[which.max(cnt)]<-0
}
key <- "Luminal B"
clinical_key<-dplyr::filter(clinical,subtype==key)         
sample_id_filter<-clinical_key$sample_id         
mutation_key<-dplyr::filter(mutation, sample_id %in% sample_id_filter)
gene_key_raw<-(mutation_key$Hugo_Symbol)
gene_key<-unique(gene_key_raw)
cnt<-rep(0,length(gene_key))
for(i in 1:length(gene_key_raw)){
  for(j in 1 : length(gene_key)){
    if(gene_key_raw[i]==gene_key[j]){
      cnt[j]<-cnt[j]+1
      break
    }
  }
}
print(paste("This is the result which subtype is ",key,". "))
for(i in 1:10){
  temp<-gene_key[which.max(cnt)]
  print(paste("Total number of patients is ",length(sample_id_filter),". The number of patients with ",temp," mutation is ",max(cnt),". In conclusion, conditional probability is ",round(max(cnt)/length(sample_id_filter),3),"."))
  cnt[which.max(cnt)]<-0
}
key <- "HER2-enriched"
clinical_key<-dplyr::filter(clinical,subtype==key)         
sample_id_filter<-clinical_key$sample_id         
mutation_key<-dplyr::filter(mutation, sample_id %in% sample_id_filter)
gene_key_raw<-(mutation_key$Hugo_Symbol)
gene_key<-unique(gene_key_raw)
cnt<-rep(0,length(gene_key))
for(i in 1:length(gene_key_raw)){
  for(j in 1 : length(gene_key)){
    if(gene_key_raw[i]==gene_key[j]){
      cnt[j]<-cnt[j]+1
      break
    }
  }
}
print(paste("This is the result which subtype is ",key,". "))
for(i in 1:10){
  temp<-gene_key[which.max(cnt)]
  print(paste("Total number of patients is ",length(sample_id_filter),". The number of patients with ",temp," mutation is ",max(cnt),". In conclusion, conditional probability is ",round(max(cnt)/length(sample_id_filter),3),"."))
  cnt[which.max(cnt)]<-0
}
gex <- readRDS("expression.rds")
clin <- readRDS("clinical.rds")

DATA <-na.omit(t(gex))
hc_average <- hclust(dist(DATA), method="average")
{plot(hc_average, main="Average Linkage", xlab="", sub="", cex=.9)
  rect.hclust(hc_average, k=4)}

hc_complete <- hclust(dist(DATA), method="complete")
{plot(hc_complete, main="Complete Linkage", xlab="", sub="", cex=.9)
  rect.hclust(hc_complete, k=4)}


fviz_nbclust(DATA, hcut, method="silhouette", hc_method="complete")

DATA_sam <- data.frame(as.character(rownames(DATA)), DATA)
colnames(DATA_sam)[1] <- "sample_id"
stage <- as.vector(clin$stage)
stage_p <- rep(0,length(stage))
for (i in 1:length(stage)){
  if(stage[i]=='stage i'){stage_p[i]=1}
  else if(stage[i]=='stage ii'){stage_p[i]=2}
  else if(stage[i]=='stage iii'){stage_p[i]=3}
  else if(stage[i]=='stage iv'){stage_p[i]=4}
  else if(stage[i]=='stage v'){stage_p[i]=5}
}
clin$stage <- stage_p
DATA_clin <- merge(clin[,c("sample_id","stage")], DATA_sam, by = "sample_id", all = FALSE)
data <- DATA_clin
data$stage <- NULL
data$sample_id <- NULL
data_hc <- cutree(hclust(dist(data), method="complete"), 3)
print(cluster_similarity(DATA_clin$stage, data_hc, similarity="jaccard"))


PCA_DATA <- t(gex)
PCA <- prcomp(na.omit(PCA_DATA), scale=T, center=T)

DATA <-PCA$x[,1:2]
fviz_nbclust(DATA, hcut, method="silhouette", hc_method="complete")
DATA_sam <- data.frame(as.character(rownames(DATA)), DATA)
colnames(DATA_sam)[1] <- "sample_id"
stage <- as.vector(clin$stage)
stage_p <- rep(0,length(stage))
for (i in 1:length(stage)){
  if(stage[i]=='stage i'){stage_p[i]=1}
  else if(stage[i]=='stage ii'){stage_p[i]=2}
  else if(stage[i]=='stage iii'){stage_p[i]=3}
  else if(stage[i]=='stage iv'){stage_p[i]=4}
  else if(stage[i]=='stage v'){stage_p[i]=5}
}
clin$stage <- stage_p
DATA_clin <- merge(clin[,c("sample_id","stage")], DATA_sam, by = "sample_id", all = FALSE)
data <- DATA_clin
data$stage <- NULL
data$sample_id <- NULL
data_hc <- cutree(hclust(dist(data), method="complete"), 2)
print(cluster_similarity(DATA_clin$stage, data_hc, similarity="jaccard"))

DATA <-PCA$x[,1:3]
fviz_nbclust(DATA, hcut, method="silhouette", hc_method="complete")
DATA_sam <- data.frame(as.character(rownames(DATA)), DATA)
colnames(DATA_sam)[1] <- "sample_id"
stage <- as.vector(clin$stage)
stage_p <- rep(0,length(stage))
for (i in 1:length(stage)){
  if(stage[i]=='stage i'){stage_p[i]=1}
  else if(stage[i]=='stage ii'){stage_p[i]=2}
  else if(stage[i]=='stage iii'){stage_p[i]=3}
  else if(stage[i]=='stage iv'){stage_p[i]=4}
  else if(stage[i]=='stage v'){stage_p[i]=5}
}
clin$stage <- stage_p
DATA_clin <- merge(clin[,c("sample_id","stage")], DATA_sam, by = "sample_id", all = FALSE)
data <- DATA_clin
data$stage <- NULL
data$sample_id <- NULL
data_hc <- cutree(hclust(dist(data), method="complete"), 3)
print(cluster_similarity(DATA_clin$stage, data_hc, similarity="jaccard"))

DATA <-PCA$x[,1:4]
fviz_nbclust(DATA, hcut, method="silhouette", hc_method="complete")
DATA_sam <- data.frame(as.character(rownames(DATA)), DATA)
colnames(DATA_sam)[1] <- "sample_id"
stage <- as.vector(clin$stage)
stage_p <- rep(0,length(stage))
for (i in 1:length(stage)){
  if(stage[i]=='stage i'){stage_p[i]=1}
  else if(stage[i]=='stage ii'){stage_p[i]=2}
  else if(stage[i]=='stage iii'){stage_p[i]=3}
  else if(stage[i]=='stage iv'){stage_p[i]=4}
  else if(stage[i]=='stage v'){stage_p[i]=5}
}
clin$stage <- stage_p
DATA_clin <- merge(clin[,c("sample_id","stage")], DATA_sam, by = "sample_id", all = FALSE)
data <- DATA_clin
data$stage <- NULL
data$sample_id <- NULL
data_hc <- cutree(hclust(dist(data), method="complete"), 3)
print(cluster_similarity(DATA_clin$stage, data_hc, similarity="jaccard"))

DATA <-PCA$x[,1:5]
fviz_nbclust(DATA, hcut, method="silhouette", hc_method="complete")
DATA_sam <- data.frame(as.character(rownames(DATA)), DATA)
colnames(DATA_sam)[1] <- "sample_id"
stage <- as.vector(clin$stage)
stage_p <- rep(0,length(stage))
for (i in 1:length(stage)){
  if(stage[i]=='stage i'){stage_p[i]=1}
  else if(stage[i]=='stage ii'){stage_p[i]=2}
  else if(stage[i]=='stage iii'){stage_p[i]=3}
  else if(stage[i]=='stage iv'){stage_p[i]=4}
  else if(stage[i]=='stage v'){stage_p[i]=5}
}
clin$stage <- stage_p
DATA_clin <- merge(clin[,c("sample_id","stage")], DATA_sam, by = "sample_id", all = FALSE)
data <- DATA_clin
data$stage <- NULL
data$sample_id <- NULL
data_hc <- cutree(hclust(dist(data), method="complete"), 4)
print(cluster_similarity(DATA_clin$stage, data_hc, similarity="jaccard"))
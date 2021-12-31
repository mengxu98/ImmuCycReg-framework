

library(ConsensusClusterPlus)
library(limma)
library(NMF)
library(pheatmap)
load('luad_mrna.Rdata')
dataset_raw <- as.matrix(data_derived)
dataset_raw <- log(dataset_raw+1,10)
gene_no <- 30
mads <- apply(dataset_raw, 1, mad)
dataset_all<-dataset_raw[rev(order(mads)),]
dataset <- dataset_all[1:gene_no,]
max_cluster_num <- 8
title <- '../results'
results <- ConsensusClusterPlus(dataset, maxK =max_cluster_num,
                                reps = 50, pItem = 0.8,
                                pFeature = 0.8,  
                                clusterAlg = "pam", 
                                distance = "euclidean",
                                seed = 135792468,
                                title = title,
                                corUse = 'complete.obs',
                                plot = "pdf")
icl = calcICL(results,title=title,plot="pdf")
for(i in 2:max_cluster_num) {
  write.table(results[[i]][["consensusClass"]],
              file=paste(title,paste(paste("sample_cluster_",as.character(i),sep=""),".csv",sep=""),sep="//"),
              na = "",
              col.names = FALSE,
              sep = ",")
  }
sample_label <- read.table(paste(title,'/sample_cluster_4.csv',sep=''),
                           header = F,
                           sep=',',
                           check.names=FALSE)
cluster1<-sample_label[which(sample_label$V2==1),1]
cluster2<-sample_label[which(sample_label$V2==2),1]
cluster3<-sample_label[which(sample_label$V2==3),1]
cluster4<-sample_label[which(sample_label$V2==4),1]
tcga_cluster1<-dataset_all[,cluster1]
tcga_cluster2<-dataset_all[,cluster2]
tcga_cluster3<-dataset_all[,cluster3]
tcga_cluster4<-dataset_all[,cluster4]
group1_no<-dim(tcga_cluster1)[2]
group2_no<-dim(tcga_cluster2)[2]
group3_no<-dim(tcga_cluster3)[2]
group4_no<-dim(tcga_cluster4)[2]
group1<-rep(c('C1'),each=group1_no)
group2<-rep(c('C2'),each=group2_no)
group3<-rep(c('C3'),each=group3_no)
group4<-rep(c('C4'),each=group4_no)
group_label<-c(group1,group2,group3,group4)
annotation_c <- data.frame(group_label)
sorted_samples<-cbind(tcga_cluster1,tcga_cluster2)
sorted_samples<-cbind(sorted_samples,tcga_cluster3)
sorted_samples<-cbind(sorted_samples,tcga_cluster4)
rownames(annotation_c) <- colnames(sorted_samples)
pheatmap(sorted_samples, 
         #scale='row',
         cluster_row =TRUE,cluster_col = FALSE,
         show_colnames = F,
         annotation_col =annotation_c)
smallset <- sorted_samples
bigset <- dataset_raw
identical_flag <- TRUE
for(i in 1:nrow(smallset)) {
  for(j in 1:ncol(smallset)) {
    row_name <- rownames(smallset)[i]
    col_name <- colnames(smallset)[j]
    if(smallset[i,j] != bigset[row_name,col_name]) {
      print('not equal!! wrong!')
      identical_flag <- FALSE
    }
    if(smallset[i,j]==bigset[row_name,col_name]) {
      print(i)
    }
  }
}
if(identical_flag) {
  print('data is OK')
}
dataset <- as.matrix(data_derived)
dataset <- log(dataset+1,2)
gene_no <- 30
mads <- apply(dataset, 1, mad)
dataset <- dataset[rev(order(mads))[1:gene_no],]
res <- nmf(dataset,2:6,nrun=10,seed=135792468)
plot(res)
plot(2:6,res$measures$cophenetic, type="b", col="purple")
aheatmap(dataset)
consensusmap(res)
res_4 <- nmf(dataset, 4, nrun=10,seed=135792468)
w <- basis(res_4)
dim(w)
h <- coef(res_4)
dim(h)
opar <- par(mfrow=c(1,2))
basismap(res_4, subsetRow=TRUE)
basismap(res_4)
coefmap(res_4)
par(opar)
consensusmap(res_4)
consensusmap(res_4, annCol=dataset, tracks=NA)


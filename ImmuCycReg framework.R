

rm(list=ls())
library(psych)
library(rtracklayer)
library(svglite)
ASPAR_Predicted <- read.table(' ASPAR Predicted Transcription Factor Targets.txt.gz',
                              header = T,
                              check.names=FALSE)
encode_tf <- read.table(' ENCODE Transcription Factor Targets.txt.gz',
                        header = T,
                        check.names=FALSE)
CHEA_tf <- read.table(' CHEA Transcription Factor Targets.txt.gz',
                      header = T,
                      check.names=FALSE)
MotifMap_Predicted <- read.table(' MotifMap Predicted Transcription Factor Targets.txt.gz',
                                 header = T,
                                 sep = '\t',
                                 check.names=FALSE)
TRANSFAC_Curated <- read.table(' TRANSFAC Curated Transcription Factor Targets.txt.gz',
                               header = T,
                               check.names=FALSE)
TRANSFAC_Predicted <- read.table(' TRANSFAC Predicted Transcription Factor Targets.txt.gz',
                                 header = T,
                                 check.names=FALSE)
load('luad_mrna.Rdata')
negative_gene <- c('CXCR5','CCL20','CCL22','CCL28','EZH2','NT5E','LAG3','RAET1G')
source('myFunctions.R')
target_gene_list <- read.table('IRgenes.txt',header = T)
Regulation_result <- c()
Regulation_result_4 <- c()
JustPositiveGene <- c()
JustNegativeGene <- c()
immune_score_results <- c()
RegulationIndex <- c()
res_frame_all <- c()
for (k in 1:length(t(target_gene_list))) {
  load('tcga_gtex_luad_mrna.RData')
  load('luad_cnv.Rdata')
  load('peak_luad.Rdata')
  load('genes_adj_peak.Rdata')
  load('IRgene_tcga_mrna.Rdata')
  load('geneinfo_df.Rdata')
  target_gene <- t(target_gene_list)[,k]
  genes_temp <- as.vector(genes_adj_peak[,1])
  genes_list <- strsplit(genes_temp,',')
  peak_around_gene <- c()
  row_names <- row.names(genes_adj_peak)
  for(i in 1:nrow(genes_adj_peak)){
    if(target_gene%in%genes_list[[i]])
      peak_around_gene <- c(peak_around_gene,row_names[i])
  }
  Y <- t(allgene_tcga_mrna[target_gene,])
  sample_name <- colnames(allgene_tcga_mrna)
  X <- t(peak_luad[peak_around_gene,sample_name])
  peak_corr_r <- c()
  peak_corr_p <- c()
  for(i in 1:length(peak_around_gene)){
    peak_corr <- corr.test(X[,i],
                           Y,
                           method = "spearman",
                           adjust = "none")
    peak_corr_r <- c(peak_corr_r,peak_corr$r)
    peak_corr_p <- c(peak_corr_p,peak_corr$p)
  }
  deg.data <- data.frame(Symbol=peak_around_gene,
                         corr=peak_corr_r,
                         pval=peak_corr_p)
  deg.data$logP <- -log10(deg.data$pval)
  deg.data$Group = "not-significant"
  deg.data$Group[which((deg.data$pval<0.05)&(deg.data$corr > 0.1))] = "up-regulated"
  deg.data$Group[which((deg.data$pval<0.05)&(deg.data$corr < -0.1))] = "down-regulated"
  deg.data$Label = ''
  deg.data <- deg.data[order(deg.data$pval),]
  up.peak <- head(deg.data$Symbol[which(deg.data$Group == 'up-regulated')],20)
  down.peak <- head(deg.data$Symbol[which(deg.data$Group=='down-regulated')],20)
  not.sig.peak <- head(deg.data$Symbol[which(deg.data$Group == 'not-significant')],20)
  deg.top.peak <- c(as.character(up.peak),as.character(down.peak),as.character(not.sig.peak))
  deg.data$Label[match(deg.top.peak,deg.data$Symbol)] <- deg.top.peak
  high_cor_peak <- as.character(up.peak[1])
  down_cor_peak <- as.character(down.peak[1])
  temp_res <- geneinfo_df[which(geneinfo_df[, 'gene_name']==target_gene),c('seqnames', 'gene_id', 'start', 'end', 'width')]
  order_index <- order(temp_res[,'start'],decreasing = F)
  most_start <- temp_res[order_index[1],'start']
  most_end <- temp_res[order_index[1],'end']
  most_mid <- (most_start+most_end) / 2
  pos_col_name <- c('chrom','chromStart','chromEnd','strand')
  peaks_pos_temp <- genes_adj_peak[peak_around_gene,pos_col_name]
  peaks_pos_temp$mid <- (peaks_pos_temp$chromStart+peaks_pos_temp$chromEnd)/2
  peaks_pos_temp$mid_distance <- peaks_pos_temp$mid-most_mid
  peaks_pos_temp$start_distance <- peaks_pos_temp$chromEnd-most_start
  distance_thres<- -8000
  nearest_peaks <- row.names(peaks_pos_temp[which(peaks_pos_temp$start_distance<0 &peaks_pos_temp$start_distance>distance_thres& peaks_pos_temp$chromEnd<most_start),])
  no_peaks <- length(nearest_peaks)
  nearest_peaks <- nearest_peaks[no_peaks]
  genes_adj_peak[nearest_peaks, pos_col_name]
  top_k <- min(1,length(peak_around_gene))
  target_sample <- 'TCGA-86-A4D0-01'
  order_index <- order(peak_luad[peak_around_gene, target_sample],decreasing = T)
  order_index <- order(peak_luad[peak_around_gene, target_sample],decreasing = T)
  order_peak_id <- peak_around_gene[order_index]
  highest_score_peak <- order_peak_id[1:top_k]
  candidate_peak_id <- c(high_cor_peak)
  candidate_peak_id <- na.omit(candidate_peak_id)
  if (length(t(candidate_peak_id)) == 0) {
    candidate_peak_id <- c(down_cor_peak)
  }
  candidate_peak_id <- na.omit(candidate_peak_id)
  gene_list <- read.table(paste('TFs_list/',target_gene,'.txt',sep = ''),
                          row.names=1)
  sample_list <-read.table('peak_luad_sample.csv',
                           row.names=1,
                           header = T)
  row.names(sample_list)<-substr(row.names(sample_list),1,15)
  col_matrix<- colnames(raw_tcga_mrna)%in%row.names(sample_list)
  allgene_tcga_mrna <- raw_tcga_mrna[,which(col_matrix)]
  row_matrix_gtex <- row.names(raw_gtex_mrna)%in%row.names(gene_list)
  row_matrix_tcga <- row.names(allgene_tcga_mrna)%in%row.names(gene_list)
  gtex_temp <- raw_gtex_mrna[which(row_matrix_gtex),]
  tcga_temp <- allgene_tcga_mrna[which(row_matrix_tcga),]
  exsit_gene_list <- row.names(tcga_temp)
  high_exp_sample <- c()
  down_exp_sample <- c()
  for(sample_i in 1:ncol(allgene_tcga_mrna)){
    target_sample <- colnames(allgene_tcga_mrna)[sample_i]
    t_res <- t.test(raw_gtex_mrna[target_gene,],mu=allgene_tcga_mrna[target_gene,target_sample])
    if(t_res$statistic < 0 & t_res$p.value < 0.01) {
      print('target gene higher expression than gtex samples !')
      high_exp_sample <- c(high_exp_sample,target_sample)
    }else{
      print('target gene lower expression than gtex samples !')
      down_exp_sample <- c(down_exp_sample,target_sample)
    }
  }
  results_summary_up <- list()
  if (length(high_exp_sample) > 0) {
    for (sample_i in 1:length(high_exp_sample)) {
      target_sample <- high_exp_sample[sample_i]
      up_reg_gene <- c()
      for (i in 1:length(exsit_gene_list)) {
        t_res_temp <- t.test(gtex_temp[i,],mu = tcga_temp[i,target_sample])
        if (t_res_temp$statistic < 0 & t_res_temp$p.value < 0.01) {
          up_reg_gene<-c(up_reg_gene,exsit_gene_list[i])
        }
      }
      TF_up_genes <- up_reg_gene[which(raw_tcga_cnv[up_reg_gene,target_sample]>0)]
      results_summary_up[[target_sample]] <- TF_up_genes
    }
  }
  results_summary_down <- list()
  if (length(down_exp_sample) > 0) {
    for (sample_i in 1:length(down_exp_sample)) {
      target_sample <- down_exp_sample[sample_i]
      up_reg_gene <- c()
      for (i in 1:length(exsit_gene_list)) {
        t_res_temp <- t.test(gtex_temp[i,],mu = tcga_temp[i,target_sample])
        if (t_res_temp$statistic < 0 & t_res_temp$p.value < 0.01) {
          up_reg_gene<-c(up_reg_gene,exsit_gene_list[i])
        }
      }
      TF_up_genes <- up_reg_gene[which(raw_tcga_cnv[up_reg_gene,target_sample]>0)]
      results_summary_down[[target_sample]] <- TF_up_genes
    }
  }
  TFs <- c()
  type_inter <- c()
  target_gene_inter <- c()
  sample_level <- c()
  if (length(results_summary_up)>0) {
    for(i in 1:length(results_summary_up)){
      j=1
      while(j <= length(results_summary_up[[i]])){
        TFs <- c(TFs,results_summary_up[[i]][j]);
        type_inter <- c(type_inter,'1');
        sample_level <- c(sample_level,names(results_summary_up)[i]);
        target_gene_inter <- c(target_gene_inter,target_gene);
        j=j+1
      }
    }
  }else{
    JustPositiveGene <- c(JustPositiveGene,target_gene)
  }
  if (length(results_summary_down)>0) {
    for(i in 1:length(results_summary_down)){
      j=1
      while(j <= length(results_summary_down[[i]])){
        TFs <- c(TFs,results_summary_down[[i]][j])
        type_inter <- c(type_inter,'2')# type 2
        sample_level <- c(sample_level,names(results_summary_down)[i])
        target_gene_inter <- c(target_gene_inter,target_gene)
        j=j+1
      }
    }
  }else{
    JustNegativeGene <- c(JustNegativeGene,target_gene)
  }
  regulation <- data.frame(cbind(TFs,type_inter,target_gene_inter,sample_level))
  colnames(regulation) <- c('TF','type','target_gene','sample')
  for(i in 1:nrow(regulation)){
    for(j in i:nrow(regulation)){
      if((regulation[i,'TF']==regulation[j,'TF']) &
         (regulation[i,'sample']==regulation[j,'sample']) &
         (regulation[i,'type']!=regulation[j,'type'])){
        print(paste(i,j,sep = '--'))
      }
    }
  }
  results_TF_up <- list()
  results_RP_down <- list()
  if (length(high_exp_sample) > 0) {
    for(sample_i in 1:length(high_exp_sample)) {
      target_sample<-high_exp_sample[sample_i]
      up_thres <- quantile(peak_luad[, target_sample])[4]
      up_reg_gene <- c()
      down_reg_gene <- c()
      if (peak_luad[candidate_peak_id, target_sample] > up_thres) {
        for (i in 1:length(exsit_gene_list)) {
          t_res_temp <- t.test(gtex_temp[i,], mu = tcga_temp[i, target_sample])
          if (t_res_temp$statistic < 0&t_res_temp$p.value < 0.01) {
            up_reg_gene <- c(up_reg_gene, exsit_gene_list[i])
          }
          if (t_res_temp$statistic > 0&t_res_temp$p.value<0.01) {
            down_reg_gene <- c(down_reg_gene, exsit_gene_list[i])
          }
        }
      }
      TF_up_genes <- up_reg_gene
      results_TF_up[[target_sample]] <- TF_up_genes
      TF_down_genes <- down_reg_gene
      results_RP_down[[target_sample]] <- TF_down_genes
    }
  }
  results_TF_down <- list()
  results_RP_up <- list()
  if (length(down_exp_sample) > 0) {
    for(sample_i in 1:length(down_exp_sample)) {
      target_sample <- down_exp_sample[sample_i]
      up_reg_gene <- c()
      down_reg_gene <- c()
      up_thres <- quantile(peak_luad[,target_sample])[4]
      if(peak_luad[candidate_peak_id,target_sample] > up_thres) {
        for(i in 1:length(exsit_gene_list)) {
          t_res_temp<-t.test(gtex_temp[i,],mu = tcga_temp[i,target_sample])
          if(t_res_temp$statistic < 0&t_res_temp$p.value < 0.01) {
            up_reg_gene<-c(up_reg_gene,exsit_gene_list[i])
          }
          if(t_res_temp$statistic > 0&t_res_temp$p.value < 0.01) {
            down_reg_gene<-c(down_reg_gene,exsit_gene_list[i])
          }
        }
      }
      TF_up_genes <- up_reg_gene
      results_RP_up[[target_sample]] <- TF_up_genes
      TF_down_genes <- down_reg_gene
      results_TF_down[[target_sample]] <- TF_down_genes
    }
  }
  res_positive_TF_up <- formatPositiveResult(results_TF_up)
  res_negative_RP_down <- formatNegativeResult(results_RP_down)
  res_positive_TF_down <- formatPositiveResult(results_TF_down)
  res_negative_RP_up <- formatNegativeResult(results_RP_up)
  freq_threshold <- 0
  high_sample_num<-length(high_exp_sample)
  if (length(res_positive_TF_up) <= 1){
    print('res_negative_RP_down = NULL list')
  }else{
    res_positive_TF_up_frame <- data.frame(table(res_positive_TF_up$TF))
    print(res_positive_TF_up_frame)
  }
  res_positive_TF_up_frame$Freq <- res_positive_TF_up_frame$Freq/high_sample_num
  res_positive_TF_up_frame <- res_positive_TF_up_frame[which(res_positive_TF_up_frame$Freq >= freq_threshold),]
  if (length(res_positive_TF_down) <= 1){
    print('res_negative_RP_down = NULL list')
    res_positive_TF_down_frame <- res_positive_TF_up_frame[1,]
    res_positive_TF_down_frame <- res_positive_TF_down_frame[-1,]
  }else{
    res_positive_TF_down_frame <- data.frame(table(res_positive_TF_down$TF))
    print(res_positive_TF_down_frame)
  }
  res_positive_TF_down_frame$Freq <- res_positive_TF_down_frame$Freq/length(down_exp_sample)
  res_positive_TF_down_frame <- res_positive_TF_down_frame[which(res_positive_TF_down_frame$Freq >= freq_threshold),]
  if (length(res_negative_RP_down) <= 1){
    print('res_negative_RP_down = NULL list')
    res_negative_RP_down_frame <- res_positive_TF_up_frame[1,]
    res_negative_RP_down_frame <- res_negative_RP_down_frame[-1,]
    }else{
    res_negative_RP_down_frame<-data.frame(table(res_negative_RP_down$TF))
    print(res_negative_RP_down_frame)
    }
  res_negative_RP_down_frame$Freq <- res_negative_RP_down_frame$Freq/high_sample_num
  res_negative_RP_down_frame <- res_negative_RP_down_frame[which(res_negative_RP_down_frame$Freq >= freq_threshold),]
  if (length(res_negative_RP_up) <= 1){
    print('res_negative_RP_up = NULL list')
    res_negative_RP_up_frame <- res_positive_TF_up_frame[1, ]
    res_negative_RP_up_frame <- res_negative_RP_up_frame[-1, ]
  }else{
    res_negative_RP_up_frame <- data.frame(table(res_negative_RP_up$TF))
    print(res_negative_RP_up_frame)
  }
  res_negative_RP_up_frame$Freq <- res_negative_RP_up_frame$Freq/length(down_exp_sample)
  res_negative_RP_up_frame <- res_negative_RP_up_frame[which(res_negative_RP_up_frame$Freq >= freq_threshold),]
  temp_positive <- setdiff(res_positive_TF_up_frame$Var1, res_negative_RP_down_frame$Var1)
  temp_negative <- setdiff(res_negative_RP_up_frame$Var1, res_positive_TF_down_frame$Var1)
  real_positive <- setdiff(temp_positive, temp_negative)
  real_negative <- setdiff(temp_negative, temp_positive)
  sel_col <- c('ASPAR','ENCODE','CHEA','MotifMap','TRANSFAC_Curated','TRANSFAC_Predicted')
  if (target_gene%in%negative_gene) {
    res_frame_positive <- data.frame(Var1 = temp_positive)
    res_frame_positive <- FramePositive(res_frame_positive)
    write.table(res_frame_positive,
                file=paste('results/',target_gene,'/', target_gene, '_TF_validation_positive.csv', sep = ''),
                sep=',',
                row.names = F,
                col.names = T)
    if (length(res_frame_positive) == 1) {
      print('no positive result')
    }else{
      hitmaxtirx_positive <- subset(res_frame_positive, select=sel_col)
      hit_tf_positive <- as.vector(res_frame_positive[which(apply(hitmaxtirx_positive, 1, sum)!=0),1])
      if (length(results_TF_up) > 0&length(hit_tf_positive)) {
        hit_res_temp_positive <- SampleHit(results_TF_up, hit_tf_positive)
        write.table(hit_res_temp_positive,
                    file=paste('results/', target_gene,'/', target_gene, '_hit_sample-CNV_positive.csv',sep = ''),
                    sep=',',
                    row.names = F,
                    col.names = T,
                    quote = T)
      }
    }
    res_frame_all <- rbind(res_frame_all,res_frame_positive)
  }else{
    res_frame_negative <- data.frame(Var1 = temp_negative)
    res_frame_negative <- FrameNegative(res_frame_negative)
    write.table(res_frame_negative,
                file=paste('results/',target_gene,'/',target_gene,'_TF_validation_negative.csv',sep = ''),
                sep=',',
                row.names = F,
                col.names = T,
                quote = T)
    if (length(res_frame_negative) == 1) {
      print('no negative result')
    }else{
      hitmaxtirx_negative <- subset(res_frame_negative,select = sel_col)
      hit_tf_negative <- as.vector(res_frame_negative[which(apply(hitmaxtirx_negative,1,sum)!=0),1])
      if (length(results_TF_down) > 0&length(hit_tf_negative)) {
        hit_res_temp_negative <- SampleHit(results_TF_down,hit_tf_negative)
        write.table(hit_res_temp_negative,
                    file=paste('results/',target_gene,'/',target_gene,'_hit_sample-CNV_negative.csv',sep = ''),
                    sep=',',
                    row.names = F,
                    col.names = T,
                    quote = T)
      }
    }
    res_frame_all <- rbind(res_frame_all,res_frame_negative)
  }
  if (target_gene%in%negative_gene) {
    res_all <- rbind(res_positive_TF_down,
                     res_positive_TF_up)
  }else{
    res_all <- rbind(res_negative_RP_up,
                     res_negative_RP_down)
  }
  if (ncol(res_all) == 4) {
    Regulation_result_4 <- rbind(Regulation_result_4,res_all)
  }
  new_TFs_list <- res_all$TF
  new_TFs_list <- as.data.frame(new_TFs_list)
  new_TFs_list[new_TFs_list == 'NULL list'] <- NA
  new_TFs_list <- na.omit(new_TFs_list)
  new_TFs_list <- new_TFs_list[!duplicated(new_TFs_list),]
  new_TFs_list <- as.data.frame(new_TFs_list)
  names(new_TFs_list) <- 'TFs'
  write.table(new_TFs_list,
              paste(target_gene,'_TFs_list.txt',sep = ''),
              quote = F,
              row.names = F,
              col.names = T,
              sep = '\t')
  Regulation_result <- rbind(Regulation_result,regulation)
  if (target_gene%in%negative_gene ) {
    RegulationIndex_res$ImmuneScore <- RegulationIndex_res[,2]*(-1)
  }
  RegulationIndex_results <- rbind.data.frame(immune_score_results,immune_score_res)
  RegulationIndex_TF_up1 <- RegulationIndex_calculation(results_TF_up)
  RegulationIndex_RP_down1 <- RegulationIndex_calculation(results_RP_down)
  RegulationIndex_TF_down1 <- RegulationIndex_calculation(results_TF_down)
  RegulationIndex_up1 <- RegulationIndex_calculation(results_RP_up)
  RegulationIndex1 <- rbind.data.frame(immune_score_TF_up1,immune_score_RP_down1,immune_score_TF_down1,immune_score_RP_up1)
  RegulationIndex <- rbind.data.frame(RegulationIndex,RegulationIndex1)
}
names(RegulationIndex) <- c('Sample','RegulationIndex')
RegulationIndex <- RegulationIndex[order(RegulationIndex$ImmuneScore,decreasing = T),]
write.csv(RegulationIndex,'RegulationIndex.csv',row.names = F)
Regulation_result <- na.omit(Regulation_result)
write.table(Regulation_result,
            'results/Regulation_result.csv',
            sep=',',
            row.names = F,
            col.names = T)


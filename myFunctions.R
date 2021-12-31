
Peak_is_open <- function(candidate_peak_id_input,target_sample_input) {
  up_thres<-quantile(peak_luad[,target_sample_input])[3]
  #print(up_thres)
  #if(TRUE %in% c(peak_luad[candidate_peak_id_input,target_sample_input]>up_thres))
  if (peak_luad[candidate_peak_id_input,target_sample_input]>up_thres) {
    return(TRUE)
  } else {
    return(FALSE)
    }
}
formatPositiveResult <- function(results_summary_input){
  TFs <- c()
  type_inter <- c()
  target_gene_inter <- c()
  sample_level <- c()
  if(length(results_summary_input)) {
    i=1
    while(i<=length(results_summary_input)) {
      j=1
      while(j <= length(results_summary_input[[i]])) {
        TFs<-c(TFs,results_summary_input[[i]][j]);
        type_inter<-c(type_inter,'1');# type 1
        sample_level<-c(sample_level,names(results_summary_input)[i]);
        target_gene_inter<-c(target_gene_inter,target_gene);
        j=j+1
      }
      i=i+1
    }
    format_result <- data.frame(cbind(TFs,type_inter,target_gene_inter,sample_level))
    colnames(format_result) <- c('TF','type','target_gene','sample')
    return(format_result)
  } else {
    print('NULL list')
  }
}

formatNegativeResult <- function(results_summary_input){
  TFs<-c()
  type_inter<-c()
  target_gene_inter<-c()
  sample_level<-c()
  if(length(results_summary_input)){
    i=1
    while(i<=length(results_summary_input)){
      j=1
      while(j <= length(results_summary_input[[i]])){
        TFs<-c(TFs,results_summary_input[[i]][j]);
        type_inter<-c(type_inter,'2');# type 1
        sample_level<-c(sample_level,names(results_summary_input)[i]);
        target_gene_inter<-c(target_gene_inter,target_gene);
        j=j+1
      }
      i=i+1
    }
    format_result<-data.frame(cbind(TFs,type_inter,target_gene_inter,sample_level))
    colnames(format_result)<-c('TF','type','target_gene','sample')
    return(format_result)
  }else{
    print('NULL list')
  }
}

SampleHit <- function(results_summary_input,hit_tf) {
  hit_count<-c()
  sample_level<-c()
  for(i in 1:length(results_summary_input)){#each sample
    j=1
    count_temp<-0
    while (j <= length(results_summary_input[[i]])) {
      if (results_summary_input[[i]][j] %in% hit_tf & names(results_summary_input)[i]%in% colnames(raw_tcga_cnv)) {
        if(raw_tcga_cnv[results_summary_input[[i]][j],names(results_summary_input)[i]]>0)
        {count_temp=count_temp+1}
      }
      j=j+1
    }
    sample_level<-c(sample_level,names(results_summary_input)[i]);
    hit_count<-c(hit_count,count_temp)
  }
  format_result<-data.frame(cbind(sample_level,hit_count))
  colnames(format_result)<-c('sample',paste('hitsby',target_gene,sep = ''))
  return(format_result)
}

SampleHitOnlyCNV <- function(high_exp_sample_input,hit_tf){
  sample_level<-c()
  #tem_ind=0
  for(i in high_exp_sample_input){#each sample
    if(hit_tf %in% row.names(raw_tcga_cnv)& i %in% colnames(raw_tcga_cnv)){
      #tem_ind=tem_ind+1
      #print(tem_ind)
      if(raw_tcga_cnv[hit_tf,i]>0){
        sample_level<-c(sample_level,i)
      }
    }
  }
  return(sample_level)
}

FramePositive <- function(res_frame_input){
  ASPAR_hit <- c()
  encode_hit <- c()
  CHEA_hit <- c()
  MotifMap_hit <- c()
  TRANSFAC_Curated_hit <- c()
  TRANSFAC_Predicted_hit <- c()
  for(i in res_frame_input$Var1){
    ASPAR_res<-ASPAR_Predicted[which(ASPAR_Predicted$source==target_gene&ASPAR_Predicted$target==i),]
    if(nrow(ASPAR_res)){
      ASPAR_hit<-c(ASPAR_hit,1)
    }else{
      ASPAR_hit<-c(ASPAR_hit,0)
    }
    encode_tf_res<-encode_tf[which(encode_tf$source==target_gene&encode_tf$target==i),]
    if(nrow(encode_tf_res)){
      encode_hit<-c(encode_hit,1)
    }else {encode_hit<-c(encode_hit,0)}
    CHEA_tf_res<-CHEA_tf[which(CHEA_tf$source==target_gene&CHEA_tf$target==i),]
    if(nrow(CHEA_tf_res)){
      CHEA_hit<-c(CHEA_hit,1)
    }else {CHEA_hit<-c(CHEA_hit,0)}
    MotifMap_Predicted_res<-MotifMap_Predicted[which(MotifMap_Predicted$source==target_gene&MotifMap_Predicted$target==i),]
    if(nrow(MotifMap_Predicted_res)){
      MotifMap_hit<-c(MotifMap_hit,1)
    }else {MotifMap_hit<-c(MotifMap_hit,0)}
    TRANSFAC_Curated_res<-TRANSFAC_Curated[which(TRANSFAC_Curated$source==target_gene&TRANSFAC_Curated$target==i),]
    if(nrow(TRANSFAC_Curated_res)){
      TRANSFAC_Curated_hit<-c(TRANSFAC_Curated_hit,1)
    }else {TRANSFAC_Curated_hit<-c(TRANSFAC_Curated_hit,0)}
    TRANSFAC_Predicted_res<-TRANSFAC_Predicted[which(TRANSFAC_Predicted$source==target_gene&TRANSFAC_Predicted$target==i),]
    if(nrow(TRANSFAC_Predicted_res)){
      TRANSFAC_Predicted_hit<-c(TRANSFAC_Predicted_hit,1)
    }else {TRANSFAC_Predicted_hit<-c(TRANSFAC_Predicted_hit,0)}
  }
  res_frame_positive$ASPAR <- ASPAR_hit
  res_frame_positive$ENCODE <- encode_hit
  res_frame_positive$CHEA <- CHEA_hit
  res_frame_positive$MotifMap <- MotifMap_hit
  res_frame_positive$TRANSFAC_Curated <- TRANSFAC_Curated_hit
  res_frame_positive$TRANSFAC_Predicted <- TRANSFAC_Predicted_hit
  return(res_frame_positive)
}

FrameNegative <- function(res_frame_input){
  ASPAR_hit <- c()
  encode_hit <- c()
  CHEA_hit <- c()
  MotifMap_hit <- c()
  TRANSFAC_Curated_hit <- c()
  TRANSFAC_Predicted_hit <- c()
  for(i in res_frame_input$Var1){
    ASPAR_res<-ASPAR_Predicted[which(ASPAR_Predicted$source==target_gene&ASPAR_Predicted$target==i),]
    if(nrow(ASPAR_res)){
      ASPAR_hit<-c(ASPAR_hit,1)
    }else{
      ASPAR_hit<-c(ASPAR_hit,0)
      }
    encode_tf_res<-encode_tf[which(encode_tf$source==target_gene&encode_tf$target==i),]
    if(nrow(encode_tf_res)){
      encode_hit<-c(encode_hit,1)
    }else {encode_hit<-c(encode_hit,0)}
    CHEA_tf_res<-CHEA_tf[which(CHEA_tf$source==target_gene&CHEA_tf$target==i),]
    if(nrow(CHEA_tf_res)){
      CHEA_hit<-c(CHEA_hit,1)
    }else {CHEA_hit<-c(CHEA_hit,0)}
    MotifMap_Predicted_res<-MotifMap_Predicted[which(MotifMap_Predicted$source==target_gene&MotifMap_Predicted$target==i),]
    if(nrow(MotifMap_Predicted_res)){
      MotifMap_hit<-c(MotifMap_hit,1)
    }else {MotifMap_hit<-c(MotifMap_hit,0)}
    TRANSFAC_Curated_res<-TRANSFAC_Curated[which(TRANSFAC_Curated$source==target_gene&TRANSFAC_Curated$target==i),]
    if(nrow(TRANSFAC_Curated_res)){
      TRANSFAC_Curated_hit<-c(TRANSFAC_Curated_hit,1)
    }else {TRANSFAC_Curated_hit<-c(TRANSFAC_Curated_hit,0)}
    TRANSFAC_Predicted_res<-TRANSFAC_Predicted[which(TRANSFAC_Predicted$source==target_gene&TRANSFAC_Predicted$target==i),]
    if(nrow(TRANSFAC_Predicted_res)){
      TRANSFAC_Predicted_hit<-c(TRANSFAC_Predicted_hit,1)
    }else {TRANSFAC_Predicted_hit<-c(TRANSFAC_Predicted_hit,0)}
  }
  res_frame_negative$ASPAR <- ASPAR_hit
  res_frame_negative$ENCODE <- encode_hit
  res_frame_negative$CHEA <- CHEA_hit
  res_frame_negative$MotifMap <- MotifMap_hit
  res_frame_negative$TRANSFAC_Curated <- TRANSFAC_Curated_hit
  res_frame_negative$TRANSFAC_Predicted <- TRANSFAC_Predicted_hit
  return(res_frame_negative)
}

RegulationIndex_calculation <- function(results_summary_input){
  immune_score_all <- c()
  if(length(results_summary_input) > 1) {
    q=1
    while (q <= length(results_summary_input)) {
      ImmuneScore <- length(results_summary_input[[q]])/nrow(gene_list)/10
      SampleSelected <- names(results_summary_input)[q]
      immune_score_sample <- cbind.data.frame(SampleSelected,ImmuneScore)
      immune_score_all <- rbind.data.frame(immune_score_all,immune_score_sample)
      q=q+1
    }
  }
  immune_score_all <- as.data.frame(immune_score_all)
  return(immune_score_all)
}


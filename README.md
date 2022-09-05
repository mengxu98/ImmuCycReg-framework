# ImmuCycReg-framework

These code is for paper: Integration of single sample and population analysis for understanding Immune evasion mechanisms of lung cancer

The core code is the file: ImmuCycReg framework.R and L0Reg framework.R

The requirement of input data:

  The TCGA RNA-seq is necessary
  
  The ATAC-seq is necessary
  
  The CNV data is not required
  
## DATA links:

  Processed data from this study are available in the reproducibility GitHub repository (https://github.com/mengxu98/ImmuCycReg-framework/tree/main/data).
  
  The Cancer Genome Atlas (TCGA) and the Genotype-Tissue Expression (GTEx) datasets were downloaded from figshare.
  
     https://figshare.com/articles/dataset/Data_record_1/5330539
    
     https://figshare.com/articles/dataset/Data_record_2/5330575
     
     https://figshare.com/articles/dataset/Data_record_3/5330593
     
  Transposase-Accessible Chromatin with high throughput sequencing (ATAC-seq) was downloaded from UCSC-Xena (https://atacseq.xenahubs.net). 
  
  Copy number variations (CNV) dataset was downloaded from GISTIC2.0 (https://api.gdc.cancer.gov/data/7d64377f-2cea-4ee3-917f-8fcfbcd999e7).
  
  Genome annotation file was downloaded from hg38.ensGene.gtf (ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/).

  
## The environment required:
  
  R 4.1.2
  
  NMF==0.24.0
  
  DESeq2==1.32.0
  
  L0Learn==2.0.3
  
  e1071==1.7-11
  
  glmnet==4.1-4
  
  timeROC==0.4
  
  rms==6.3-0
  
  survival==3.3-1
  
  ggplot2==3.3.6

  bedtools v2.27.167
  
  PROMO: http://alggen.lsi.upc.es/cgi-bin/promo_v3/promo/promoinit.cgi?dirDB=TF_8.3

  Cytoscape 3.8.2
  
  ClueGo v2.5.9

  JAVA v18.0.1.1
  
  GeneNetWeaver: http://gnw.sourceforge.net/
  

# &#x1F4D7;ImmuCycReg-framework
<img src="https://img.shields.io/badge/-R-green"/> <img src="https://img.shields.io/badge/-Immune%20escape%20analysis-blue"/> <img src="https://img.shields.io/badge/-Gene%20Regulatory%20Network-blue"/> <img src="https://img.shields.io/eclipse-marketplace/last-update/mengxu98?style=flat-square"/><br/>
The code repository is for paper: Integration of single sample and population analysis for understanding Immune evasion mechanisms of lung cancer (In revision).<br/>
<img src="https://github.com/mengxu98/ImmuCycReg-framework/blob/main/Workflow.png"/><br/>
## &#x1F537;DATA links:
  &#x1F538;Processed data from this study are available in the reproducibility GitHub repository:<br/>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```html https://github.com/mengxu98/ImmuCycReg-framework/tree/main/data```<br/>
  &#x1F538;The Cancer Genome Atlas (TCGA) and the Genotype-Tissue Expression (GTEx) datasets were downloaded from figshare:<br/>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`https://figshare.com/articles/dataset/Data_record_1/5330539`<br/>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`https://figshare.com/articles/dataset/Data_record_2/5330575`<br/>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`https://figshare.com/articles/dataset/Data_record_3/5330593`<br/>
  &#x1F538;Transposase-Accessible Chromatin with high throughput sequencing (ATAC-seq) was downloaded from UCSC-Xena:<br/>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`https://atacseq.xenahubs.net`<br/>
  &#x1F538;Copy number variations (CNV) dataset was downloaded from GISTIC2.0:<br/>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`https://api.gdc.cancer.gov/data/7d64377f-2cea-4ee3-917f-8fcfbcd999e7`<br/>
  &#x1F538;Genome annotation file was downloaded from hg38.ensGene.gtf:<br/>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/`<br/>
## &#x1F537;The input data:
  &#x1F538;The TCGA RNA-seq is necessary<br/>
  &#x1F538;The ATAC-seq is necessary<br/>
  &#x1F538;The CNV data is not required<br/>
## &#x1F537;The environment and softwares required:
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
<!--
## If you encounter any problems when use these code, please contact me by Wechat or QQ: 
Wechat: <img src="https://github.com/mengxu98/scGRN-L0/blob/master/contact/Wechat.jpg" width="100" height="100" alt="Wechat"/> QQ: <img src="https://github.com/mengxu98/scGRN-L0/blob/master/contact/QQ.PNG" width="100" height="100" alt="QQ"/><br/>
-->

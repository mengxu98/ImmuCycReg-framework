

library(clusterProfiler)
library(DOSE)
library(dplyr)
library(gage)
library(gageData)
library(org.Hs.eg.db)
library(pathview)
library(stringr)
workplace <- paste(getwd(),'/results',sep = '')
setwd(workplace)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs =  kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs,3)

for (j in 1:4) {
  if (dir.exists(sprintf('workplace/cluster%s',j)) == F) {
    dir.create(file.path(workplace, paste('cluster', j, sep = '')))
  }
  dir <- paste(workplace,'/cluster',j,sep = '')
  setwd(dir)
  sig.gene <- read.table("DESeq2_res.txt",
                         sep="\t",
                         head=T,
                         comment.char = "")
  sig.gene$otu_id[which(sig.gene$otu_id=='C10orf54')]='VSIR'
  gene.df <- bitr(sig.gene$otu_id, fromType = "SYMBOL",
                  toType = c(
                    "ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  boole_value <- row.names(sig.gene) %in% row.names(gene.df)
  sig.gene_filter <- sig.gene[which(boole_value),]
  ENTREZID <- gene.df$ENTREZID
  foldchanges <- sig.gene_filter$log2FoldChange
  names(foldchanges) = ENTREZID
  keggres = gage(foldchanges, gsets = kegg.sets.hs, same.dir = TRUE)
  lapply(keggres, head)
  keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
    tbl_df() %>% 
    filter(row_number()<=10) %>% 
    .$id %>% 
    as.character()
  write.table(keggrespathways,'keggrespathways.txt',
              row.names = F,
              quote = F)
  keggresids = substr(keggrespathways, start=1, stop=8)
  keggresids
  plot_pathway = function(pid) pathview(gene.data=foldchanges,
                                        pathway.id=pid,
                                        species="hsa",
                                        new.signature=FALSE)
  tmp1 = sapply(keggresids,
               function(pid) pathview(gene.data=foldchanges,
                                      pathway.id=pid,
                                      species="hsa"))
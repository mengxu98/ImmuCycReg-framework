

library(openxlsx)
des <- read.csv("DESeq2_log2FC_res.csv")
imcy <- read.table("../data/immunecycle.txt", sep = "\t", header = T)
des <- des[, -1]
des <- des[, -6]
des <- na.omit(des)
des$reg <- ""
genes <- des$gene
for (i in 1:length(genes)) {
  gene <- genes[i]
  if (gene %in% imcy$GeneSymbol) {
    attr <- imcy[which(imcy$GeneSymbol == gene), ]
    if (attr$Direction == "positive") {
      des2 <- des[which(des$gene == gene), ]

      if (des2$cluster1 >= 1 &
        des2$cluster2 < 1 &
        des2$cluster3 < 1 &
        des2$cluster4 < 1) {
        des$reg[i] <- "cluster2,3,4 Down"
      }
      if (des2$cluster1 < 1 &
        des2$cluster2 >= 1 &
        des2$cluster3 < 1 &
        des2$cluster4 < 1) {
        des$reg[i] <- "cluster1,3,4 Down"
      }
      if (des2$cluster1 < 1 &
        des2$cluster2 < 1 &
        des2$cluster3 >= 1 &
        des2$cluster4 < 1) {
        des$reg[i] <- "cluster1,2,4 Down"
      }
      if (des2$cluster1 < 1 &
        des2$cluster2 < 1 &
        des2$cluster3 < 1 &
        des2$cluster4 >= 1) {
        des$reg[i] <- "cluster1,2,4 Down"
      }
      if (des2$cluster1 >= 1 &
        des2$cluster2 >= 1 &
        des2$cluster3 < 1 &
        des2$cluster4 < 1) {
        des$reg[i] <- "cluster3,4 Down"
      }
      if (des2$cluster1 >= 1 &
        des2$cluster2 >= 1 &
        des2$cluster3 >= 1 &
        des2$cluster4 < 1) {
        des$reg[i] <- "cluster4 Down"
      }
      if (des2$cluster1 >= 1 &
        des2$cluster2 >= 1 &
        des2$cluster3 < 1 &
        des2$cluster4 >= 1) {
        des$reg[i] <- "cluster3 Down"
      }
      if (des2$cluster1 >= 1 &
        des2$cluster2 < 1 &
        des2$cluster3 >= 1 &
        des2$cluster4 < 1) {
        des$reg[i] <- "cluster2,4 Down"
      }
      if (des2$cluster1 >= 1 &
        des2$cluster2 < 1 &
        des2$cluster3 >= 1 &
        des2$cluster4 >= 1) {
        des$reg[i] <- "cluster2 Down"
      }
      if (des2$cluster1 >= 1 &
        des2$cluster2 < 1 &
        des2$cluster3 < 1 &
        des2$cluster4 >= 1) {
        des$reg[i] <- "cluster2,3 Down"
      }
      if (des2$cluster1 < 1 &
        des2$cluster2 >= 1 &
        des2$cluster3 >= 1 &
        des2$cluster4 < 1) {
        des$reg[i] <- "cluster1,4 Down"
      }
      if (des2$cluster1 < 1 &
        des2$cluster2 >= 1 &
        des2$cluster3 >= 1 &
        des2$cluster4 >= 1) {
        des$reg[i] <- "cluster1 Down"
      }
      if (des2$cluster1 < 1 &
        des2$cluster2 >= 1 &
        des2$cluster3 < 1 &
        des2$cluster4 >= 1) {
        des$reg[i] <- "cluster1,3 Down"
      }
      if (des2$cluster1 < 1 &
        des2$cluster2 < 1 &
        des2$cluster3 >= 1 &
        des2$cluster4 >= 1) {
        des$reg[i] <- "cluster1,2 Down"
      }
      #
    } else {
      des2 <- des[which(des$gene == gene), ]
      #
      if (des2$cluster1 >= 1 &
        des2$cluster2 < 1 &
        des2$cluster3 < 1 &
        des2$cluster4 < 1) {
        des$reg[i] <- "cluster1 UP"
      }
      if (des2$cluster1 < 1 &
        des2$cluster2 >= 1 &
        des2$cluster3 < 1 &
        des2$cluster4 < 1) {
        des$reg[i] <- "cluster2 UP"
      }
      if (des2$cluster1 < 1 &
        des2$cluster2 < 1 &
        des2$cluster3 >= 1 &
        des2$cluster4 < 1) {
        des$reg[i] <- "cluster3 UP"
      }
      if (des2$cluster1 < 1 &
        des2$cluster2 < 1 &
        des2$cluster3 < 1 &
        des2$cluster4 >= 1) {
        des$reg[i] <- "cluster4 UP"
      }
      if (des2$cluster1 >= 1 &
        des2$cluster2 >= 1 &
        des2$cluster3 < 1 &
        des2$cluster4 < 1) {
        des$reg[i] <- "cluster1,2 UP"
      }
      if (des2$cluster1 >= 1 &
        des2$cluster2 >= 1 &
        des2$cluster3 >= 1 &
        des2$cluster4 < 1) {
        des$reg[i] <- "cluster1,2,3 UP"
      }
      if (des2$cluster1 >= 1 &
        des2$cluster2 >= 1 &
        des2$cluster3 < 1 &
        des2$cluster4 >= 1) {
        des$reg[i] <- "cluster1,2,4 UP"
      }
      if (des2$cluster1 >= 1 &
        des2$cluster2 < 1 &
        des2$cluster3 >= 1 &
        des2$cluster4 < 1) {
        des$reg[i] <- "cluster1,3 UP"
      }
      if (des2$cluster1 >= 1 &
        des2$cluster2 < 1 &
        des2$cluster3 >= 1 &
        des2$cluster4 >= 1) {
        des$reg[i] <- "cluster1,3,4 UP"
      }
      if (des2$cluster1 >= 1 &
        des2$cluster2 < 1 &
        des2$cluster3 < 1 &
        des2$cluster4 >= 1) {
        des$reg[i] <- "cluster1,4 UP"
      }
      if (des2$cluster1 < 1 &
        des2$cluster2 >= 1 &
        des2$cluster3 >= 1 &
        des2$cluster4 < 1) {
        des$reg[i] <- "cluster2,3 UP"
      }
      if (des2$cluster1 < 1 &
        des2$cluster2 >= 1 &
        des2$cluster3 >= 1 &
        des2$cluster4 >= 1) {
        des$reg[i] <- "cluster2,3,4 UP"
      }
      if (des2$cluster1 < 1 &
        des2$cluster2 >= 1 &
        des2$cluster3 < 1 &
        des2$cluster4 >= 1) {
        des$reg[i] <- "cluster2,4 UP"
      }
      if (des2$cluster1 < 1 &
        des2$cluster2 < 1 &
        des2$cluster3 >= 1 &
        des2$cluster4 >= 1) {
        des$reg[i] <- "cluster3,4 UP"
      }
    }
  }
}
sheets2 <- list(des)
write.xlsx(sheets2, "DESeq2_log2FC_res.xlsx")

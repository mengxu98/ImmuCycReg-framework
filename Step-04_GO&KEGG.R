

library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("GOplot")
library("ggrepel")
library(ggnewscale)

for (j in 1:4) {
       if (dir.exists(paste("D:/test/immune landscape_LUAD/immune cycle/results/cluster", j, sep = "")) == F) {
              dir.create(file.path("D:/test/immune landscape_LUAD/immune cycle/results", paste("cluster", j, sep = "")))
       }
       dir <- paste("D:/test/immune landscape_LUAD/immune cycle/results/cluster", j, sep = "")
       # print(dir)
       setwd(dir)

       if (dir.exists(file.path(getwd(), paste("GO", sep = ""))) == F) {
              dir.create(file.path(getwd(), paste("GO", sep = "")))
       }
       if (dir.exists(file.path(getwd(), paste("KEGG", sep = ""))) == F) {
              dir.create(file.path(getwd(), paste("KEGG", sep = "")))
       }

       sig_gene <- read.table("DESeq2_res_sort.csv",
              header = T,
              # row.names = 1,
              sep = ","
       )
       sig_gene <- na.omit(sig_gene)
       sig_gene$otu_id[which(sig_gene$otu_id == "C10orf54")] <- "VSIR"
       row.names(sig_gene) <- sig_gene$otu_id
       genes <- as.vector(row.names(sig_gene))
       logFC <- sig_gene$log2FoldChange

       entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG,
              ifnotfound = NA
       )

       entrezIDs <- as.character(entrezIDs)

       symbol2id <- cbind(
              gene = genes,
              logFC,
              entrezID = entrezIDs
       )

       write.table(symbol2id,
              file = "symbol2id.txt",
              sep = "\t",
              quote = F,
              row.names = F
       )
       #### read data####
       pathway_data <- read.table("symbol2id.txt", sep = "\t", header = T, check.names = F)
       pathway_data <- pathway_data[is.na(pathway_data[, "entrezID"]) == F, ]
       pathway_gene <- pathway_data$entrezID

       #### GO####
       GO_kk <- enrichGO(
              gene = pathway_gene,
              OrgDb = org.Hs.eg.db,
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05,
              ont = "all",
              readable = T
       )

       write.table(GO_kk,
              file = "GO/GO_kk.txt",
              sep = "\t",
              quote = F,
              row.names = F
       )


       # png(file="GO/GO_barplot.png",width = 800,height = 600)
       # p1 <- barplot(GO_kk,
       #               drop = TRUE,
       #               showCategory =10,
       #               split="ONTOLOGY")+
       #   facet_grid(ONTOLOGY~., scale='free')
       # print(p1)
       # dev.off()

       p1 <- barplot(GO_kk,
              drop = TRUE,
              showCategory = 10,
              split = "ONTOLOGY"
       ) +
              facet_grid(ONTOLOGY ~ ., scale = "free")
       ggsave("GO/GO_barplot.png",
              p1,
              width = 12,
              height = 8
       )

       # png(file="GO/GO_bubble.png",width = 800,height = 600)
       # dotplot(GO_kk,showCategory = 10,split="ONTOLOGY",orderBy="x") + facet_grid(ONTOLOGY~., scale='free')
       # dev.off()
       p2 <- dotplot(GO_kk,
              showCategory = 10,
              split = "ONTOLOGY",
              orderBy = "x"
       ) +
              facet_grid(ONTOLOGY ~ .,
                     scale = "free"
              )
       ggsave("GO/GO_bubble.png",
              p2,
              width = 12,
              height = 8
       )
       ####


       GO_ego <- read.table("GO/GO_kk.txt", header = T, sep = "\t", check.names = F)

       go <- data.frame(
              Category = "All",
              ID = GO_ego$ID,
              Term = GO_ego$Description,
              Genes = gsub(
                     "/", ", ",
                     GO_ego$geneID
              ),
              adj_pval = GO_ego$p.adjust
       )
       ####
       symbol2id_go <- read.table("symbol2id.txt", sep = "\t", header = T, check.names = F)
       symbol2id_go <- symbol2id_go[is.na(symbol2id_go[, "entrezID"]) == F, ]
       ####

       GO_circ_data <- data.frame(ID = symbol2id_go$gene, logFC = symbol2id_go$logFC)

       row.names(GO_circ_data) <- GO_circ_data[, 1]

       GO_circ <- circle_dat(go, GO_circ_data)

       termNum <- 5
       geneNum <- nrow(GO_circ_data)

       chord <- chord_dat(GO_circ, GO_circ_data[1:geneNum, ], go$Term[1:termNum])
       # png(file="GO/GO_circ.png",width = 800,height = 600)
       # GOChord(chord,
       #         space = 0.001,
       #         gene.order = 'logFC',
       #         gene.space = 0.25,
       #         gene.size = 5,
       #         border.size = 0.1,
       #         process.label = 6.5)
       # dev.off()
       p3 <- GOChord(chord,
              space = 0.001,
              gene.order = "logFC",
              gene.space = 0.25,
              gene.size = 5,
              border.size = 0.1,
              process.label = 6.5
       )
       ggsave("GO/GO_circ.png",
              p3,
              width = 12,
              height = 8
       )

       termCol <- c(
              "#223D6C",
              "#D20A13",
              "#FFD121",
              "#088247",
              "#58CDD9",
              "#7A142C",
              "#5D90BA",
              "#431A3D",
              "#91612D",
              "#6E568C",
              "#E0367A",
              "#D8D155",
              "#64495D",
              "#7CC767"
       )
       # png(file="GO/GO_cluster.png",width = 1000,height = 800)
       # GOCluster(circ.gsym,
       #           go$Term[1:termNum],
       #           lfc.space = 0.2,
       #           lfc.width = 1,
       #           term.col = termCol[1:termNum],
       #           term.space = 0.2,
       #           term.width = 1)
       # dev.off()
       p4 <- GOCluster(circ.gsym,
              go$Term[1:termNum],
              lfc.space = 0.2,
              lfc.width = 1,
              term.col = termCol[1:termNum],
              term.space = 0.2,
              term.width = 1
       )
       ggsave("GO/GO_cluster.png",
              p4,
              width = 12,
              height = 8
       )
       #### KEGG####
       kegg_kk <- enrichKEGG(
              gene = pathway_gene,
              organism = "hsa",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05
       )

       write.table(kegg_kk, file = "KEGG/kegg_kk.txt", sep = "\t", quote = F, row.names = F)

       p5 <- barplot(kegg_kk,
              drop = TRUE,
              showCategory = 30
       )
       ggsave("KEGG/kegg_barplot.png",
              p5,
              width = 12,
              height = 8
       )

       p6 <- dotplot(kegg_kk,
              showCategory = 30,
              orderBy = "x"
       )
       ggsave("KEGG/kegg_bubble.png",
              p6,
              width = 12,
              height = 8
       )
}

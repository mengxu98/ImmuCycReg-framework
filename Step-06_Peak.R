

rm(list = ls())

library("psych")
library("ggpubr")
library("ggthemes")
library("rtracklayer")
library("corrgram")
library("GGally")
library("ggplot2")
library("patchwork")

pathRead <- "data/"
pathSave <- "../Results/"

# 1:select peak
# 2:make bed file

# Load data ---------------------------------------------------------------
load(paste0(pathSave, "TCGA-LUAD.Rdata"))
load(paste0(pathSave, "Peak-LUAD.Rdata"))
load(paste0(pathSave, "genes_adj_peak.Rdata"))
row_names <- row.names(genes_adj_peak)
genes_temp <- as.vector(genes_adj_peak[, 1])
genes_list <- strsplit(genes_temp, ",")
load("../data ATAC-seq/allgene_tcga_mrna.Rdata")
sample_name <- colnames(allgene_tcga_mrna)
sample_name <- intersect(colnames(tcga_luad), colnames(peak_luad))
# load(paste0(pathSave, "Geneinfo_df.Rdata"))
target_genes_list <- read.table(paste0(pathRead, "Genes_17.txt"), header = T)
target_genes_list <- target_genes_list$gene

# Run ---------------------------------------------------------------------
peak_seq_pos_all <- c()
volcano_plot_pvalue_list <- list()
for (j in 1:length(target_genes_list)) {
  if (!dir.exists(paste0(pathSave, "BED"))) {
    dir.create(paste0(pathSave, "BED"))
  }
  target_gene <- target_genes_list[j]
  peak_around_gene <- c()
  for (i in 1:nrow(genes_adj_peak)) {
    bool_exsit <- target_gene %in% genes_list[[i]]
    if (bool_exsit) {
      peak_around_gene <- c(peak_around_gene, row_names[i])
    }
  }

  # choice 1: correlation ---------------------------------------------------
  if (!dir.exists(paste0("../results ATAC-seq/corrgram/", target_gene))) {
    dir.create(paste0("../results ATAC-seq/corrgram/", target_gene), recursive = TRUE)
  }

  # select peaks according to distance and score correlation (population level)
  Y <- t(allgene_tcga_mrna[target_gene, ]) %>% as.data.frame()
  X <- t(peak_luad[peak_around_gene, sample_name]) %>% as.data.frame()
  peak_name <- colnames(X)

  for (i in 1:length(peak_name)) {
    peak <- X[, i]
    peak_data <- cbind.data.frame(Y, peak)
    names(peak_data) <- c(target_gene, peak_name[i])
    p <- ggpairs(peak_data,
      title = target_gene,
      upper = list(continuous = wrap("cor",
        method = "spearman"
      ))
    )
    print(p)
    ggsave(paste("../results ATAC-seq/corrgram/", target_gene, "/", target_gene, "_", peak_name[i], "_corrgram.png", sep = ""),
      width = 4,
      height = 3
    )
  }
  peak_corr_r <- c()
  peak_corr_p <- c()
  for (i in 1:length(peak_around_gene)) {
    peak_corr <- corr.test(X[, i],
      Y,
      method = "spearman",
      adjust = "none"
    )
    peak_corr_r <- c(peak_corr_r, peak_corr$r)
    peak_corr_p <- c(peak_corr_p, peak_corr$p)
  }

  deg.data <- data.frame(
    Symbol = peak_around_gene,
    corr = peak_corr_r,
    pval = peak_corr_p
  )
  deg.data$logP <- -log10(deg.data$pval)
  deg.data$Group <- "not-significant"
  deg.data$Group[which((deg.data$pval < 0.05) & (deg.data$corr > 0.1))] <- "up-regulated"
  deg.data$Group[which((deg.data$pval < 0.05) & (deg.data$corr < -0.1))] <- "down-regulated"

  deg.data$Label <- ""
  deg.data <- deg.data[order(deg.data$pval), ]

  up.peak <- head(deg.data$Symbol[which(deg.data$Group == "up-regulated")], 20)
  down.peak <- head(deg.data$Symbol[which(deg.data$Group == "down-regulated")], 20)
  not.sig.peak <- head(deg.data$Symbol[which(deg.data$Group == "not-significant")], 20)
  deg.top.peak <- c(as.character(up.peak), as.character(down.peak), as.character(not.sig.peak))
  deg.data$Label[match(deg.top.peak, deg.data$Symbol)] <- deg.top.peak

  write.csv(deg.data,
    file = paste("../results ATAC-seq/corrgram/", target_gene, "_peak_corr.csv", sep = ""),
    row.names = F
  )

  high_cor_peak <- as.character(up.peak[1])
  down_cor_peak <- as.character(down.peak[1])

  # choice2: distance----------------------------------------------------------
  # http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/

  pos_col_name <- c("chrom", "chromStart", "chromEnd", "strand")
  # temp_res <- geneinfo_df[which(geneinfo_df[, "gene_name"] == target_gene), c("seqnames", "type", "gene_id", "start", "end", "width")]
  temp_res <- genes_adj_peak[which(genes_adj_peak[, "gene_name"] == target_gene), c("seqnames", "type", "gene_id", "start", "end", "width")]
  temp_res <- temp_res[which(temp_res[, "type"] == "transcript"), ]
  # select the NEAREST long region
  order_index <- order(temp_res[, "start"], decreasing = FALSE)
  most_start <- temp_res[order_index[1], "start"]
  most_end <- temp_res[order_index[1], "end"]
  most_mid <- (most_start + most_end) / 2

  peaks_pos_temp <- genes_adj_peak[peak_around_gene, pos_col_name]
  peaks_pos_temp$mid <- (peaks_pos_temp$chromStart + peaks_pos_temp$chromEnd) / 2
  peaks_pos_temp$mid_distance <- peaks_pos_temp$mid - most_mid
  peaks_pos_temp$start_distance <- peaks_pos_temp$chromEnd - most_start
  distance_thres <- -8000
  nearest_peaks <- row.names(peaks_pos_temp[which(peaks_pos_temp$start_distance < 0 &
    peaks_pos_temp$start_distance > distance_thres &
    peaks_pos_temp$chromEnd < most_start), ])
  no_peaks <- length(nearest_peaks)
  nearest_peak <- nearest_peaks[no_peaks] # the most nearest peak

  # choice 3:get top peaks in score BUT SCORE IS SAMPLE LEVEL
  top_k <- min(1, length(peak_around_gene))
  target_sample <- "TCGA-73-A9RS-01"

  order_peak_id <- peak_around_gene[order(peak_luad[peak_around_gene, target_sample], decreasing = TRUE)]
  highest_score_peak <- order_peak_id[1:top_k]

  # get the sequence position of peak
  if (length(nearest_peak) != 0) {
    if (nearest_peak %in% c(high_cor_peak, down_cor_peak)) {
      candidate_peak_id <- nearest_peak
    }
  } else {
    if (highest_score_peak %in% c(high_cor_peak, down_cor_peak)) {
      candidate_peak_id <- highest_score_peak
    } else {
      if (!is.na(high_cor_peak)) {
        candidate_peak_id <- c(high_cor_peak)
      } else {
        candidate_peak_id <- c(down_cor_peak)
      }
    }
  }
  candidate_peak_id <- na.omit(candidate_peak_id)

  peak_seq_pos <- genes_adj_peak[candidate_peak_id, pos_col_name]
  peak_seq_pos[, "chromStart"] <- peak_seq_pos[, "chromStart"] - 1

  if (dir.exists("../results ATAC-seq/volcano") == FALSE) {
    dir.create("../results ATAC-seq/volcano")
  }

  if (nrow(peak_seq_pos) > 0) {
    volcano_plot_pvalue <- ggscatter(deg.data,
      x = "corr",
      y = "logP",
      color = "Group",
      palette = c("#006699", "#008B00", "#EE3A8C"),
      label = deg.data$Label,
      font.label = 8,
      repel = T,
      xlab = paste(target_gene, "correlation"),
      ylab = "-log10(P-val)",
      size = 3
    ) +
      theme_bw() +
      theme(
        text = element_text(size = 16),
        legend.position = "none"
      ) +
      geom_hline(yintercept = 1.30, linetype = "dashed") +
      geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed")
    volcano_plot_pvalue
    volcano_plot_pvalue_list[[j]] <- volcano_plot_pvalue
    ggsave(paste("../results ATAC-seq/volcano/", target_gene, "_volcano.pdf", sep = ""),
      volcano_plot_pvalue,
      width = 6.5,
      height = 5
    )

    write.table(peak_seq_pos[, 1:3],
      file = paste("../results ATAC-seq/BED/", target_gene, "_around_peak.bed", sep = ""),
      sep = "\t",
      row.names = F,
      col.names = F,
      quote = F
    )
    peak_seq_pos$gene <- target_gene
    peak_seq_pos_all <- rbind(peak_seq_pos_all, peak_seq_pos)
  }
}

volcano_plot_pvalue_list[[1]] +
  volcano_plot_pvalue_list[[2]] +
  volcano_plot_pvalue_list[[3]] +
  volcano_plot_pvalue_list[[4]] +
  volcano_plot_pvalue_list[[5]] +
  volcano_plot_pvalue_list[[6]] +
  volcano_plot_pvalue_list[[7]] +
  volcano_plot_pvalue_list[[8]] +
  plot_layout(ncol = 4)
ggsave(paste0("../results ATAC-seq/results/boxplot/", "peaks1.png"), width = 12, height = 6, dpi = 600)
ggsave(paste0("../results ATAC-seq/results/boxplot/", "peaks1.pdf"), width = 12, height = 6)

volcano_plot_pvalue_list[[9]] +
  volcano_plot_pvalue_list[[10]] +
  volcano_plot_pvalue_list[[11]] +
  volcano_plot_pvalue_list[[12]] +
  volcano_plot_pvalue_list[[13]] +
  volcano_plot_pvalue_list[[14]] +
  volcano_plot_pvalue_list[[15]] +
  volcano_plot_pvalue_list[[16]] +
  volcano_plot_pvalue_list[[17]] +
  plot_layout(ncol = 4)

ggsave(paste0("../results ATAC-seq/results/boxplot/", "peaks2.png"), width = 12, height = 9, dpi = 600)
ggsave(paste0("../results ATAC-seq/results/boxplot/", "peaks2.pdf"), width = 12, height = 9)

peak_seq_pos_all$peak <- row.names(peak_seq_pos_all)
write.table(peak_seq_pos_all[, 1:3],
  "../results ATAC-seq/allBED.bed",
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)
write.table(peak_seq_pos_all[c(5:6, 1:4)], "../results ATAC-seq/all_peaks.csv",
  sep = ",",
  quote = F,
  row.names = F
)

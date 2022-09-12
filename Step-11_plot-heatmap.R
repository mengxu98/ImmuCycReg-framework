
# select_genes <- up_reg_gene
select_genes <- target_gene_list$gene

intersect_samples <- intersect(colnames(raw_tcga_cnv), common_samples)
cnv_data <- raw_tcga_cnv[select_genes, high_exp_sample]
expmat <- raw_tcga[select_genes, high_exp_sample]

mat <- as.matrix(cnv_data)
mat[which(mat == -2)] <- "AMP"
mat[which(mat == -1)] <- "AMP"
mat[which(mat == 0)] <- ""
mat[which(mat == 1)] <- "AMP1"
mat[which(mat == 2)] <- "MUT1"

Heatmap(log10(expmat + 1),
  name = "expr",
  show_column_names = F,
  col = circlize::colorRamp2(c(0, 3, 5), c("#2c58af", "white", "#d43a35")),
  width = unit(8, "cm"), cluster_rows = T, cluster_columns = T, show_row_dend = F, show_column_dend = T
) + oncoPrint(
  mat,
  alter_fun = alter_fun,
  col = col,
  column_title = column_title,
  heatmap_legend_param = heatmap_legend_param,
  show_column_names = F,
  width = unit(8, "cm")
)

cnv_data <- raw_tcga_cnv[select_genes, intersect_samples]
expmat <- raw_tcga[select_genes, intersect_samples]
Heatmap(cnv_data)
colors <- structure(circlize::rand_color(5), names = c("-2", "-1", "0", "1", "2"))

colors <- structure(1:5, names = c("-2", "-1", "0", "1", "2"))
colors <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Heatmap(cnv_data, col = colors)



mat <- as.matrix(cnv_data)
mat[which(mat == -2)] <- "HOMDEL"
mat[which(mat == -1)] <- "AMP"
mat[which(mat == 0)] <- ""
mat[which(mat == 1)] <- "AMP1"
mat[which(mat == 2)] <- "MUT1"



col <- c("HOMDEL" = "#330099", "AMP" = "#0099cc", "AMP1" = "#ff6600", "MUT1" = "#ff0066")
alter_fun <- list(
  background = alter_graphic("rect", fill = "#CCCCCC"),
  HOMDEL = alter_graphic("rect", fill = col["HOMDEL"]),
  AMP = alter_graphic("rect", fill = col["AMP"]),
  AMP1 = alter_graphic("rect", fill = col["AMP1"]),
  MUT1 = alter_graphic("rect", fill = col["MUT1"])
)
column_title <- "OncoPrint for TCGA Lung Adenocarcinoma, genes in Ras Raf MEK JNK signalling"
heatmap_legend_param <-
  list(
    title = "Alternations",
    at = c("HOMDEL", "AMP", "AMP1", "MUT1"),
    labels = c("-2", "-1", "1", "2")
  )

oncoPrint(
  mat,
  alter_fun = alter_fun, col = col,
  remove_empty_columns = TRUE,
  remove_empty_rows = TRUE,
  column_title = column_title,
  heatmap_legend_param = heatmap_legend_param
)

oncoPrint(
  mat,
  alter_fun = alter_fun,
  col = col,
  remove_empty_columns = TRUE,
  remove_empty_rows = TRUE,
  top_annotation = HeatmapAnnotation(
    cbar = anno_oncoprint_barplot(),
    foo1 = 1:172,
    bar1 = anno_points(1:172)
  ),
  left_annotation = rowAnnotation(foo2 = 1:26),
  right_annotation = rowAnnotation(bar2 = anno_barplot(1:26)),
  column_title = column_title,
  heatmap_legend_param = heatmap_legend_param
)

oncoPrint(
  mat,
  alter_fun = alter_fun,
  col = col,
  top_annotation = HeatmapAnnotation(
    column_barplot = anno_oncoprint_barplot(
      "HOMDEL",
      border = TRUE,
      show_fraction = TRUE,
      height = unit(4, "cm")
    )
  ),
  right_annotation = rowAnnotation(
    row_barplot = anno_oncoprint_barplot(
      c("AMP", "HOMDEL"),
      border = TRUE,
      height = unit(4, "cm"),
      axis_param = list(side = "bottom", labels_rot = 90)
    )
  ),
  remove_empty_columns = TRUE,
  remove_empty_rows = TRUE,
  column_title = column_title,
  heatmap_legend_param = heatmap_legend_param
)


sample_label <- read.table("../results/2191/NMF/cluster-nk-ligands-k=8_raw/sample_cluster_4.csv",
  header = F,
  sep = ",",
  check.names = FALSE
)

cluster1 <- sample_label[which(sample_label$V2 == 1), 1]
cluster1 <- intersect(cluster1, intersect_samples)

cluster2 <- sample_label[which(sample_label$V2 == 2), 1]
cluster2 <- intersect(cluster2, intersect_samples)

cluster3 <- sample_label[which(sample_label$V2 == 3), 1]
cluster3 <- intersect(cluster3, intersect_samples)

cluster4 <- sample_label[which(sample_label$V2 == 4), 1]
cluster4 <- intersect(cluster4, intersect_samples)

expmat1 <- raw_tcga[select_genes, cluster1]
expmat2 <- raw_tcga[select_genes, cluster2]
expmat3 <- raw_tcga[select_genes, cluster3]
expmat4 <- raw_tcga[select_genes, cluster4]
expmat5 <- raw_gtex_mrna[select_genes, ]
mat1 <- mat[, cluster1]
mat2 <- mat[, cluster2]
mat3 <- mat[, cluster3]
mat4 <- mat[, cluster4]

ht_list <- Heatmap(log10(expmat1 + 1),
  # matrix(rnorm(nrow(expmat) * ncol(expmat)), ncol = ncol(expmat)),
  name = "expr",
  show_column_names = F,
  # clustering_distance_rows = "spearman",
  col = circlize::colorRamp2(c(0, 3, 5), c("#2c58af", "white", "#d43a35")),
  width = unit(5, "cm")
) + Heatmap(log10(expmat2 + 1),
  # matrix(rnorm(nrow(expmat) * ncol(expmat)), ncol = ncol(expmat)),
  name = "expr",
  show_column_names = F,
  col = circlize::colorRamp2(c(0, 3, 5), c("#2c58af", "white", "#d43a35")),
  width = unit(5, "cm")
) + Heatmap(log10(expmat3 + 1),
  # matrix(rnorm(nrow(expmat) * ncol(expmat)), ncol = ncol(expmat)),
  name = "expr",
  show_column_names = F,
  col = circlize::colorRamp2(c(0, 3, 5), c("#2c58af", "white", "#d43a35")),
  width = unit(5, "cm")
) + Heatmap(log10(expmat4 + 1),
  # matrix(rnorm(nrow(expmat) * ncol(expmat)), ncol = ncol(expmat)),
  name = "expr",
  show_column_names = F,
  col = circlize::colorRamp2(c(0, 3, 5), c("#2c58af", "white", "#d43a35")),
  width = unit(5, "cm")
) + Heatmap(log10(expmat5 + 1),
  # matrix(rnorm(nrow(expmat) * ncol(expmat)), ncol = ncol(expmat)),
  name = "expr",
  show_column_names = F,
  col = circlize::colorRamp2(c(0, 3, 5), c("#2c58af", "white", "#d43a35")),
  width = unit(5, "cm")
) + oncoPrint(
  mat,
  alter_fun = alter_fun,
  col = col,
  column_title = column_title,
  heatmap_legend_param = heatmap_legend_param
)

draw(ht_list, n = 2)

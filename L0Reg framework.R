

library(L0Learn)
workplace <- paste(getwd(),'/results',sep = '')
setwd(workplace)
load('luad_mrna.Rdata')
sample_label <- read.table('sample_cluster_4.csv',
                           header = F,
                           sep=',',
                           check.names=FALSE)
for (j in 1:4) {
  if (dir.exists(sprintf('workplace/cluster%s',j)) == F) {
    dir.create(file.path(workplace, paste('cluster', j, sep = '')))
  }
  dir <- paste(workplace,'/cluster',j,sep = '')
  setwd(dir)
  genes_list <- read.table(sprintf('genes_list_cluster%s.txt',j),
                           header = T,
                           row.names = 1)
  cluster <- sample_label[which(sample_label$V2 == j), 1]
  tcga_cluster <- tcga_luad[, cluster]
  Y_raw <- tcga_cluster[which(row.names(tcga_cluster)%in%row.names(genes_list)), ]
  Y_raw <- as.data.frame(Y_raw[match(row.names(genes_list),row.names(Y_raw)),])
  maxSize = 100
  for(i in (1:nrow(Y_raw))){
    if (file.exists(paste(rownames(Y_raw[i,]), '_TFs_list.txt', sep = '')) == T) {
      Y <- t(tcga_cluster[which(row.names(tcga_cluster)%in%row.names(genes_list)), ][i, ])
      TFs_list  <- read.table(sprintf('%s_TFs_list.txt',rownames(Y_raw[i, ])),
                              header = T,
                              row.names = 1)
      row_matrix_TFs <- row.names(tcga_cluster)%in%row.names(TFs_list)
      X <- tcga_cluster[which(row_matrix_TFs), ]
      X_t <- t(X)
      gene_no = ncol(X_t)
      X = as.matrix(X_t[, 1:gene_no])
      row.names(X) = row.names(X_t)
      target_gene <- rownames(Y_raw[i, ])
      cvfit_L0 <- L0Learn.fit(X, Y,
                              penalty="L0",
                              maxSuppSize=100)
      lambda_value <- as.data.frame(cvfit_L0$lambda)
      names(lambda_value) <- 'lambda'
      lambda_value[nrow(lambda_value)/2,]
      y_cat = predict(cvfit_L0,
                      newx=X,
                      lambda=lambda_value[nrow(lambda_value)/2, ],
                      gamma=0)
      y_hat=as.vector(y_cat)
      plot(cvfit_L0, gamma = cvfit_L0$gamma, showLines=TRUE)
      FuData = coef(cvfit_L0, lambda=lambda_value[nrow(lambda_value)/2, ], gamma=cvfit_L0$gamma)
      FuData = as.vector(FuData)
      FuData = FuData[-1]
      FuData = which(FuData!=0)
      FuData = colnames(X)[FuData]
      if (length(FuData) == 1) {
        X_Y=cbind(X[, FuData], Y)
        colnames(X_Y)[1] <- FuData
        X_Y_frame = as.data.frame(X_Y)
      }else{
        X_Y=cbind(X[, FuData],Y)
        X_Y_frame = as.data.frame(X_Y)
      }
      lmfit = lm(Y~., data=X_Y_frame)
      fit_FuData=summary(lmfit)
      res_data <- as.matrix(fit_FuData$coefficients)
      res_data
      L0_Rsquare = 1-mse(Y,y_hat)/var(Y)
      evaluate_L0 <- data.frame(L0_Rsquare = 1-mse(Y,y_hat)/var(Y), L0_RMSE = RMSE(Y,y_hat))
      res_data_f <- as.data.frame(res_data)
      res_data_f <- res_data_f[which(res_data_f$`Pr(>|t|)` <= 0.05),]
      if (nrow(res_data_f)>0) {
        for (k in 1:nrow(res_data_f)) {
          if (res_data_f$Estimate[k] < 0) {
            res_data_f$reg[k] <- '2'
          }else
            res_data_f$reg[k] <- '1'
        }
        write.csv(res_data_f,sprintf('%s_TFs_Selected_L0.csv',target_gene))
        res_TRN <- cbind(row.names(res_data_f),res_data_f$reg, rownames(Y_raw[i,]), paste('cluster', j, sep = ''))
        res_TRN_cluster <- rbind(res_TRN_cluster, res_TRN)
      }
    }
  }
}

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    
    sample_cluster = pd.read_csv(os.path.join(path_read, "sample_cluster_4.csv"), header=None)
    
    genes_17_file = os.path.join(path_read, "Genes_17.txt")
    if os.path.exists(genes_17_file):
        gene_list = pd.read_table(genes_17_file, header=0).iloc[:, 0].astype(str).tolist()
    else:
        print("Warning: Genes_17.txt not found")
        return 1
    
    tcga_cnv = pd.read_csv(os.path.join(path_save, "CNV-LUAD.csv"), index_col=0)
    tcga_luad = pd.read_csv(os.path.join(path_save, "TCGA-LUAD.csv"), index_col=0)
    
    records = []
    for k in range(1, 5):
        cls_samples = sample_cluster.loc[sample_cluster.iloc[:, 1] == k, 0].astype(str).tolist()
        cls_samples = [s for s in cls_samples if s in tcga_cnv.columns and s in tcga_luad.columns]
        if len(cls_samples) == 0:
            continue
        
        cnv_data = tcga_cnv.reindex(index=gene_list, columns=cls_samples).T
        mrna_data = tcga_luad.reindex(index=gene_list, columns=cls_samples).T
        mrna_data = np.log2(mrna_data + 1.0)
        
        inter_cor_matrix = np.zeros((len(gene_list), len(gene_list)))
        for i in range(len(gene_list)):
            if gene_list[i] not in cnv_data.columns:
                continue
            for j in range(len(gene_list)):
                if gene_list[j] not in mrna_data.columns:
                    continue
                x = cnv_data[gene_list[i]].astype(float).to_numpy()
                y = mrna_data[gene_list[j]].astype(float).to_numpy()
                mask = np.isfinite(x) & np.isfinite(y)
                if mask.sum() < 3:
                    continue
                r, _ = spearmanr(x[mask], y[mask])
                inter_cor_matrix[i, j] = r
        
        for i in range(len(gene_list)):
            select_gene = gene_list[i]
            if select_gene not in cnv_data.columns or select_gene not in mrna_data.columns:
                continue
            correction_value = inter_cor_matrix[i, i]
            print(f"The correlation coefficient of {select_gene}: {correction_value}")
            records.append({"Gene": select_gene, "CNVs": correction_value, "Cluster": f"Cluster{k}"})
    
    out = pd.DataFrame.from_records(records, columns=["Gene", "CNVs", "Cluster"]).dropna()
    os.makedirs(path_save, exist_ok=True)
    out.to_csv(os.path.join(path_save, "CNV_correlation_value.csv"), index=False)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

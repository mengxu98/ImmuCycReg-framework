import os
import sys
import numpy as np
import pandas as pd
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr


def run_deseq2(tcga_cluster_df: pd.DataFrame, gtex_df: pd.DataFrame) -> pd.DataFrame:
    deseq2 = importr('DESeq2')
    base = importr('base')
    
    tcga_gtex = pd.concat([tcga_cluster_df, gtex_df], axis=1)
    group = ["C1"] * tcga_cluster_df.shape[1] + ["C2"] * gtex_df.shape[1]
    group_df = pd.DataFrame({"groupInfo": group}, index=tcga_gtex.columns)
    
    with localconverter(default_converter + pandas2ri.converter):
        r_counts = ro.conversion.py2rpy(tcga_gtex.round().astype(int))
        r_coldata = ro.conversion.py2rpy(group_df)
    
    dds = deseq2.DESeqDataSetFromMatrix(countData=r_counts, colData=r_coldata, design=ro.Formula("~ groupInfo"))
    dds = deseq2.DESeq(dds)
    res = deseq2.results(dds, contrast=ro.Vector(["groupInfo", "C1", "C2"]))
    
    with localconverter(default_converter + pandas2ri.converter):
        res_df = ro.conversion.rpy2py(res)
    
    res_df.index.name = 'gene'
    return res_df


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    os.makedirs(os.path.join(path_save, 'DESeq2'), exist_ok=True)
    os.makedirs(os.path.join(path_save, 'Figures'), exist_ok=True)
    
    genes_249_file = os.path.join(path_read, "Genes_249.txt")
    if os.path.exists(genes_249_file):
        genes_249 = pd.read_table(genes_249_file, header=0).iloc[:, 0].astype(str).tolist()
    else:
        print("Warning: Genes_249.txt not found, using all genes")
        genes_249 = None
    
    samples_cluster = pd.read_csv(os.path.join(path_read, "sample_cluster_4.csv"), header=None)
    
    tcga = pd.read_csv(os.path.join(path_save, "TCGA-LUAD.csv"), index_col=0)
    gtex = pd.read_csv(os.path.join(path_save, "GTEx-LUAD.csv"), index_col=0)
    
    if genes_249 is not None:
        tcga = tcga.reindex(index=genes_249)
        gtex = gtex.reindex(index=genes_249)
    
    VolcanoPlotList = []
    
    for j in range(1, 5):
        print(f"Cluster {j} DESeq2 start......")
        cls_samples = samples_cluster.loc[samples_cluster.iloc[:, 1] == j, 0].astype(str).tolist()
        cls_samples = [s for s in cls_samples if s in tcga.columns]
        if len(cls_samples) == 0:
            print(f"Warning: No samples found for cluster {j}")
            continue
        
        tcga_cluster = tcga.loc[:, cls_samples].copy().dropna()
        gtex_cluster = gtex.copy().dropna()
        
        common_genes = list(set(tcga_cluster.index) & set(gtex_cluster.index))
        if len(common_genes) == 0:
            print(f"Warning: No common genes for cluster {j}")
            continue
        
        tcga_cluster = tcga_cluster.loc[common_genes, :]
        gtex_cluster = gtex_cluster.loc[common_genes, :]
        
        res_df = run_deseq2(tcga_cluster, gtex_cluster)
        res_df = res_df.sort_values(by='padj', na_last=True)
        
        padj_threshold = 0.05
        fold_change = 1.0
        res_df['reg'] = 'not DE'
        res_df.loc[(res_df['log2FoldChange'] >= 0.5), 'reg'] = 'min-up-regulated'
        res_df.loc[(res_df['log2FoldChange'] <= -0.5), 'reg'] = 'min-down-regulated'
        res_df.loc[(res_df['log2FoldChange'] >= fold_change) & (res_df['padj'] <= padj_threshold), 'reg'] = 'up-regulated'
        res_df.loc[(res_df['log2FoldChange'] <= -fold_change) & (res_df['padj'] <= padj_threshold), 'reg'] = 'down-regulated'
        
        if genes_249 is not None:
            res_df = res_df.reindex(index=genes_249)
        
        res_df = res_df.reset_index()
        res_df = res_df.rename(columns={'index': 'gene'})
        
        col_order = ['gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'reg']
        res_df = res_df[[col for col in col_order if col in res_df.columns]]
        
        res_df.to_csv(os.path.join(path_save, f"DESeq2_results_cluster{j}.csv"), index=False)
        print(f"Cluster {j} DESeq2 done......")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

import os
import sys
import numpy as np
import pandas as pd
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr


def main():
    ccp = importr('ConsensusClusterPlus')
    stats = importr('stats')
    
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    os.makedirs(os.path.join(path_save, 'NMF'), exist_ok=True)
    
    tcga_luad = pd.read_csv(os.path.join(path_save, 'TCGA-LUAD.csv'), index_col=0)
    
    genes2230_file = os.path.join(path_read, 'genes2230.csv')
    if os.path.exists(genes2230_file):
        genes2230 = pd.read_csv(genes2230_file)
        gene_col = 'Gene' if 'Gene' in genes2230.columns else genes2230.columns[0]
        genes_list = genes2230[gene_col].astype(str).tolist()
        tcga_luad = tcga_luad.reindex(index=genes_list)
    
    data_nmf = np.log10(tcga_luad.to_numpy().astype(float) + 1.0)
    mads = np.nanmedian(np.abs(data_nmf - np.nanmedian(data_nmf, axis=1, keepdims=True)), axis=1)
    order = np.argsort(mads)[::-1]
    dataset = data_nmf[order[:30], :]
    dataset_df = pd.DataFrame(dataset, index=tcga_luad.index[order[:30]], columns=tcga_luad.columns)
    
    with localconverter(default_converter + pandas2ri.converter):
        r_dataset = ro.conversion.py2rpy(dataset_df)
    
    title = os.path.join(path_save, 'NMF') + '/'
    os.makedirs(title, exist_ok=True)
    
    results = ccp.ConsensusClusterPlus(r_dataset,
                                       maxK=6,
                                       reps=50,
                                       pItem=0.8,
                                       pFeature=0.8,
                                       clusterAlg="pam",
                                       distance="euclidean",
                                       seed=ro.StrVector(["20210101"]),
                                       title=title,
                                       corUse="complete.obs",
                                       plot="pdf")
    
    for k in range(2, 7):
        cls = results.rx2(k).rx2("consensusClass")
        cls_np = np.array(cls, dtype=int)
        cls_df = pd.DataFrame({'cluster': cls_np}, index=dataset_df.columns)
        cls_df.to_csv(os.path.join(title, f"sample_cluster_{k}.csv"), header=False, index=True)
    
    Kvec = list(range(2, 7))
    x1 = 0.1
    x2 = 0.9
    PAC = {}
    
    for i in Kvec:
        M = np.array(results.rx2(i).rx2("consensusMatrix"))
        lower_tri = np.tril(M, k=-1)
        M_vals = lower_tri[lower_tri != 0]
        Fn_x1 = np.mean(M_vals <= x1)
        Fn_x2 = np.mean(M_vals <= x2)
        PAC[f"K={i}"] = Fn_x2 - Fn_x1
    
    optK = min(Kvec, key=lambda k: PAC[f"K={k}"])
    print(f"Optimal K: {optK}")
    print(f"PAC values: {PAC}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

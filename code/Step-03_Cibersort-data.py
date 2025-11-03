import os
import sys
import pandas as pd


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    os.makedirs(path_save, exist_ok=True)
    
    tcga_luad = pd.read_csv(os.path.join(path_save, "TCGA-LUAD-tpm-tcga-t_unnormlized.csv"), index_col=0)
    gtex_luad = pd.read_csv(os.path.join(path_save, "GTEx-LUAD-tpm-gtex_unnormlized.csv"), index_col=0)
    
    sample_label = pd.read_csv(os.path.join(path_read, "sample_cluster_4.csv"), 
                               header=None, names=["sample", "cluster"])
    
    data_list = []
    for i in range(1, 5):
        cluster_samples = sample_label.loc[sample_label["cluster"] == i, "sample"].astype(str).tolist()
        cluster_samples = [s for s in cluster_samples if s in tcga_luad.columns]
        if len(cluster_samples) > 0:
            data_list.append(tcga_luad.loc[:, cluster_samples])
        else:
            data_list.append(pd.DataFrame())
    
    data_list.append(gtex_luad)
    
    cibersort_data = pd.concat(data_list, axis=1)
    cibersort_data.to_csv(os.path.join(path_save, "Cibersort_data_ALL.txt"), sep="\t", index=True)
    
    cibersort_data_label_list = []
    for i in range(5):
        if len(data_list[i]) == 0:
            continue
        cluster = data_list[i].T
        if i < 4:
            cluster_label = f"Cluster{i+1}"
        else:
            cluster_label = "GTEx"
        cluster['group'] = cluster_label
        cibersort_data_label_list.append(cluster.T)
    
    cibersort_data_label = pd.concat(cibersort_data_label_list, axis=1)
    cibersort_data_label.to_csv(os.path.join(path_save, "Data_ALL_group.txt"), sep="\t", index=True)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())


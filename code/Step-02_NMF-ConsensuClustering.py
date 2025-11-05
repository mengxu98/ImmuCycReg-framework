import os
import sys
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import random


def pam_clustering(data, n_clusters, distance='euclidean', max_iter=100):
    """PAM (Partitioning Around Medoids) clustering"""
    n_samples = data.shape[0]
    if n_samples < n_clusters:
        n_clusters = n_samples
    
    if distance == 'precomputed':
        dist_matrix = data
    else:
        dist_matrix = pairwise_distances(data, metric=distance)
    
    medoids = np.random.choice(n_samples, size=n_clusters, replace=False)
    
    for iteration in range(max_iter):
        clusters = np.argmin(dist_matrix[:, medoids], axis=1)
        
        new_medoids = np.zeros(n_clusters, dtype=int)
        for k in range(n_clusters):
            cluster_mask = clusters == k
            if cluster_mask.sum() == 0:
                new_medoids[k] = medoids[k]
            else:
                cluster_dist = dist_matrix[np.ix_(cluster_mask, cluster_mask)]
                new_medoids[k] = np.where(cluster_mask)[0][np.argmin(cluster_dist.sum(axis=0))]
        
        if np.array_equal(np.sort(medoids), np.sort(new_medoids)):
            break
        medoids = new_medoids
    
    return clusters


def consensus_cluster_plus(data, maxK=6, reps=50, pItem=0.8, pFeature=0.8, 
                          clusterAlg='pam', distance='euclidean', seed=None):
    """Python implementation of ConsensusClusterPlus"""
    if seed is not None:
        np.random.seed(int(seed))
        random.seed(int(seed))
    
    n_samples, n_features = data.shape
    
    results = {}
    
    for k in range(2, maxK + 1):
        print(f"Clustering for K={k}...")
        
        consensus_matrix = np.zeros((n_samples, n_samples))
        count_matrix = np.zeros((n_samples, n_samples))
        
        for rep in range(reps):
            n_items = max(1, int(n_samples * pItem))
            n_feat = max(1, int(n_features * pFeature))
            
            item_indices = np.random.choice(n_samples, size=n_items, replace=False)
            feature_indices = np.random.choice(n_features, size=n_feat, replace=False)
            
            sub_data = data[np.ix_(item_indices, feature_indices)]
            
            if clusterAlg == 'pam':
                clusters = pam_clustering(sub_data, k, distance=distance)
            elif clusterAlg == 'hc':
                dists = pairwise_distances(sub_data, metric=distance)
                linkage_matrix = linkage(squareform(dists, checks=False), method='average')
                clusters = fcluster(linkage_matrix, k, criterion='maxclust') - 1
            else:
                kmeans = KMeans(n_clusters=k, random_state=seed, n_init=10)
                clusters = kmeans.fit_predict(sub_data)
            
            for i, idx_i in enumerate(item_indices):
                for j, idx_j in enumerate(item_indices):
                    if clusters[i] == clusters[j]:
                        consensus_matrix[idx_i, idx_j] += 1
                    count_matrix[idx_i, idx_j] += 1
        
        consensus_matrix = consensus_matrix / np.maximum(count_matrix, 1)
        np.fill_diagonal(consensus_matrix, 1.0)
        
        if clusterAlg == 'hc':
            dists = 1 - consensus_matrix
            np.fill_diagonal(dists, 0)
            linkage_matrix = linkage(squareform(dists, checks=False), method='average')
            consensus_class = fcluster(linkage_matrix, k, criterion='maxclust')
        elif clusterAlg == 'pam':
            dists = 1 - consensus_matrix
            np.fill_diagonal(dists, 0)
            consensus_class = pam_clustering(dists, k, distance='precomputed')
            consensus_class = consensus_class + 1
        else:
            kmeans = KMeans(n_clusters=k, random_state=seed, n_init=10)
            consensus_class = kmeans.fit_predict(consensus_matrix) + 1
        
        results[k] = {
            'consensusMatrix': consensus_matrix,
            'consensusClass': consensus_class
        }
    
    return results


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    os.makedirs(os.path.join(path_save, 'NMF'), exist_ok=True)
    
    tcga_luad = pd.read_csv(os.path.join(path_save, 'TCGA-LUAD.csv'), index_col=0)
    
    genes2230_file = os.path.join(path_read, 'data/Genes_2230.csv')
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
    
    title = os.path.join(path_save, 'NMF') + '/'
    os.makedirs(title, exist_ok=True)
    
    results = consensus_cluster_plus(
        dataset_df.T.values,
        maxK=6,
        reps=50,
        pItem=0.8,
        pFeature=0.8,
        clusterAlg='pam',
        distance='euclidean',
        seed=20210101
    )
    
    for k in range(2, 7):
        cls = results[k]['consensusClass']
        cls_df = pd.DataFrame({'cluster': cls}, index=dataset_df.columns)
        cls_df.to_csv(os.path.join(title, f"sample_cluster_{k}.csv"), header=False, index=True)
    
    Kvec = list(range(2, 7))
    x1 = 0.1
    x2 = 0.9
    PAC = {}
    
    for i in Kvec:
        M = results[i]['consensusMatrix']
        lower_tri = np.tril(M, k=-1)
        M_vals = lower_tri[lower_tri != 0]
        if len(M_vals) > 0:
            Fn_x1 = np.mean(M_vals <= x1)
            Fn_x2 = np.mean(M_vals <= x2)
            PAC[f"K={i}"] = Fn_x2 - Fn_x1
        else:
            PAC[f"K={i}"] = np.nan
    
    optK = min(Kvec, key=lambda k: PAC.get(f"K={k}", np.inf))
    print(f"Optimal K: {optK}")
    print(f"PAC values: {PAC}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

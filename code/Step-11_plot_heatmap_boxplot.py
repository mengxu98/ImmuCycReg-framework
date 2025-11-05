import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches


def check_dir(dirpath):
    os.makedirs(dirpath, exist_ok=True)
    return dirpath


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    os.makedirs(os.path.join(path_save, "Figures", "boxplot"), exist_ok=True)
    
    tcga_luad = pd.read_csv(os.path.join(path_save, "TCGA-LUAD.csv"), index_col=0)
    gtex_luad = pd.read_csv(os.path.join(path_save, "GTEx-LUAD.csv"), index_col=0)
    
    samples_20_file = os.path.join(path_read, "Samples_20.csv")
    if os.path.exists(samples_20_file):
        samples_20 = pd.read_csv(samples_20_file).iloc[:, 0].tolist()
    else:
        samples_20 = tcga_luad.columns[:20].tolist()
    
    genes_17_file = os.path.join(path_read, "Genes_17.txt")
    if os.path.exists(genes_17_file):
        gene_list = pd.read_table(genes_17_file, header=0).iloc[:, 0].astype(str).tolist()
    else:
        print("Error: Genes_17.txt not found")
        return 1
    
    sample_label = pd.read_csv(os.path.join(path_read, "sample_cluster_4.csv"),
                               header=None, names=["sample", "cluster"])
    
    tcga = tcga_luad.reindex(index=gene_list).dropna()
    gtex = gtex_luad.reindex(index=gene_list).dropna()
    
    tcga_scaled = tcga.T.apply(lambda x: (x - x.mean()) / x.std(), axis=0).T
    tcga_scaled = tcga_scaled.clip(-2, 2)
    
    gtex_scaled = gtex.T.apply(lambda x: (x - x.mean()) / x.std(), axis=0).T
    gtex_scaled = gtex_scaled.clip(-2, 2)
    
    samples_20_in_tcga = [s for s in samples_20 if s in tcga_scaled.columns]
    if len(samples_20_in_tcga) > 0:
        fig, ax = plt.subplots(figsize=(10, 4))
        sns.heatmap(tcga_scaled.loc[:, samples_20_in_tcga], cmap='RdBu_r', center=0,
                   vmin=-2, vmax=2, cbar_kws={'label': 'Z-score'},
                   xticklabels=False, yticklabels=True, ax=ax)
        ax.set_title('TCGA Expression Heatmap')
        plt.tight_layout()
        plt.savefig(os.path.join(path_save, "Figures", "TCGA_heatmap_sample20.pdf"))
        plt.close()
    
    exp_mat_list = []
    for i in range(1, 5):
        cluster_samples = sample_label[sample_label['cluster'] == i]['sample'].tolist()
        cluster_samples = [s for s in cluster_samples if s in tcga_scaled.columns]
        if len(cluster_samples) > 0:
            exp_mat_list.append(tcga_scaled.loc[:, cluster_samples])
        else:
            exp_mat_list.append(pd.DataFrame())
    
    exp_mat_list.append(gtex_scaled)
    
    sample_colors = ["#2874c5", "#008a00", "#c6524a", "#eabf00", "#696969"]
    
    fig, axes = plt.subplots(1, len(exp_mat_list), figsize=(14, 4.5), gridspec_kw={'width_ratios': [len(df.columns) if len(df) > 0 else 1 for df in exp_mat_list]})
    
    if len(exp_mat_list) == 1:
        axes = [axes]
    
    for i, exp_mat in enumerate(exp_mat_list):
        if len(exp_mat) == 0:
            continue
        
        if i < 4:
            label = f"TCGA-Cluster{i+1}"
        else:
            label = "GTEx"
        
        sample_color = sample_colors[i]
        
        sns.heatmap(exp_mat, cmap='RdBu_r', center=0, vmin=-2, vmax=2,
                   cbar_kws={'label': 'Z-score'}, xticklabels=False, yticklabels=True,
                   ax=axes[i])
        
        axes[i].set_title(label)
        
        rect = Rectangle((0, 0), len(exp_mat.columns), 0.02 * len(exp_mat.index),
                        transform=axes[i].transAxes, facecolor=sample_color, clip_on=False)
        axes[i].add_patch(rect)
    
    plt.tight_layout()
    plt.savefig(os.path.join(path_save, "SuppFig5.pdf"), dpi=300)
    plt.close()
    
    for target_gene in gene_list:
        if target_gene not in tcga_scaled.index or target_gene not in gtex_scaled.index:
            continue
        
        data_tcga = pd.DataFrame({
            'Group': ['TCGA'] * len(tcga_scaled.columns),
            'Expression': tcga_scaled.loc[target_gene, :].values
        })
        
        data_gtex = pd.DataFrame({
            'Group': ['GTEx'] * len(gtex_scaled.columns),
            'Expression': gtex_scaled.loc[target_gene, :].values
        })
        
        data_plot = pd.concat([data_tcga, data_gtex], ignore_index=True)
        
        fig, ax = plt.subplots(figsize=(2.5, 4))
        
        positions = [0, 1]
        groups = ['TCGA', 'GTEx']
        data_groups = [data_plot[data_plot['Group'] == g]['Expression'].values for g in groups]
        
        bp = ax.boxplot(data_groups, positions=positions, patch_artist=True,
                       labels=groups, widths=0.6)
        
        colors = ["#e1181f", "#253870"]
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax.set_ylabel('Expression')
        ax.set_title(target_gene)
        ax.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        
        plt.savefig(os.path.join(path_save, "Figures", "boxplot", f"{target_gene}.pdf"),
                   bbox_inches='tight')
        plt.close()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

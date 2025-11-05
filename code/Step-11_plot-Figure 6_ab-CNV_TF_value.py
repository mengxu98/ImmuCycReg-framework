import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    os.makedirs(os.path.join(path_save, "Figures"), exist_ok=True)
    
    evaluate_r2_file = os.path.join(path_save, "TF_regulatory_value.csv")
    if not os.path.exists(evaluate_r2_file):
        print("Warning: TF_regulatory_value.csv not found")
        return 1
    
    evaluate_r2 = pd.read_csv(evaluate_r2_file)
    evaluate_r2 = evaluate_r2[['Gene', 'Cluster', 'L0_Rsquare']].copy()
    
    evaluate_cnvs_file = os.path.join(path_save, "CNV_correlation_value.csv")
    if not os.path.exists(evaluate_cnvs_file):
        print("Warning: CNV_correlation_value.csv not found")
        return 1
    
    evaluate_cnvs = pd.read_csv(evaluate_cnvs_file)
    evaluate_cnvs = evaluate_cnvs[['Gene', 'Cluster', 'CNVs']].copy()
    
    cluster_color = ["#2874c5", "#008a00", "#c6524a", "#eabf00"]
    
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    
    sns.barplot(data=evaluate_cnvs, x='Gene', y='CNVs', hue='Cluster', 
               palette=cluster_color, ax=axes[0])
    axes[0].set_ylabel('CNV correlation value')
    axes[0].set_ylim(-0.5, 0.5)
    axes[0].set_xlabel('Gene')
    axes[0].tick_params(axis='x', rotation=45)
    axes[0].invert_yaxis()
    
    sns.barplot(data=evaluate_r2, x='Gene', y='L0_Rsquare', hue='Cluster',
               palette=cluster_color, ax=axes[1])
    axes[1].set_ylabel('TF regulatory value')
    axes[1].set_ylim(0, 1)
    axes[1].set_xlabel('')
    axes[1].set_xticks([])
    axes[1].set_xticklabels([])
    axes[1].invert_yaxis()
    
    axes[0].legend(loc='lower left')
    axes[1].legend(loc='lower left')
    
    plt.tight_layout()
    plt.savefig(os.path.join(path_save, "Figures", "Figure 6_ab-CNV_TF_value.pdf"),
               bbox_inches='tight')
    plt.close()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
